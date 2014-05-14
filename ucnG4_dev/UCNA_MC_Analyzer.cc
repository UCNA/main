#include <TMath.h>
#include <cfloat>
#include <string.h>
#include <algorithm>

#include "UCNA_MC_Analyzer.hh"

UCNA_MC_Analyzer::UCNA_MC_Analyzer(const std::string& outfname): ucnG4_analyzer(outfname), saveAllEvents(false), undeadLayer(false), calcCathCharge(false) {
	// set up cathode wire positions array, 0.1 inch wire spacing in cm
	for(int i = 0; i<kWiresPerCathode*kMaxCathodes; i++)
		cathWirePos[i] = (i-0.5*(kWiresPerCathode*kMaxCathodes-1.))*(0.1*2.54);
	// set up circle offsets
	for(int i=0; i<N_CHG_CIRC_PTS; i++) {
		double th = 2*M_PI*i/float(N_CHG_CIRC_PTS);
		circ_pts[X_DIRECTION][i] = cos(th);
		circ_pts[Y_DIRECTION][i] = sin(th);
	}
}

void UCNA_MC_Analyzer::setupOutputTree() {
	printf("Adding branches for UCNA_MC_Analyzer...\n");
	anaTree->Branch("Edep",DE0.Edep,"EdepE/D:EdepW/D");
	anaTree->Branch("EdepQ",DE0.EdepQ,"EdepQE/D:EdepQW/D");
	anaTree->Branch("EdepAll",&DE0.EdepAll,"EdepAll/D");
	anaTree->Branch("MWPCEnergy",DE0.fMWPCEnergy,"MWPCEnergyE/D:MWPCEnergyW/D");
	anaTree->Branch("time",DE0.hitTime,"timeE/D:timeW/D");
	anaTree->Branch("tmTime",DE0.trapMonTime,"tmTimeE/D:tmTimeW/D");
	anaTree->Branch("subEvt",&DE0.subEvt,"subEvt/I");
	
	char tmp[1024];
	sprintf(tmp,"EdepSD[%i]/D",N_SD);
	anaTree->Branch("EdepSD",DE0.EdepSD,tmp);
	sprintf(tmp,"thetaInSD[%i]/D",N_SD);
	anaTree->Branch("thetaInSD",DE0.thetaInSD,tmp);
	sprintf(tmp,"thetaOutSD[%i]/D",N_SD);
	anaTree->Branch("thetaOutSD",DE0.thetaOutSD,tmp);
	sprintf(tmp,"keInSD[%i]/D",N_SD);
	anaTree->Branch("keInSD",DE0.keInSD,tmp);
	sprintf(tmp,"keOutSD[%i]/D",N_SD);
	anaTree->Branch("keOutSD",DE0.keOutSD,tmp);
	sprintf(tmp,"hitCountSD[%i]/I",N_SD);
	anaTree->Branch("hitCountSD",DE0.hitCountSD,tmp);
	
	anaTree->Branch("MWPCPos",DE0.MWPCPos,"MWPCPosE[3]/D:MWPCPosW[3]/D");
	anaTree->Branch("ScintPos",DE0.ScintPos,"ScintPosE[3]/D:ScintPosW[3]/D");
	anaTree->Branch("MWPCPosSigma",DE0.MWPCPosSigma,"MWPCPosSigmaE[3]/D:MWPCPosSigmaW[3]/D");
	anaTree->Branch("ScintPosSigma",DE0.ScintPosSigma,"ScintPosSigmaE[3]/D:ScintPosSigmaW[3]/D");
	
	if(calcCathCharge) {
		for(Side sd = EAST; sd <= WEST; ++sd) {
			for(AxisDirection d=X_DIRECTION; d<=Y_DIRECTION; ++d) {
				std::string cname = sideSubst("Cath_%c",sd)+(d==X_DIRECTION?"X":"Y");
				anaTree->Branch(cname.c_str(),DE0.cathCharge[sd][d],(cname+"["+itos(kMaxCathodes)+"]/F").c_str());
			}
		}
	}
}

void ucnaDAQEvt::reset(bool calcCathCharge) {
	subEvt = 0;
	for(Side sd = EAST; sd <= WEST; ++sd) {
		Edep[sd] = EdepQ[sd] = 0;
		fMWPCEnergy[sd]=0;
		hitTime[sd]=trapMonTime[sd]=FLT_MAX;
		for(AxisDirection d=X_DIRECTION; d<=Z_DIRECTION; ++d)
			MWPCPos[sd][d] = ScintPos[sd][d] = MWPCPosSigma[sd][d] = ScintPosSigma[sd][d] = 0;
		if(calcCathCharge)
			for(AxisDirection d=X_DIRECTION; d<=Y_DIRECTION; ++d)
				memset(cathCharge[sd][d],0,sizeof(cathCharge[sd][d][0])*kMaxCathodes);
	}
	EdepAll = 0;
	for(size_t ii=0; ii<N_SD; ii++) {
		hitCountSD[ii] = keInSD[ii] = keOutSD[ii] = thetaInSD[ii] = thetaOutSD[ii] = EdepSD[ii] = 0.;
		hitTimeSD[ii]=FLT_MAX;
	}
}

void UCNA_MC_Analyzer::resetAnaEvt() {
	DE0.reset(calcCathCharge);
	subEvts.clear();
	scintHitTimes.clear();
	triggerTimes.clear();
	pass_number = 2;
}

void UCNA_MC_Analyzer::processTrack() {
	
	// detector ID numbers
	const int ID_scint[2] = {0,10};
	const int ID_dscint[2] = {1,11};
	const int ID_mwpc[2] = {6,16};
	const int ID_deadmwpc[2] = {8,18};
	const int ID_winOut[2] = {3,13};
	const int ID_winIn[2] = {4,14};
	
	const int ID_trapmon[2] = {21,22};
	// particle ID numbers
	const int PDG_ELECTRON = 11;
	const int PDG_POSITRON = -11;
	
	const double kMe = 510.998910; // electron mass, keV
	
	// first pass to collect timing structure
	if(pass_number == 2) {
		for(Side sd = EAST; sd <= WEST; ++sd)
			if(detectorID==ID_scint[sd] || (undeadLayer && detectorID==ID_dscint[sd]))
				if(trackinfo->Edep > 1.0) scintHitTimes.push_back(trackinfo->hitTime * 1.e-9);
		return;
	}
	
	// select appropriate event to fill
	ucnaDAQEvt* DE = &DE0;
	if(triggerTimes.size()>2) {
		// determine nearest trigger
		double t = trackinfo->hitTime * 1e-9;
		std::vector<double>::iterator it = std::lower_bound(triggerTimes.begin(), triggerTimes.end(), t); // first element >= t
		assert(it != triggerTimes.end());
		unsigned int i = it-triggerTimes.begin();
		assert(i>0);
		if(*it - t > t - triggerTimes[i-1]) i--;
		double t0 = triggerTimes[i];
		
		if(i>0 && i < triggerTimes.size()-1) DE = &subEvts[i-1];	// by default, assign deposition to nearest-time trigger
		if( detectorID==ID_scint[EAST] || detectorID==ID_scint[WEST] || (undeadLayer && (detectorID==ID_dscint[EAST] || detectorID==ID_dscint[WEST])) ) {
			// PMT QADC gate, 140ns (with 5ns slop for prior events)
			if(!(t0-5e-9 < t && t < t0 + 140.e-9)) DE = &DE0;
		}
		if(detectorID==ID_mwpc[EAST] || detectorID==ID_mwpc[WEST] || detectorID==ID_deadmwpc[EAST] || detectorID==ID_deadmwpc[WEST]) {
			// MWPC PADC gate, 6us (with 1us slop for prior events)
			if(!(t0-1e-6 < t && t < t0 + 6e-6)) DE = &DE0;
		}
	}
	if(!saveAllEvents && DE==&DE0) return;
	
	// total deposited energy in all sensitive volumes
	DE->EdepAll += trackinfo->Edep;
	if(detectorID < N_SD) {
		DE->EdepSD[detectorID] += trackinfo->Edep;
		if(trackinfo->isEntering && pID==PDG_ELECTRON)
			DE->hitCountSD[detectorID]++;
		if(pID==PDG_ELECTRON && trackinfo->hitTime < DE->hitTimeSD[detectorID]) {
			DE->hitTimeSD[detectorID] = trackinfo->hitTime;
			double pin_x = trackinfo->pIn[0];
			double pin_y = trackinfo->pIn[1];
			double pin_z = trackinfo->pIn[2];
			double pout_x = trackinfo->pOut[0];
			double pout_y = trackinfo->pOut[1];
			double pout_z = trackinfo->pOut[2];
			if(pin_x || pin_y || pin_z) {
				double magpin2 = pin_x*pin_x+pin_y*pin_y+pin_z*pin_z;
				DE->keInSD[detectorID] = sqrt(magpin2+kMe*kMe)-kMe;
				DE->thetaInSD[detectorID] = acos(fabs(pin_z/sqrt(magpin2)));
			}
			if(pout_x || pout_y || pout_z) {
				double magpout2 = pout_x*pout_x+pout_y*pout_y+pout_z*pout_z;
				DE->keOutSD[detectorID] = sqrt(magpout2+kMe*kMe)-kMe;
				DE->thetaOutSD[detectorID] = acos(fabs(pout_z/sqrt(magpout2)));
			}
		}
	}
	
	for(Side sd = EAST; sd <= WEST; ++sd) {
		
		// scintillator deposited energy, position, hit time
		if(detectorID==ID_scint[sd] || (undeadLayer && detectorID==ID_dscint[sd])) {
			DE->Edep[sd] += trackinfo->Edep;
			DE->EdepQ[sd] += trackinfo->EdepQuenched;
			for(AxisDirection d=X_DIRECTION; d<=Z_DIRECTION; ++d) {
				DE->ScintPos[sd][d] += trackinfo->edepPos[d];
				DE->ScintPosSigma[sd][d] += trackinfo->edepPos2[d];
			}
			//earliest hit time with possibly detectable (not really) Edep; hitTime[EAST,WEST] have been initialized to FLT_MAX
			if(trackinfo->Edep>5.0 && trackinfo->hitTime < DE->hitTime[sd])
				DE->hitTime[sd] = trackinfo->hitTime;
		}
		
		// trap monitor timing for electrons, initialized to FLT_MAX
		if(detectorID==ID_trapmon[sd] && (pID==PDG_ELECTRON || pID==PDG_POSITRON)){
			if(trackinfo->hitTime < DE->trapMonTime[sd])
				DE->trapMonTime[sd] = trackinfo->hitTime; //earliest hit time
		}
		
		// wirechamber deposited energy and position
		double anodeScale = 0;
		if(detectorID==ID_deadmwpc[sd]) // contribution from "dead" gas
			anodeScale = 0.5;
		if(detectorID==ID_mwpc[sd])
			anodeScale = 1.0;
		if(anodeScale && trackinfo->Edep) {
			Double_t Ew = anodeScale*trackinfo->Edep;
			
			DE->fMWPCEnergy[sd] += Ew;
			for(AxisDirection d=X_DIRECTION; d<=Z_DIRECTION; ++d) {
				DE->MWPCPos[sd][d] += anodeScale*trackinfo->edepPos[d];
				DE->MWPCPosSigma[sd][d] += anodeScale*trackinfo->edepPos2[d];
			}
			if(calcCathCharge) {
				double c[Y_DIRECTION+1]; // track center
				for(AxisDirection d=X_DIRECTION; d<=Y_DIRECTION; ++d) c[d] = trackinfo->edepPos[d]/trackinfo->Edep;
				// charge radius
				Double_t cr = sqrt((trackinfo->edepPos2[X_DIRECTION]+trackinfo->edepPos2[Y_DIRECTION])/trackinfo->Edep
								   -c[X_DIRECTION]*c[X_DIRECTION]-c[Y_DIRECTION]*c[Y_DIRECTION]);
				if(!(cr==cr)) cr=0; // avoid NAN
				
				const double dx = 0.254;		// distance between wires, cm
				const double h = 1.0;			// plane spacing, cm
				const double dlambda = dx/h;	// wire spacing in plane spacing units
				
				// Matthieson 5.41+
				static const double K3 = 0.18; // read off chart for h/sd=3.93, r_a/sd=0.002
				static const double K2 = M_PI/2*(1-0.5*sqrt(K3));
				static const double K1 = K2*sqrt(K3)/(4.*atan(sqrt(K3))*N_CHG_CIRC_PTS);
				
				for(AxisDirection d=X_DIRECTION; d<=Y_DIRECTION; ++d) {
					for(int i = 0; i<kWiresPerCathode*kMaxCathodes; i++) {
						for(int j=0; j<N_CHG_CIRC_PTS; j++) {
							int flip = (sd==EAST && d==X_DIRECTION)?-1:1; // flip for East X because of mirrored position coordinates
							const double lambda = (c[d]+cr*circ_pts[d][j]-cathWirePos[i]*flip)/h;
							const double tkl = tanh(K2*lambda);
							DE->cathCharge[sd][d][i/kWiresPerCathode] += Ew * dlambda * K1 * (1-tkl*tkl) / (1+K3*tkl*tkl);
						}
					}
				}
			}
		}
	}
}

void ucnaDAQEvt::calcPositions() {
	for(Side sd = EAST; sd <= WEST; ++sd) {
		for(AxisDirection d=X_DIRECTION; d<=Z_DIRECTION; ++d) {
			if(Edep[sd]>0) {
				ScintPos[sd][d] /= Edep[sd];
				ScintPosSigma[sd][d] = sqrt(ScintPosSigma[sd][d]/Edep[sd]-ScintPos[sd][d]*ScintPos[sd][d]);
			}
			if(fMWPCEnergy[sd] > 0) {
				MWPCPos[sd][d] /= fMWPCEnergy[sd];
				MWPCPosSigma[sd][d] = sqrt(MWPCPosSigma[sd][d]/fMWPCEnergy[sd]-MWPCPos[sd][d]*MWPCPos[sd][d]);
			}
		}
	}
}

void UCNA_MC_Analyzer::processEvent() {
	
	if(pass_number==2) {
		// determine trigger timing structure
		const double busy_blackout = 12.e-6;	// DAQ "busy logic" blackout time [s]
		std::sort(scintHitTimes.begin(),scintHitTimes.end());
		triggerTimes.push_back(-1e6);
		for(unsigned int i=0; i<scintHitTimes.size(); i++) {
			if(scintHitTimes[i] > triggerTimes.back() + busy_blackout) {
				triggerTimes.push_back(scintHitTimes[i]);
				subEvts.push_back(DE0);
			}
		}
		triggerTimes.push_back(1e6);
		return;
	}
}

void UCNA_MC_Analyzer::writeEvent() {
	if(saveAllEvents) {
		DE0.calcPositions();
		anaTree->Fill();
	}
	for(unsigned int i=0; i<subEvts.size(); i++) {
		DE0 = subEvts[i];
		if(saveAllEvents || DE0.Edep[EAST]+DE0.Edep[WEST] > 0) {
			DE0.calcPositions();
			DE0.subEvt = i+1;
			anaTree->Fill();
		}
	}
}

int main(int argc, char** argv) {
	
	if(argc<3) {
		cout<<"Syntax: "<<argv[0]<<" <filename containing list of raw root files> <output root file name> [saveall|undead]"<<endl;
		exit(1);
	}
	
	UCNA_MC_Analyzer UMA(argv[2]);
	
	for(int i=3; i<argc; i++) {
		std::string arg = argv[i];
		if(arg=="saveall") {
			UMA.saveAllEvents = true;
		} else if(arg=="undead") {
			UMA.undeadLayer = true;
		} else if(arg=="cathodes") {
			UMA.calcCathCharge = true;
		} else {
			cout<<"Unknown argument: "<<argv[0]<<endl;
			exit(1);
		}
	}
	
	UMA.analyzeFileList(argv[1]);
	
	return 0;
}
