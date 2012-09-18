#include "UCNA_MC_Analyzer.hh"
#include <TMath.h>
#include <cfloat>

void UCNA_MC_Analyzer::setupOutputTree() {
	printf("Adding branches for UCNA_MC_Analyzer...\n");
	anaTree->Branch("Edep",Edep,"EdepE/D:EdepW/D");
	anaTree->Branch("EdepQ",EdepQ,"EdepQE/D:EdepQW/D");
	anaTree->Branch("EdepAll",&EdepAll,"EdepAll/D");
	anaTree->Branch("MWPCEnergy",fMWPCEnergy,"MWPCEnergyE/D:MWPCEnergyW/D");
	anaTree->Branch("time",hitTime,"timeE/D:timeW/D");
	anaTree->Branch("tmTime",trapMonTime,"tmTimeE/D:tmTimeW/D");
	
	char tmp[1024];
	sprintf(tmp,"EdepSD[%i]/D",N_SD);
	anaTree->Branch("EdepSD",EdepSD,tmp);
	sprintf(tmp,"thetaInSD[%i]/D",N_SD);
	anaTree->Branch("thetaInSD",thetaInSD,tmp);
	sprintf(tmp,"thetaOutSD[%i]/D",N_SD);
	anaTree->Branch("thetaOutSD",thetaOutSD,tmp);
	sprintf(tmp,"keInSD[%i]/D",N_SD);
	anaTree->Branch("keInSD",keInSD,tmp);
	sprintf(tmp,"keOutSD[%i]/D",N_SD);
	anaTree->Branch("keOutSD",keOutSD,tmp);
	sprintf(tmp,"hitCountSD[%i]/I",N_SD);
	anaTree->Branch("hitCountSD",hitCountSD,tmp);
	
	anaTree->Branch("MWPCPos",MWPCPos,"MWPCPosE[3]/D:MWPCPosW[3]/D");
	anaTree->Branch("ScintPos",ScintPos,"ScintPosE[3]/D:ScintPosW[3]/D");
	anaTree->Branch("MWPCPosSigma",MWPCPosSigma,"MWPCPosSigmaE[3]/D:MWPCPosSigmaW[3]/D");
	anaTree->Branch("ScintPosSigma",ScintPosSigma,"ScintPosSigmaE[3]/D:ScintPosSigmaW[3]/D");
}

void UCNA_MC_Analyzer::resetAnaEvt() {
	for(Side s = EAST; s <= WEST; ++s) {
		Edep[s] = EdepQ[s] = 0;
		fMWPCEnergy[s]=0;
		hitTime[s]=trapMonTime[s]=FLT_MAX;
		for(AxisDirection d=X_DIRECTION; d<=Z_DIRECTION; ++d)
			MWPCPos[s][d] = ScintPos[s][d] = MWPCPosSigma[s][d] = ScintPosSigma[s][d] = 0;
	}
	EdepAll = 0;
	for(size_t ii=0; ii<N_SD; ii++) {
		hitCountSD[ii] = keInSD[ii] = keOutSD[ii] = thetaInSD[ii] = thetaOutSD[ii] = EdepSD[ii] = 0.;
		hitTimeSD[ii]=FLT_MAX;
	}
}

void UCNA_MC_Analyzer::processTrack() {
	
	// detector ID numbers
	const int ID_scint[2] = {0,10};
	const int ID_mwpc[2] = {6,16};
	const int ID_deadmwpc[2] = {8,18};
	const int ID_winOut[2] = {3,13};
	const int ID_winIn[2] = {4,14};
	
	const int ID_trapmon[2] = {21,22};
	// particle ID numbers
	const int PDG_ELECTRON = 11;
	const int PDG_POSITRON = -11;
	
	const double kMe = 510.998910; // electron mass, keV
	
	// total deposited energy in all sensitive volumes
	EdepAll += trackinfo->Edep;
	if(detectorID < N_SD) {
		EdepSD[detectorID] += trackinfo->Edep;
		if(trackinfo->isEntering && pID==PDG_ELECTRON)
			hitCountSD[detectorID]++;
		if(pID==PDG_ELECTRON && trackinfo->hitTime<hitTimeSD[detectorID]) {
			hitTimeSD[detectorID] = trackinfo->hitTime;
			double pin_x = trackinfo->pIn[0];
			double pin_y = trackinfo->pIn[1];
			double pin_z = trackinfo->pIn[2];
			double pout_x = trackinfo->pOut[0];
			double pout_y = trackinfo->pOut[1];
			double pout_z = trackinfo->pOut[2];
			if(pin_x || pin_y || pin_z) {
				double magpin2 = pin_x*pin_x+pin_y*pin_y+pin_z*pin_z;
				keInSD[detectorID] = sqrt(magpin2+kMe*kMe)-kMe;
				thetaInSD[detectorID] = acos(fabs(pin_z/sqrt(magpin2)));
			}
			if(pout_x || pout_y || pout_z) {
				double magpout2 = pout_x*pout_x+pout_y*pout_y+pout_z*pout_z;
				keOutSD[detectorID] = sqrt(magpout2+kMe*kMe)-kMe;
				thetaOutSD[detectorID] = acos(fabs(pout_z/sqrt(magpout2)));
			}
		}
	}
	
	for(Side s = EAST; s <= WEST; ++s) {
		
		// scintillator deposited energy, position, hit time
		if(detectorID==ID_scint[s]) {
			Edep[s] += trackinfo->Edep;
			EdepQ[s] += trackinfo->EdepQuenched;
			for(AxisDirection d=X_DIRECTION; d<=Z_DIRECTION; ++d) {
				ScintPos[s][d] += trackinfo->edepPos[d];
				ScintPosSigma[s][d] += trackinfo->edepPos2[d];
			}
			//earliest hit time with possibly detectable (not really) Edep; hitTime[EAST,WEST] have been initialized to FLT_MAX
			if(trackinfo->Edep>5.0 && trackinfo->hitTime<hitTime[s]) 
				hitTime[s] = trackinfo->hitTime;
		}
		
		// trap monitor timing for electrons, initialized to FLT_MAX
		if(detectorID==ID_trapmon[s] && (pID==PDG_ELECTRON||pID==PDG_POSITRON)){
			if(trackinfo->hitTime<trapMonTime[s]) 
				trapMonTime[s] = trackinfo->hitTime; //earliest hit time 
		}
		
		// wirechamber deposited energy and position
		double anodeScale = 0;
		if(detectorID==ID_deadmwpc[s]) // contribution from "dead" gas
			anodeScale = 0.5;
		if(detectorID==ID_mwpc[s])
			anodeScale = 1.0;
		if(anodeScale) {
			fMWPCEnergy[s] += anodeScale*trackinfo->Edep;
			for(AxisDirection d=X_DIRECTION; d<=Z_DIRECTION; ++d) {
				MWPCPos[s][d] += anodeScale*trackinfo->edepPos[d];
				MWPCPosSigma[s][d] += anodeScale*trackinfo->edepPos2[d];
			}
		}
	}
}

void UCNA_MC_Analyzer::processEvent() {
	// normalize position variables
	for(Side s = EAST; s <= WEST; ++s) {
		for(AxisDirection d=X_DIRECTION; d<=Z_DIRECTION; ++d) {
			if(EdepQ[s]>0) {
				ScintPos[s][d] /= Edep[s];
				ScintPosSigma[s][d] = sqrt(ScintPosSigma[s][d]/EdepQ[s]-ScintPos[s][d]*ScintPos[s][d]);
			}
			if(fMWPCEnergy[s] > 0) {
				MWPCPos[s][d] /= fMWPCEnergy[s];
				MWPCPosSigma[s][d] = sqrt(MWPCPosSigma[s][d]/fMWPCEnergy[s]-MWPCPos[s][d]*MWPCPos[s][d]);
			}
		}
	}
}

int main(int argc, char** argv) {
	
	if(argc<3) {
		cout<<"Syntax: "<<argv[0]<<" <filename containing list of raw root files> <output root file name> [save all evts]"<<endl;
		exit(1);
	}
	
	UCNA_MC_Analyzer UMA(argv[2]);
	UMA.saveAllEvents = (argc==4);
	UMA.analyzeFileList(argv[1]);
	
	return 0;
}
