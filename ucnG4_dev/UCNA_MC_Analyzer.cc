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
		for(unsigned int i=0; i<3; i++)
			MWPCPos[s][i] = ScintPos[s][i] = MWPCPosSigma[s][i] = ScintPosSigma[s][i] = 0;
	}
	EdepAll = 0;
	for(size_t ii=0; ii<N_SD; ii++) keInSD[ii] = keOutSD[ii] = thetaInSD[ii] = thetaOutSD[ii] = EdepSD[ii] = 0.;
}

void UCNA_MC_Analyzer::processTrack() {
	
	// detector ID numbers
	const int ID_scint[2] = {0,10};
	const int ID_mwpc[2] = {6,16};
		
	const int ID_trapmon[2] = {21,22};
	// particle ID numbers
	const int PDG_ELECTRON = 11;
	const int PDG_POSITRON = -11;
	
	const double kMe = 510.998910; // electron mass, keV
	
	// total deposited energy in all sensitive volumes
	EdepAll += trackinfo->Edep;
	if(detectorID < N_SD) {
		EdepSD[detectorID] += trackinfo->Edep;
		if(trackID==1 && pID==PDG_ELECTRON) {
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
			for(unsigned int j=0; j<3; j++) ScintPos[s][j] += trackinfo->edepPos[j];
			for(unsigned int j=0; j<3; j++) ScintPosSigma[s][j] += trackinfo->edepPos2[j];
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
		if(detectorID==ID_mwpc[s]) {
			fMWPCEnergy[s] += trackinfo->Edep;
			for(unsigned int j=0; j<3; j++) MWPCPos[s][j] += trackinfo->edepPos[j];
			for(unsigned int j=0; j<3; j++) MWPCPosSigma[s][j] += trackinfo->edepPos2[j];
		}
	}
}

void UCNA_MC_Analyzer::processEvent() {
	// normalize position variables
	for(Side s = EAST; s <= WEST; ++s) {
		for(int i = 0; i < 3; i++) {
			if(EdepQ[s]>0) {
				ScintPos[s][i] /= Edep[s];
				ScintPosSigma[s][i] = sqrt(ScintPosSigma[s][i]/EdepQ[s]-ScintPos[s][i]*ScintPos[s][i]);
			}
			if(fMWPCEnergy[s] > 0) {
				MWPCPos[s][i] /= fMWPCEnergy[s];
				MWPCPosSigma[s][i] = sqrt(MWPCPosSigma[s][i]/fMWPCEnergy[s]-MWPCPos[s][i]*MWPCPos[s][i]);
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
