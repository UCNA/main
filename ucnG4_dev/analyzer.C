#include <vector>
#include <map>
#include <stdlib.h>
#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "Rtypes.h"
#include "TString.h"
#include "TRandom.h"
#include "TF1.h"
#include "TMath.h"
#include "TCanvas.h"
#include <math.h>
#include "bmMCEvent.hh"
#include <iostream>
#include <TH1.h>
#include <TVector3.h>
#include <fstream>
#include <cassert>
#include <cfloat>
#include "Sides.hh"

using namespace std;

const double kMe = 510.998910; //keV

Int_t fTrapped;				// whether this particle is trapped in mag field (lives too long)
Double_t fCompTime;			// computer time needed to track this event

Double_t Edep[2];			// scintillator deposited energy
Double_t EdepQ[2];			// quenched energy in scintillator
Double_t MWPCPos[2][3];		// MWPC deposited energy weighted position
Double_t ScintPos[2][3];	// scintillator quenched energy weighted position	
Double_t MWPCPosSigma[2][3];	// MWPC deposited energy weighted position variance
Double_t ScintPosSigma[2][3];	// scintillator quenched energy weighted position variance

Double_t fMWPCEnergy[2];	// MWPC deposited energy
Double_t EdepSD[19];		// array for energy deposition in all SDs
Double_t EdepAll;			// total edep in all SDs
Double_t hitTime[2];		// timing info for hits on each side
Double_t trapMonTime[2];	// timing for trap monitor hits
Double_t primKE;			// primary event kinetic energy
Double_t primTheta;			// primary event emission angle
Double_t primPos[4];		// primary event position (4th coordinate = radius)
Double_t thetaIn[2];		// particle angle entering wirechamber
Double_t thetaOut[2];		// particle angle exiting wirechamber
Double_t kEIn[2];			// particle energy entering scintillator
Double_t kEInTrapMon[2];	// particle energy entering decay trap monitor volume
Double_t kEOut[2];			// particle energy leaving scintillator
Long_t seed;				// random seed used for this event


void GetListOfFiles(const TString filename, vector<TString>& filelist) {
	ifstream file;
	TString stmp;
	file.open(filename.Data());
	while(file){
		file>>stmp;
		if(!stmp.IsNull()) filelist.push_back(stmp);
		else break;
	}
	return;
}


void ResetAnaEvt() {
	for(Side s = EAST; s <= WEST; ++s) {
		Edep[s] = EdepQ[s] = 0;
		fMWPCEnergy[s]=0;
		hitTime[s]=trapMonTime[s]=FLT_MAX;
		thetaIn[s]=0;
		thetaOut[s]=0;
		kEIn[s]=0;
		kEInTrapMon[s]=0;
		kEOut[s]=0;
		for(unsigned int i=0; i<3; i++)
			MWPCPos[s][i] = ScintPos[s][i] = MWPCPosSigma[s][i] = ScintPosSigma[s][i] = 0;
	}
	fTrapped=0;
	primKE=0;
	primTheta=0;
	EdepAll = 0;
	for(size_t ii=0; ii<19; ii++) EdepSD[ii] = 0.;
	seed=0;
}

int main(int argc, char **argv){
	
	if(argc<3) {
		cout<<"Syntax: "<<argv[0]<<" <filename containing list of raw root files> <output root file name> [save all evts]"<<endl;
		exit(1);
	}
	TString inputfile = TString(argv[1]);
	TString outputfile = TString(argv[2]);
	bool saveAllEvents = (argc==4);
	
	gRandom->SetSeed(2011);	// note consistently set random seed
	
	TFile* fout = new TFile(outputfile,"RECREATE");
	TTree* anaTree = new TTree("anaTree", "tree for analysis");
	
	anaTree->Branch("Edep",Edep,"EdepE/D:EdepW/D");
	anaTree->Branch("EdepQ",EdepQ,"EdepQE/D:EdepQW/D");
	anaTree->Branch("EdepAll",&EdepAll,"EdepAll/D");
	anaTree->Branch("MWPCEnergy",fMWPCEnergy,"MWPCEnergyE/D:MWPCEnergyW/D");
	anaTree->Branch("time",hitTime,"timeE/D:timeW/D");
	anaTree->Branch("tmTime",trapMonTime,"tmTimeE/D:tmTimeW/D");
	anaTree->Branch("primKE",&primKE,"primKE/D");
	anaTree->Branch("primTheta",&primTheta,"primTheta/D");
	anaTree->Branch("primPos",primPos,"primPos[4]/D");
	//angle in and out of the wire chamber
	anaTree->Branch("thetaIn",thetaIn,"thetaInE/D:thetaInW/D");
	anaTree->Branch("thetaOut",thetaOut,"thetaOutE/D:thetaOutW/D");
	//Kinetic energy in and out of the wire chamber
	anaTree->Branch("kEIn",kEIn,"kEInE/D:kEInW/D");
	anaTree->Branch("kEOut",kEOut,"kEOutE/D:kEOutW/D");
	anaTree->Branch("kEInTrapMon",kEInTrapMon,"kEInTrapMonE/D:kEInTrapMonW/D");
	//energy loss
	anaTree->Branch("EdepSD",EdepSD,"EdepSD[19]/D");
	//energy weighted positions
	anaTree->Branch("MWPCPos",MWPCPos,"MWPCPosE[3]/D:MWPCPosW[3]/D");
	anaTree->Branch("ScintPos",ScintPos,"ScintPosE[3]/D:ScintPosW[3]/D");
	anaTree->Branch("MWPCPosSigma",MWPCPosSigma,"MWPCPosSigmaE[3]/D:MWPCPosSigmaW[3]/D");
	anaTree->Branch("ScintPosSigma",ScintPosSigma,"ScintPosSigmaE[3]/D:ScintPosSigmaW[3]/D");
	
	anaTree->Branch("trapped",&fTrapped,"trapped/I");
	anaTree->Branch("compTime",&fCompTime,"compTime/D");
	anaTree->Branch("seed",&seed,"seed/I");
	
	vector<TString> sfiles;
	//Get List of Files.
	GetListOfFiles(inputfile,sfiles);
	
	for(size_t nn=0;nn<sfiles.size();nn++){
		
		cout<<"opening "<<sfiles[nn].Data()<<endl;
		TFile *f = new TFile(sfiles[nn],"read");
		if(f->IsZombie()) {
			cout << "*** Error opening" << sfiles[nn].Data()<<endl;
			continue;
		}
		
		bmMCEvent* myevt = new bmMCEvent();
		
		TTree *tree = (TTree*) f->Get("EventTree");
		tree->SetBranchAddress("MC_event_output",&myevt);
		
		for(int ii=0;ii<tree->GetEntriesFast();ii++){
			
			tree->GetEntry(ii);
			
			TClonesArray &priminfo = *(myevt->primaryInfo);    
			TClonesArray &trackinfo = *(myevt->trackInfo);
			
			ResetAnaEvt(); //reset values from the last entry
			
			// detector ID numbers
			const int ID_scint[2] = {0,8};
			const int ID_mwpc[2] = {6,14};
			const int ID_mwpc_frontwin[2] = {4,12};
			const int ID_trapmon[2] = {17,18};
			// particle ID numbers
			const int PDG_ELECTRON = 11;
			const int PDG_POSITRON = -11;
			
			//get the size of array for this event first
			Int_t size = trackinfo.GetEntries();			
			fTrapped = myevt->trapped;
			fCompTime = myevt->compTime;
			
			//const double jitter = 2.5;
			double pin_x(0), pin_y(0), pin_z(0);
			double pout_x(0), pout_y(0), pout_z(0);

			for(int nn=0;nn<size;nn++){
				
				int pID = ((bmTrackInfo*)trackinfo[nn])->pID;
				int detectorID = ((bmTrackInfo*)trackinfo[nn])->hcID;
				int trackID = ((bmTrackInfo*)trackinfo[nn])->trackID;
				
				// total deposited energy in all sensitive volumes
				EdepAll += ((bmTrackInfo*)trackinfo[nn])->Edep;
				if(detectorID < 19)
					EdepSD[detectorID] += ((bmTrackInfo*)trackinfo[nn])->Edep;
				
				for(Side s = EAST; s <= WEST; ++s) {
					
					// scintillator deposited energy, position, hit time
					if(detectorID==ID_scint[s]) {
						Edep[s] += ((bmTrackInfo*)trackinfo[nn])->Edep;
						EdepQ[s] += ((bmTrackInfo*)trackinfo[nn])->EdepQuenched;
						for(unsigned int j=0; j<3; j++) ScintPos[s][j] += ((bmTrackInfo*)trackinfo[nn])->edepPos[j];
						for(unsigned int j=0; j<3; j++) ScintPosSigma[s][j] += ((bmTrackInfo*)trackinfo[nn])->edepPos2[j];
						//earliest hit time with possibly detectable (not really) Edep; hitTime[EAST,WEST] have been initialized to FLT_MAX
						if(((bmTrackInfo*)trackinfo[nn])->Edep>5.0 && ((bmTrackInfo*)trackinfo[nn])->hitTime<hitTime[s]) 
							hitTime[s] =((bmTrackInfo*)trackinfo[nn])->hitTime;
					}
					
					// trap monitor timing for electrons, initialized to FLT_MAX
					if(detectorID==ID_trapmon[s] && (pID==PDG_ELECTRON||pID==PDG_POSITRON)){
						if(((bmTrackInfo*)trackinfo[nn])->hitTime<trapMonTime[s]) 
							trapMonTime[s] = ((bmTrackInfo*)trackinfo[nn])->hitTime; //earliest hit time 
					}
					
					
					// wirechamber deposited energy and position
					if(detectorID==ID_mwpc[s]) {
						fMWPCEnergy[s] += ((bmTrackInfo*)trackinfo[nn])->Edep;
						for(unsigned int j=0; j<3; j++) MWPCPos[s][j] += ((bmTrackInfo*)trackinfo[nn])->edepPos[j];
						for(unsigned int j=0; j<3; j++) MWPCPosSigma[s][j] += ((bmTrackInfo*)trackinfo[nn])->edepPos2[j];
					}
					
					// entrance, exit variables for electron
					if(trackID==1 && pID==PDG_ELECTRON) {
						pin_x = ((bmTrackInfo*)trackinfo[nn])->pIn[0];
						pin_y = ((bmTrackInfo*)trackinfo[nn])->pIn[1];
						pin_z = ((bmTrackInfo*)trackinfo[nn])->pIn[2];
						pout_x = ((bmTrackInfo*)trackinfo[nn])->pOut[0];
						pout_y = ((bmTrackInfo*)trackinfo[nn])->pOut[1];
						pout_z = ((bmTrackInfo*)trackinfo[nn])->pOut[2];
						if(detectorID==ID_scint[s]) {
							if(pin_x*pin_y*pin_z)
								kEIn[s] = sqrt(pin_x*pin_x+pin_y*pin_y+pin_z*pin_z+kMe*kMe)-kMe;
							if(pout_x*pout_y*pout_z)
								kEOut[s] = sqrt(pout_x*pout_x+pout_y*pout_y+pout_z*pout_z+kMe*kMe)-kMe;
						} 
						if(detectorID==ID_trapmon[s] && pin_x*pin_y*pin_z) {
							kEInTrapMon[s] = sqrt(pin_x*pin_x+pin_y*pin_y+pin_z*pin_z+kMe*kMe)-kMe;
						}						
						if(detectorID==ID_mwpc_frontwin[s]) {
							if(pin_x*pin_y*pin_z)
								thetaIn[s] = acos(fabs(pin_z/sqrt(pin_x*pin_x+pin_y*pin_y+pin_z*pin_z)))*180./TMath::Pi();
							if(pout_x*pout_y*pout_z)
								thetaOut[s] = acos(fabs(pout_z/sqrt(pout_x*pout_x+pout_y*pout_y+pout_z*pout_z)))*180./TMath::Pi();	 
						}
					}
					
				}
				
			}    //finish looping hits in this event
			
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
			
			// primary event info
			primKE = ((bmPrimaryInfo*)priminfo[0])->KE;
			for(unsigned int i=0; i<3; i++)
				primPos[i] = ((bmPrimaryInfo*)priminfo[0])->vertex[i];
			primPos[3] = sqrt(primPos[0]*primPos[0]+primPos[1]*primPos[1]);
			TVector3 mom(((bmPrimaryInfo*)priminfo[0])->p[0],((bmPrimaryInfo*)priminfo[0])->p[1],((bmPrimaryInfo*)priminfo[0])->p[2]);
			primTheta = mom.Theta();
			seed = ((bmPrimaryInfo*)priminfo[0])->seed;
			
			// only fill events with a chance of triggering unless otherwise specified
			if(saveAllEvents || Edep[EAST]+Edep[WEST]>0) 
				anaTree->Fill();
			
		}//finished looping one run 
		
		delete myevt; myevt = NULL;
		f->Close();
	}//finished looping over runs
	
	//write analyzed tree
	fout->cd();
	anaTree->Write();
	fout->Close();
	    
	return 0;
}
