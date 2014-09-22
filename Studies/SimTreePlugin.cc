#include "SimTreePlugin.hh"
#include "BetaSpectrum.hh"
#include "GraphicsUtils.hh"
#include <TProfile.h>
#include <TGraph.h>
#include <iostream>


SimTreePlugin::SimTreePlugin(OctetAnalyzer* OA): OctetAnalyzerPlugin(OA,"simasymmetry") {
	
	makeSimTree();
	
}
void SimTreePlugin::makeSimTree() {
        SimFile = new TFile("/extern/mabrow05/ucna/SimAnalysisOutput/SimAnalyzedData.root","RECREATE"); 
        SimOut = new TTree("SimOut","SimAnalysis");
	MWPCE = SimOut->Branch("MWPCE", &mwpce, "mwpce/F");
	MWPCW = SimOut->Branch("MWPCW", &mwpcw, "mwpcw/F");
	Erecon = SimOut->Branch("Erecon", &erecon, "erecon/F");
	//Calside = SimOut->Branch("Side",&calside,"calside/I");
	Type = SimOut->Branch("Type", &type, "type/I");
	MWPCPosE = SimOut->Branch("MWPCPosE",&mwpcposE, "mwpcposE[2]/F");  
	MWPCPosW = SimOut->Branch("MWPCPosW",&mwpcposW, "mwpcposW[2]/F"); 
	//primarySide = SimOut->Branch("primarySide",&primside,"primside/I");
	//primCosTheta = SimOut->Branch("primCosTheta", &costheta,"costheta/F");
	//QuenchedE = SimOut->Branch("QuenchedE",&equenched,"equenched[2]/F");
	//PrimaryE = SimOut->Branch("PrimaryE",&primE,"primE/F");
	//Edep = SimOut->Branch("Edep", &edep,"edep[2]/F");
	//Make a struct which holds calibrated and simulated values. Then
	//just have two branches!
	//CalBranch = SimOut->Branch("RevCalData",
}

/*void SimTreePlugin::fillSimTree(ProcessedDataScanner& PDS, double weight) {
  //Fill all the simulation branches 
        smassert(PDS.isSimulated());
	Sim2PMT& S2P = (Sim2PMT&)PDS;
	Side s = S2P.fSide;
	if(!(s==EAST || s==WEST)) return;
	if(S2P.fPID != PID_BETA) return;
	if(S2P.passesPositionCut(s)) {
	        mwpce = S2P.mwpcEnergy[EAST];
		mwpcw = S2P.mwpcEnergy[WEST];
		erecon = S2P.getErecon();
		mwpcposE[0] = S2P.mwpcPos[EAST][X_DIRECTION];
		mwpcposE[1] = S2P.mwpcPos[EAST][Y_DIRECTION];
		mwpcposW[0] = S2P.mwpcPos[WEST][X_DIRECTION];
		mwpcposW[1] = S2P.mwpcPos[WEST][Y_DIRECTION];
		type = S2P.fType==TYPE_0_EVENT?0:(S2P.fType==TYPE_I_EVENT?1:(S2P.fType==TYPE_II_EVENT?2:(S2P.fType==TYPE_III_EVENT?3:4)));
		SimOut->Fill();
	    }
	  	  
	    }*/

void SimTreePlugin::fillCoreHists(ProcessedDataScanner& PDS, double weight) {
  //Fill all the simulation branches 
        smassert(PDS.isSimulated());
	Sim2PMT& S2P = (Sim2PMT&)PDS;
	Side s = S2P.fSide;
	if(!(s==EAST || s==WEST)) return;
	if(S2P.fPID != PID_BETA) return;
	if(S2P.passesPositionCut(s)) {
	        mwpce = S2P.mwpcEnergy[EAST];
		mwpcw = S2P.mwpcEnergy[WEST];
		erecon = S2P.getErecon();
		mwpcposE[0] = S2P.mwpcPos[EAST][X_DIRECTION];
		mwpcposE[1] = S2P.mwpcPos[EAST][Y_DIRECTION];
		mwpcposW[0] = S2P.mwpcPos[WEST][X_DIRECTION];
		mwpcposW[1] = S2P.mwpcPos[WEST][Y_DIRECTION];
		//costheta = S2P.costheta;
		//equenched[0] = S2P.eQ[0];
		//equenched[1] = S2P.eQ[1];
		//edep[0] = S2P.eDep[0];
		//edep[1] = S2P.eDep[1];
		//primE = S2P.ePrim;
		type = S2P.fType==TYPE_0_EVENT?0:(S2P.fType==TYPE_I_EVENT?1:(S2P.fType==TYPE_II_EVENT?2:(S2P.fType==TYPE_III_EVENT?3:4)));
		//if (S2P.primSide==EAST) primside=0;
		//else primside=1;
		//if (s==EAST) calside=0;
		//else calside=1;
		SimOut->Fill();
	    }
	  	  
}
  


