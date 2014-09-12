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
	Type = SimOut->Branch("Type", &type, "type/I");
	MWPCPosE = SimOut->Branch("MWPCPosE",&mwpcposE, "mwpcposE[0]/F:mwpcposE[1]");  
	MWPCPosW = SimOut->Branch("MWPCPosW",&mwpcposW, "mwpcposW[0]/F:mwpcposW[1]"); 
}

void SimTreePlugin::fillSimTree(ProcessedDataScanner& PDS, double weight) {
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
	  	  
}
  


