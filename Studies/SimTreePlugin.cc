#include "SimTreePlugin.hh"
#include "BetaSpectrum.hh"
#include "GraphicsUtils.hh"
#include <TProfile.h>
#include <TGraph.h>
#include <iostream>


SimTreePlugin::SimTreePlugin(OctetAnalyzer* OA): OctetAnalyzerPlugin(OA,"simasymmetry") {
	 
  //makeSimTree(18436);
	
}
void SimTreePlugin::makeSimTree(int rn) {
        char temp[200];
        sprintf(temp,"/extern/mabrow05/ucna/SimAnalysisOutput/hists/%i_sim.root",rn);
        SimFile = new TFile(temp,"RECREATE"); 
	//SimFile = new TFile("/home/mabrow05/SimAnalysisOutput/SimAnalyzedDataStruct.root","RECREATE"); 
        SimOut = new TTree("SimOut","SimAnalysis");
	
	SimDataBranch = (TBranch*)SimOut->Branch("SimData",&SD,"eQ[2]/D:eDep[2]:eW[2]:scintPos[2][3]:mwpcPos[2][3]:primPos[4]:time[2]:costheta:ePrim:edepFoils[2]:edepWinOut[2]:edepWinIn[2]:edepDeadMWPC[2]:edepKevlar[2]:edepWires[2]:edepDeadScint[2]"); // 256000
	RevCalSimDataBranch = (TBranch*)SimOut->Branch("RevCalSimData",&SCD,"Erecon/D:physicsWeight:mwpcEnergyE/F:mwpcEnergyW:type/I:side:PID");
}


void SimTreePlugin::fillCoreHists(ProcessedDataScanner& PDS, double weight) {
  //Fill all the simulation branches 
        smassert(PDS.isSimulated());
	Sim2PMT& S2P = (Sim2PMT&)PDS;
	Side s = S2P.fSide;
	if(!(s==EAST || s==WEST)) return;
	if(S2P.fPID != PID_BETA) return;
	if(S2P.passesPositionCut(s)) {	  
	  // For Reverse Calibrated branch
	        SCD.mwpcEnergyE = S2P.mwpcEnergy[EAST];
		SCD.mwpcEnergyW = S2P.mwpcEnergy[WEST];
		SCD.Erecon = S2P.getErecon();
		SCD.type = S2P.fType==TYPE_0_EVENT?0:(S2P.fType==TYPE_I_EVENT?1:(S2P.fType==TYPE_II_EVENT?2:(S2P.fType==TYPE_III_EVENT?3:4)));
		SCD.side = s;
		SCD.PID = S2P.fPID;
		SCD.physicsWeight = S2P.physicsWeight;

		// For pure simulated variables branch
		SD.costheta = S2P.costheta;
	  	SD.ePrim = S2P.ePrim;
		
		for (Side ss = EAST; ss<=WEST; ++ss) {
		  for (AxisDirection d=X_DIRECTION; d<=Z_DIRECTION; ++d) {
		    SD.mwpcPos[ss][d] = S2P.mwpcPos[ss][d];
		    SD.scintPos[ss][d] = S2P.scintPos[ss][d];
		  }
		  SD.eQ[ss] = S2P.eQ[ss];
		  SD.eDep[ss] = S2P.eDep[ss];
		  SD.eW[ss] = S2P.eW[ss];	  
		  SD.time[ss] = S2P.time[ss];	
		  SD.edepFoils[ss] = S2P.edepFoils[ss];
		  SD.edepWinOut[ss] = S2P.edepWinOut[ss];
		  SD.edepWinIn[ss] = S2P.edepWinIn[ss];
		  SD.edepDeadMWPC[ss] = S2P.edepDeadMWPC[ss];
		  SD.edepKevlar[ss] = S2P.edepKevlar[ss];
		  SD.edepWires[ss] = S2P.edepWires[ss];
		  SD.edepDeadScint[ss] = S2P.edepDeadScint[ss];
		  
		}
		for (unsigned int i=0; i<4; i++) {
		  SD.primPos[i] = S2P.primPos[i];
		}
		//ofile.open("~/dataCheck.txt", std::ios_base::out | std::ios_base::app);
		//std::cout << S2P.edepDeadScint[0] << "\t" << S2P.edepDeadScint[1] << std::endl;
		//ofile.close();
		//printf(&ofile,"%f\t%f\n",&S2P.edepWinOut[0],&S2P.edepWinOut[1]);
		//if (S2P.primSide==EAST) SD.primSide=0;
		//else SD.primSide=1;
		
		SimOut->Fill();
	}
	  	  
}
  


