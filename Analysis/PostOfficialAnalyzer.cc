#include "PostOfficialAnalyzer.hh"
#include "PathUtils.hh"
#include "CalDBSQL.hh"
#include <utility>

std::string PostOfficialAnalyzer::locateRun(RunNum r) {
	std::string fname = getEnvSafe("UCNAOUTPUTDIR")+"/hists/spec_"+itos(r)+".root";
	if(fileExists(fname))
		return fname;
	printf("*** Replay of run %i not found! ***\n",r);
	return "";
}

void PostOfficialAnalyzer::setReadpoints() {
	
	// reconstructed energy
	SetBranchAddress("Erecon",&Erecon);
	// event ID
	SetBranchAddress("PID",&fPID);
	SetBranchAddress("Side",&fSide);
	SetBranchAddress("Type",&fType);
	if(Tch->GetBranch("ProbIII"))
		SetBranchAddress("ProbIII",&fProbIII);
	else
		fProbIII = 0;
	// clock
	SetBranchAddress("TimeE",&runClock[EAST]);
	SetBranchAddress("TimeW",&runClock[WEST]);
	runClock[BOTH]=runClock[NOSIDE]=0.0;
	
	SetBranchAddress("EvnbGood",&fEvnbGood);
	SetBranchAddress("BkhfGood",&fBkhfGood);

	for(Side s=EAST; s<=WEST; ++s) {
		
		// wirechamber planes
		for(int p=X_DIRECTION; p<=Y_DIRECTION; p++) {
			SetBranchAddress(((p==X_DIRECTION?"x":"y")+sideSubst("%cmpm",s)).c_str(),&wires[s][p]);
			SetBranchAddress((sideSubst("Cathodes_%c",s)+(p==X_DIRECTION?"x":"y")).c_str(),cathodes[s][p]);
		}
		
		// beta scintillators
		SetBranchAddress(sideSubst("Scint%c",s).c_str(),&scints[s]);
		
		// MWPC totals
		SetBranchAddress(sideSubst("Anode%c",s).c_str(),&mwpcs[s].anode);
		SetBranchAddress(sideSubst("CathSum%c",s).c_str(),&mwpcs[s].cathodeSum);
		SetBranchAddress(sideSubst("EMWPC_%c",s).c_str(),&mwpcEnergy[s]);
		
		/// muon tags
		SetBranchAddress(sideSubst("TaggedBack%c",s).c_str(),&fTaggedBack[s]);
	}
}
