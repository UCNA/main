#include "PostOfficialAnalyzer.hh"
#include "PathUtils.hh"
#include "CalDBSQL.hh"
#include <utility>

std::string PostOfficialAnalyzer::locateRun(RunNum r) {
	std::vector<std::string> fpaths;
	fpaths.push_back(getEnvSafe("UCNAOUTPUTDIR")+"/hists/spec_"+itos(r)+".root");
	for(std::vector<std::string>::iterator it = fpaths.begin(); it != fpaths.end(); it++)
		if(fileExists(*it))
			return *it;
	printf("*** Replay of run %i not found! ***\n",r);
	return "";
}

unsigned int PostOfficialAnalyzer::addRun(RunNum r) {
	std::string f = locateRun(r);
	if(f.size() && addFile(f)) {
		ProcessedDataScanner::addRun(r);
		if(withCals)
			PCals.insert(std::make_pair(r,new PMTCalibrator(r,CDB)));
		BlindTime b = CalDBSQL::getCDB()->fiducialTime(r);
		if(!b.t[BOTH])
			printf("**** WARNING: Run %i has zero fiducial time!\n",r);
		totalTime += b;
		runTimes.add(r,b.t[BOTH]);
		return 1;
	}
	printf("**** FAILED TO LOCATE analyzed data for run %i at '%s'! *****\n",r,f.c_str());
	return 0;
}

void PostOfficialAnalyzer::setReadpoints() {
	
	// reconstructed energy
	Tch->SetBranchAddress("Etrue",&Etrue);
	// event ID
	Tch->SetBranchAddress("PID",&fPID);
	Tch->SetBranchAddress("Side",&fSide);
	Tch->SetBranchAddress("Type",&fType);
	// clock
	Tch->SetBranchAddress("TimeE",&runClock.t[EAST]);
	Tch->SetBranchAddress("TimeW",&runClock.t[WEST]);
	runClock.t[BOTH]=runClock.t[NONE]=0.0;
	
	for(Side s=EAST; s<=WEST; s=nextSide(s)) {
		
		// wirechamber planes
		for(int p=X_DIRECTION; p<=Y_DIRECTION; p++) {
			Tch->SetBranchAddress((std::string(p==X_DIRECTION?"x":"y")+sideSubst("%cmpm",s)).c_str(),&wires[s][p]);
			Tch->SetBranchAddress((sideSubst("Cathodes_%c",s)+(p==X_DIRECTION?"x":"y")).c_str(),cathodes[s][p]);
		}
		
		// beta scintillators
		Tch->SetBranchAddress(sideSubst("Scint%c",s).c_str(),&scints[s]);
		
		// MWPC totals
		Tch->SetBranchAddress(sideSubst("Anode%c",s).c_str(),&mwpcs[s].anode);
		Tch->SetBranchAddress(sideSubst("CathSum%c",s).c_str(),&mwpcs[s].cathodeSum);
		Tch->SetBranchAddress(sideSubst("EMWPC_%c",s).c_str(),&mwpcEnergy[s]);
	}
}
