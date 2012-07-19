#include "RunSetScanner.hh"
#include "PathUtils.hh"
#include "CalDBSQL.hh"
#include "CalDBFake.hh"
#include "SMExcept.hh"
#include <cassert>
#include <stdio.h>
#include <stdlib.h>

RunSetScanner::RunSetScanner(const std::string& treeName, bool withCalibrators):
TChainScanner(treeName), ActiveCal(NULL), totalTime(0), withCals(withCalibrators) {
	nAFP[0]=nAFP[1]=0;
}

RunSetScanner::~RunSetScanner() { 
	for(std::map<RunNum,PMTCalibrator*>::iterator it=PCals.begin(); it !=PCals.end(); it++)
		delete it->second;
}

unsigned int RunSetScanner::addRuns(const std::vector<RunNum>& rns) {
	printf("\n------------------- Assembling %i runs into TChain... ",(int)rns.size()); fflush(stdout);
	unsigned int n = 0;
	for(std::vector<RunNum>::const_iterator it = rns.begin(); it != rns.end(); it++) {
		n+=addRun(*it);
		printf("*"); fflush(stdout);
	}
	printf("------------------- %i Runs, %i events found, %.2fh running.\n",getnFiles(),nEvents,totalTime[BOTH]/3600.0);
	return n;
}

void RunSetScanner::display() {
	printf("RunSetScanner: %i runs, %i events [%i Off, %i On]\n",getnFiles(),nEvents,nAFP[0],nAFP[1]);
	if(runlist.size()>=getnFiles()) {
		for(unsigned int i=0; i<getnFiles(); i++) {
			RunInfo R = CalDBSQL::getCDB()->getRunInfo(runlist[i]);
			printf("\tRun %i: %i events\t",runlist[i],nnEvents[i]);
			R.display();
		}
	}
}

void RunSetScanner::writeCalInfo(QFile& qout, std::string tag) {
	for(std::map<RunNum,PMTCalibrator*>::iterator it=PCals.begin(); it !=PCals.end(); it++)
		qout.insert(tag,it->second->calSummary());
}

void RunSetScanner::speedload(unsigned int e) {
	if(e < noffset || e-noffset >= nLocalEvents) {
		Tch->LoadTree(e);
		Tch->GetTree()->LoadBaskets();
		nLocalEvents = Tch->GetTree()->GetEntries();
		noffset = Tch->GetChainOffset();
		if((int)runlist.size()>Tch->GetTreeNumber())
			evtRun = runlist[Tch->GetTreeNumber()];
		else
			evtRun = Tch->GetTreeNumber();
		if(withCals) {
			std::map<RunNum,PMTCalibrator*>::iterator it = PCals.find(evtRun);
			if(it == PCals.end()) {
				SMExcept e("missingCalibration");
				e.insert("runNum",evtRun);
				throw(e);
			}			
			ActiveCal = it->second;
		}
		loadNewRun(evtRun);
	}
	Tch->GetTree()->GetEvent(e-noffset);
}

bool RunSetScanner::addRun(RunNum r) {
	std::string f = locateRun(r);
	if(f.size() && addFile(f)) {
		runlist.push_back(r);
		RunInfo R = CalDBSQL::getCDB()->getRunInfo(r);
		if(R.afpState <= AFP_ON) nAFP[R.afpState] += nnEvents.back();
		if(withCals)
			PCals.insert(std::make_pair(r,new PMTCalibrator(r)));
		BlindTime b = CalDBSQL::getCDB()->fiducialTime(r);
		if(!b[BOTH])
			printf("**** WARNING: Run %i has zero fiducial time!\n",r);
		totalTime += b;
		runTimes.add(r,b[BOTH]);
		return true;
	}
	printf("**** FAILED TO LOCATE analyzed data for run %i at '%s'! *****\n",r,f.c_str());
	return false;
}


