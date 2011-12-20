#include "ProcessedDataScanner.hh"
#include "PathUtils.hh"
#include "CalDBSQL.hh"
#include "CalDBFake.hh"
#include <cassert>
#include <stdio.h>
#include <stdlib.h>

ProcessedDataScanner::ProcessedDataScanner(const std::string& treeName, bool withCalibrators):
TChainScanner(treeName), ActiveCal(NULL), totalTime(0),
anChoice(ANCHOICE_A), fiducialRadius(50.0), withCals(withCalibrators) {
	CDB = CalDBSQL::getCDB();
	if(!CDB->isValid(15926)) {
		printf("\n**** WARNING: Fake Calibrations in use!!! ****\n\n");
		assert(IGNORE_DEAD_DB);
		CDB = new CalDBFake();
	}
	runClock = 0;
	nAFP[0]=nAFP[1]=0;
}

ProcessedDataScanner::~ProcessedDataScanner() { 
	for(std::map<RunNum,PMTCalibrator*>::iterator it=PCals.begin(); it !=PCals.end(); it++)
		delete it->second;
}

unsigned int ProcessedDataScanner::addRuns(const std::vector<RunNum>& rns) {
	printf("\n------------------- Assembling %i runs into TChain... ",(int)rns.size()); fflush(stdout);
	unsigned int n = 0;
	for(std::vector<RunNum>::const_iterator it = rns.begin(); it != rns.end(); it++) {
		n+=addRun(*it);
		printf("*"); fflush(stdout);
	}
	printf("------------------- %i Runs, %i events found, %.2fh running.\n",getnFiles(),nEvents,totalTime.t[BOTH]/3600.0);
	return n;
}

void ProcessedDataScanner::display() {
	printf("ProcessedDataScanner: %i runs, %i events [%i Off, %i On]\n",getnFiles(),nEvents,nAFP[0],nAFP[1]);
	if(runlist.size()>=getnFiles()) {
		for(unsigned int i=0; i<getnFiles(); i++) {
			RunInfo R = CalDBSQL::getCDB()->getRunInfo(runlist[i]);
			printf("\tRun %i: %i events\t",runlist[i],nnEvents[i]);
			R.display();
		}
	}
}

void ProcessedDataScanner::writeCalInfo(QFile& qout, std::string tag) {
	for(std::map<RunNum,PMTCalibrator*>::iterator it=PCals.begin(); it !=PCals.end(); it++)
		qout.insert(tag,it->second->calSummary());
}

void ProcessedDataScanner::recalibrateEnergy() {
	assert(ActiveCal);
	for(Side s = EAST; s<=WEST; s=nextSide(s)) {
		ActiveCal->calibrateEnergy(s, wires[s][X_DIRECTION].center, wires[s][Y_DIRECTION].center, scints[s], runClock.t[BOTH]);
		mwpcEnergy[s] = ActiveCal->calibrateAnode(mwpcs[s].anode,s,wires[s][X_DIRECTION].center, wires[s][Y_DIRECTION].center, runClock.t[BOTH]);
	}
	calcEventFlags();
}

bool ProcessedDataScanner::passesPositionCut(Side s) { return radius(s)<fiducialRadius; }

float ProcessedDataScanner::getEtrue() {
	assert(ActiveCal);
	return ActiveCal->Etrue(fSide,fType,scints[EAST].energy.x,scints[WEST].energy.x);
}

void ProcessedDataScanner::speedload(unsigned int e) {
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
				printf("\n*** Bad run number for calibrations: %i at %f\n\n",evtRun,runClock.t[BOTH]);
				assert(false);
			}			
			ActiveCal = it->second;
		}
	}
	Tch->GetTree()->GetEvent(e-noffset);
}

unsigned int ProcessedDataScanner::addRun(RunNum rn) {
	runlist.push_back(rn);
	RunInfo R = CalDBSQL::getCDB()->getRunInfo(rn);
	if(R.afpState == AFP_OFF) nAFP[0] += nnEvents.back();
	else if(R.afpState == AFP_ON) nAFP[1] += nnEvents.back();
	return 1;
}

float ProcessedDataScanner::radius2(Side s) const {
	return (wires[s][X_DIRECTION].center*wires[s][X_DIRECTION].center
			+wires[s][Y_DIRECTION].center*wires[s][Y_DIRECTION].center);
}
