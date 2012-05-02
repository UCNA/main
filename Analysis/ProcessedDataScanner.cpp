#include "ProcessedDataScanner.hh"
#include "PathUtils.hh"
#include "CalDBSQL.hh"
#include "CalDBFake.hh"
#include "UCNAException.hh"
#include <cassert>
#include <stdio.h>
#include <stdlib.h>

ProcessedDataScanner::ProcessedDataScanner(const std::string& treeName, bool withCalibrators):
RunSetScanner(treeName,withCalibrators), runClock(0), physicsWeight(1.0), anChoice(ANCHOICE_A), fiducialRadius(45.0) { }

Stringmap ProcessedDataScanner::evtInfo() {
	Stringmap m;
	m.insert("EvtN",currentEvent);
	for(Side s = EAST; s <= WEST; ++s)
		m.insert(sideSubst("EQ%c",s),scints[s].energy.x);
	m.insert("PID",fPID);
	m.insert("side",sideSubst("%c",fSide));
	m.insert("type",fType);
	if(fSide <= WEST) {
		m.insert("x",wires[fSide][X_DIRECTION].center);
		m.insert("y",wires[fSide][Y_DIRECTION].center);
	}
	return m;
}

void ProcessedDataScanner::recalibrateEnergy() {
	assert(ActiveCal);
	for(Side s = EAST; s<=WEST; ++s) {
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

float ProcessedDataScanner::radius2(Side s) const {
	return (wires[s][X_DIRECTION].center*wires[s][X_DIRECTION].center
			+wires[s][Y_DIRECTION].center*wires[s][Y_DIRECTION].center);
}
