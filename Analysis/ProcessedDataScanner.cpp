#include "ProcessedDataScanner.hh"
#include "PathUtils.hh"
#include "CalDBSQL.hh"
#include "CalDBFake.hh"
#include "SMExcept.hh"
#include <cassert>
#include <stdio.h>
#include <stdlib.h>

ProcessedDataScanner::ProcessedDataScanner(const std::string& treeName, bool withCalibrators):
RunSetScanner(treeName,withCalibrators), redoPositions(false), runClock(0),
physicsWeight(1.0), anChoice(ANCHOICE_A), fiducialRadius(50.0) { }

Stringmap ProcessedDataScanner::evtInfo() {
	Stringmap m;
	m.insert("EvtN",currentEvent);
	for(Side s = EAST; s <= WEST; ++s)
		m.insert(sideSubst("EQ%c",s),scints[s].energy.x);
	for(Side s = EAST; s <= WEST; ++s)
		m.insert(sideSubst("MWPC%c",s),mwpcEnergy[s]);
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
		if(redoPositions && fPID==PID_BETA && fSide==s) {
			std::vector<float> nopeds;
			for(AxisDirection d = X_DIRECTION; d <= Y_DIRECTION; ++d) {
				std::vector<float> caths(cathodes[s][d],cathodes[s][d]+kMaxCathodes);
				wires[s][d] = ActiveCal->calcHitPos(s,d,caths,nopeds);
				ActiveCal->tweakPosition(s,d,wires[s][d],scints[s].energy.x);
			}
		}
		ActiveCal->calibrateEnergy(s, wires[s][X_DIRECTION].center, wires[s][Y_DIRECTION].center, scints[s], runClock[s]);
		mwpcEnergy[s] = ActiveCal->calibrateAnode(mwpcs[s].anode,s,wires[s][X_DIRECTION].center, wires[s][Y_DIRECTION].center, runClock[s]);
	}
	calcEventFlags();
}

bool ProcessedDataScanner::passesPositionCut(Side s) {
	return radius(s)<fiducialRadius;
}

float ProcessedDataScanner::getEtrue() {
	assert(ActiveCal);
	return ActiveCal->Etrue(fSide,fType,scints[EAST].energy.x,scints[WEST].energy.x);
}

float ProcessedDataScanner::radius2(Side s) const {
	return (wires[s][X_DIRECTION].center*wires[s][X_DIRECTION].center
			+wires[s][Y_DIRECTION].center*wires[s][Y_DIRECTION].center);
}
