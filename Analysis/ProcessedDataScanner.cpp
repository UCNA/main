#include "ProcessedDataScanner.hh"
#include "PathUtils.hh"
#include "CalDBSQL.hh"
#include "SMExcept.hh"
#include <stdio.h>
#include <stdlib.h>

bool ProcessedDataScanner::redoPositions = false;

ProcessedDataScanner::ProcessedDataScanner(const std::string& treeName, bool withCalibrators):
RunSetScanner(treeName,withCalibrators), runClock(0), physicsWeight(1.0), anChoice(ANCHOICE_A), fiducialRadius(45.0) {
	for(Side s = EAST; s<=WEST; ++s)
		for(AxisDirection d = X_DIRECTION; d <= Y_DIRECTION; ++d)
				wires[s][d].center = 0;
}

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

float ProcessedDataScanner::probTrig(Side s, unsigned int t) {
	smassert(PMTActiveCal);
	smassert(s<=WEST && t<=nBetaTubes);
	return PMTActiveCal->trigEff(s, t, scints[s].adc[0]);
}

void ProcessedDataScanner::recalibrateEnergy() {
	smassert(ActiveCal);
	for(Side s = EAST; s<=WEST; ++s) {
		if(redoPositions && fPID==PID_BETA && fSide==s) {
			for(AxisDirection d = X_DIRECTION; d <= Y_DIRECTION; ++d) {
				wires[s][d] = ActiveCal->calcHitPos(s,d,cathodes[s][d]);
				ActiveCal->tweakPosition(s,d,wires[s][d],scints[s].energy.x);
			}
		}
		ActiveCal->calibrateEnergy(s, wires[s][X_DIRECTION].center, wires[s][Y_DIRECTION].center, scints[s], runClock[s]);
		mwpcEnergy[s] = ActiveCal->wirechamberEnergy(s, wires[s][X_DIRECTION], wires[s][Y_DIRECTION], mwpcs[s]);
	}
}

bool ProcessedDataScanner::passesPositionCut(Side s) {
	return radius(s)<fiducialRadius;
}

float ProcessedDataScanner::getErecon() const {
	smassert(ActiveCal);
	return ActiveCal->Erecon(fSide,fType,scints[EAST].energy.x,scints[WEST].energy.x);
}

float ProcessedDataScanner::radius2(Side s) const {
	return (wires[s][X_DIRECTION].center*wires[s][X_DIRECTION].center
			+wires[s][Y_DIRECTION].center*wires[s][Y_DIRECTION].center);
}
