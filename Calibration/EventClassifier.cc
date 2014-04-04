#include "EventClassifier.hh"
#include "WirechamberCalibrator.hh"

EventClassifier::EventClassifier(): fEvnbGood(true), fBkhfGood(true), SIS00(0) {
	for(Side s = EAST; s <= WEST; ++s) {
		fTaggedBack[s] = fTaggedDrift[s] = fTaggedTop[s] = false;
		fPassedScint[s] = fPassedAnode[s] = fPassedCath[s] = fPassedCathMax[s] = fPassedCathMaxSum[s] = false;
	}
}

std::string EventClassifier::eventClassDescription() const {
	return sideWords(fSide)+(" "+typeWords(fType)+" "+pidWords(fPID));
}

void EventClassifier::classifyEvent() {
	// particle ID: "beta-like" with scintillator and wirechamber, or "gamma-like"
	fPID = (Is2fold(EAST)||Is2fold(WEST)) ? (taggedMuon() ? PID_MUON:PID_BETA) : PID_SINGLE;

	// type, side
	fType = TYPE_IV_EVENT;
	fSide = NOSIDE;
	for(Side s = EAST; s<=WEST; ++s) {
		if(Is2fold(s)) {
			if(fPassedScint[otherSide(s)])
				fType = TYPE_I_EVENT;
			else
				fType = passedMWPC(otherSide(s))?TYPE_II_EVENT:TYPE_0_EVENT;
		}
		if(fPassedScint[s] && !fPassedScint[otherSide(s)]) fSide = s;
	}
	if(fPassedScint[WEST] && fPassedScint[EAST]) fSide = getFirstScint();
	
	// Type II/III separation
	fProbIII = (fType==TYPE_II_EVENT) ? getProbIII() : 0.;
	if(fProbIII>0.5) fType = TYPE_III_EVENT;
}
