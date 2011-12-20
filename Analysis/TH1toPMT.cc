#include "TH1toPMT.hh"
#include <climits>

TH1toPMT::TH1toPMT(TH1* h): ProcessedDataScanner("",true), mySpectrum(h), stochasticEnergy(true), randomPositionRadius(-1),
nToSim(0), nSimmed(0) {
	for(Side s = EAST; s <= WEST; s = nextSide(s)) {
		PGen[s].setSide(s);
		PGen[s].larmorField = 0;
	}
	fSide = BOTH;
}

bool TH1toPMT::nextPoint() {	
	assert(false); // TODO put energy simulation here
	assert(fSide==EAST || fSide==WEST);
	assert(mySpectrum);
	if(randomPositionRadius>=0) {
		assert(false);
	}
	float en = 0;
	if(stochasticEnergy) {
		en = mySpectrum->GetRandom();
	} else {
		assert(false);
	}
	scints[fSide] = PGen[fSide].generate(en);
	if(PGen[fSide].triggered()) {
		fPID = PID_BETA;
		fType = TYPE_0_EVENT;
	} else {
		fPID = PID_SINGLE;
		fType = TYPE_IV_EVENT;
	}
	nSimmed++;
	return nSimmed<nToSim;
}

void TH1toPMT::startScan(unsigned int startRandom) {
	if(startRandom) nToSim = startRandom;
	else nToSim = INT_MAX;
	nSimmed = 0;
}

void TH1toPMT::setCalibrator(PMTCalibrator& PCal) {
	for(Side s = EAST; s <= WEST; s = nextSide(s))
		PGen[s].setCalibrator(&PCal);
	ActiveCal = &PCal;
}
