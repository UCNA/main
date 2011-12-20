#include "WirechamberCalibrator.hh"


WirechamberCalibrator::WirechamberCalibrator(RunNum rn, CalDB* cdb): anodeP(cdb->getAnodePositioningCorrector(rn)) {
	assert(anodeP);
	for(Side s = EAST; s <= WEST; s = nextSide(s))
		anodeGainCorr[s]=cdb->getAnodeGain(rn,s);
}

float WirechamberCalibrator::wirechamberGainCorr(Side s, float) const {  
	assert(s==EAST||s==WEST);
	return anodeGainCorr[s];
}

float WirechamberCalibrator::calibrateAnode(float adc, Side s, float x, float y, float t) const {
	assert(s==EAST||s==WEST);
	return adc*wirechamberGainCorr(s,t)/anodeP->eval(s,nBetaTubes,x,y,false);
}

void WirechamberCalibrator::printWirecalSummary() const {
	printf("Anode:\t\tcE=%.2f\tcW=%.2f\n",wirechamberGainCorr(EAST,0),wirechamberGainCorr(WEST,0));
}
