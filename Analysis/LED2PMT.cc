#include "LED2PMT.hh"
#include <math.h>
#include <climits>

LED2PMT::LED2PMT(): Sim2PMT("") {
	brightToEnergy[EAST] = brightToEnergy[WEST] = 2500.0;
	nEvents = INT_MAX;
}

double LED2PMT::ledBright() const {
	const double period = 60;
	const double flatprop = 0.10;
	double l = fmod(runClock[BOTH]+0.75*period,period)/period;
	if(l<flatprop) return 1e-12;
	l = (l-flatprop)/(1-flatprop);
	return l;
	//return exp(-5*(1.-l)+0*l);
}

void LED2PMT::startScan(bool) {
	runClock = BlindTime(0);
	currentEvent = nCounted = 0;
	if(nToSim<INT_MAX)
		printf("Scanning synthesized data for %i points... ",nToSim);
	else
		printf("Scanning unlimited synthesized data... ");
	fflush(stdout);
}

void LED2PMT::gotoEvent(unsigned int e) {
	currentEvent = e-1;
	nextPoint();
}

bool LED2PMT::nextPoint() {
	currentEvent++;
	reverseCalibrate();
	fPID = PID_LED;
	nSimmed++;
	nCounted+=simEvtCounts();
	if(nToSim<INT_MAX && !(currentEvent%(nToSim/20))) {
		printf("*"); fflush(stdout);
	}
	if(currentEvent >= nToSim) {
		printf("\n");
		return false;
	}
	return true;
}

void LED2PMT::doUnits() {
	double lt = ledBright();
	for(Side s = EAST; s <= WEST; ++s) {
		eW[s] = time[s] = 0;
		for(AxisDirection d = X_DIRECTION; d <= Z_DIRECTION; ++d)
			scintPos[s][d] = mwpcPos[s][d] = primPos[d] = 0;
		eQ[s] = eDep[s] = brightToEnergy[s]*lt;
	}
	costheta = ePrim = 0;
	primSide = BOTH;
}
