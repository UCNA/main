#include "CoTracker.hh"

CoTracker::CoTracker(RunManager* T, Side s): Subsystem(T,std::string("Co60_Tracker_")+ctos(sideNames(s)),s) {
	
	
	if(mySide == EAST) {
		if(thisRun->RI.runNum < 7000)
			loadData("Pdc21","ADCRefELED","reftube");	// GMS reference phototube
		else
			loadData("Pdc32","ADCRefELED","reftube");	// GMS reference phototube
	} else {
		if(thisRun->RI.runNum < 7000)
			loadData("Pdc216","ADCRefWLED","reftube");	// GMS reference phototube
		else			
			loadData("Pdc36","ADCRefWLED","reftube");	// GMS reference phototube		
	}
	
	loadData("Sis00","","Sis00");		// trigType
	loadData("S83028","","runclock");	// runclock
	loadData("S8200","","beamclock");	// protonClock	
	
	f_refadc = getData("reftube");
	f_rclock = getData("runclock");
	f_bclock = getData("beamclock");
	
	fixTimes();
	
	assert(verifyPedestal(0)); // need to have pedestals already generated for this run
	
	// pre-subtract pedestals
	for(UInt_t e=0; e<nEvents; e++)
		f_refadc[e] -= PC.getPedestal(sensorNames[0],f_rclock[e]);
	
}

void CoTracker::fixTimes() {
	
	// fix scaler reset jumps
	float deltaT = 0;
	int njumps = 0;
	for(unsigned int e=1; e<nEvents; e++) {
		if(f_rclock[e] < f_rclock[e-1]-deltaT-1e9) {
			deltaT += 4294967296.0;
			njumps++;
		}
		f_rclock[e] += deltaT;
	}
	
	if(njumps)
		printf("Fixed %i timing scaler overflows.\n",njumps);
	
	// convert us to s
	for(unsigned int e=0; e<nEvents; e++) {
		f_rclock[e] *= 1.0e-6;
		f_bclock[e] *= 1.0e-6;
	}
}
