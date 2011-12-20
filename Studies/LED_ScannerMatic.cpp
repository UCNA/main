#include "RunManager.hh"
#include "LEDTracker.hh"

// .L LED_ScannerMatic.cpp++
// LED_ScannerMatic(

void LED_ScannerMatic(RunNum R, unsigned int nSteps, float stepSize, float stepTime) {
	printf("Processing LED scan for %i (%i x %is %.2fdB steps)\n",R,nSteps,int(stepTime),stepSize);
	RunManager TM(R,RUNMODE_FULLPLOTS);		
	Trigger TG(&TM);
	LEDTracker LT(&TM,&TG,nSteps,stepSize,stepTime);
	TM.write();
}
