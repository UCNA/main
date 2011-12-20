/// \file BetaScint.hh \brief GMS/Calibration LED run tracker
#ifndef LEDTRACKER_HH
/// make sure this file is only included once
#define LEDTRACKER_HH 1

#include "Subsystem.hh"
#include "RunManager.hh"
#include "Trigger.hh"
#include <TF1.h>
#include <vector>

/// Class for a beta scintillator detector
class LEDTracker: public Subsystem {
public:
	/// constructor
	LEDTracker(RunManager* T,Trigger* TG,unsigned int scanSteps = 16, float stepSize = 2.0, float stepTime = 60.0);
	std::vector<float*> tubedat[2];	//< shortcut to tube ADC data
	Trigger* Trig;		//< event triggers
	Side s;				//< current data side
protected:
	/// run-specific configuration
	void specialize();
};

#endif

