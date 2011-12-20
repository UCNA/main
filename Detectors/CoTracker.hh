#ifndef COTRACKER_HH
#define COTRACKER_HH 1

#include "Subsystem.hh"
#include "Enums.hh"
#include <stdlib.h>

/// UCNA data triggers
class CoTracker: public Subsystem {
public:
	/// constructor
	CoTracker(RunManager* T, Side s);
	
	/// get sis00 trigger flags for an event
	unsigned int sis00(unsigned int e) const { return (int)getData("Sis00")[e]; }
	/// if this event is quality cut (during beam, no beam, etc.)
	bool isCrud(UInt_t e) const { return f_bclock[e] < 0.08; }	
	/// if this is a reference Co60 source decay on my side
	bool gmsCo(UInt_t e) const { return (mySide == EAST && sis00(e) == 32) || (mySide == WEST && sis00(e) == 64); }

	float* f_refadc;	//< reference tube ped-subtracted adc value
	float* f_rclock;	//< time since run start, s
	float* f_bclock;	//< time since previous beam pulse, s
	
private:
		
	/// fix time wrap-around problem, convert times to seconds
	void fixTimes();

};

#endif
