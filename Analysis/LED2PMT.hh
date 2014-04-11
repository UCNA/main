#ifndef LED2PMT_HH
#define LED2PMT_HH

#include "Sim2PMT.hh"

/// Simulates LED scan output
class LED2PMT: public Sim2PMT {
public:
	/// constructor
	LED2PMT();
	
	/// unit conversions to produce simulated data
	virtual void doUnits();
	/// overrides Sim2PMT::startScan to skip meddling with tree
	virtual void startScan(bool startRandom = false);
	/// overrides Sim2PMT::nextPoint to grab data from internal source
	virtual bool nextPoint();
	/// jump scanner to specified event
	virtual void gotoEvent(unsigned int e);
	/// whether to count this event as successfully generated
	virtual double simEvtCounts() const { return physicsWeight; }

	/// calculate LED brightness for current event
	virtual double ledBright() const;
	/// update event clock
	virtual void updateClock() { runClock = BlindTime(0.01*currentEvent); }

	double brightToEnergy[2];	///< conversion factor from LED brightness to equivalent energy on each side
};

#endif
