/// \file BetaScint.hh \brief Beta Scintillator detector subsystem
#ifndef BETASCINT_HH
/// make sure this file is only included once
#define BETASCINT_HH 1

#include "Subsystem.hh"
#include "Trigger.hh"
#include <TF1.h>
#include <vector>

/// maximum number of tubes in beta scintillator (4)
const unsigned int maxBetaTubes = 4;

/// Class for a beta scintillator detector
class BetaScint: public Subsystem {
public:
	/// constructor
	BetaScint(RunManager* T, Trigger* TG, Side s);
	/// desctructor
	virtual ~BetaScint() { delete(events); }
	
	/// register branches in output tree
	virtual void addOutBranches(TTree* T);
	/// fill output tree readout location
	virtual void fillEvent(UInt_t e);
		
	/// write final info to Stringmap after all processing is done
	virtual Stringmap finalWords();
	
	
	/// check whether event is a Chris Pulser event
	bool isPulserTrigger(unsigned int e);
	
	ScintEvent anEvent;				//< data location for output tree
	Float_t led_pd;					//< reference LED monitor photodiode readout for output tree
	ScintEvent* events;				//< event data for all events	
	std::vector<float*> tubedat;	//< shortcut to tube ADC data
	Float_t qMax[maxBetaTubes];		//< maximum non-overflow QADC readout after ped. subtr.

	Trigger* Trig;		//< event triggers
	
	/// find tube trigger efficiency thresholds
	void triggerThresholds(bool allPlots = true, const std::vector<bool>& otherCuts = std::vector<bool>());

private:
	
	/// run-specific configuration
	void specialize();
	/// Locate GMS peaks
	void gmsPlots();
};

#endif

