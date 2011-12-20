/// \file MuonVeto.hh \brief Muon-veto detector subsystem
#ifndef MUONVETO_HH
/// make sure this file is only included once
#define MUONVETO_HH 1

#include "Subsystem.hh"
#include "Trigger.hh"

/// muon veto subsystem
class MuonVeto: public Subsystem {
public:
	/// constructor
	MuonVeto(RunManager* T, Trigger* tg);
	/// destructor
	virtual ~MuonVeto() { for(Side s = EAST; s <= WEST; s = nextSide(s)) if(striggers[s]) delete[](striggers[s]); }
	/// whether this event triggered a muon veto
	bool muonTrigger(UInt_t e) const { return striggers[EAST][e] || striggers[WEST][e]; }
	/// whether this event triggered the vetos on the given side
	bool muonTrigger(UInt_t e, Side s) const { assert(s<=WEST); return striggers[s][e]; }
	/// which side this event triggered muons on
	Side muonSide(UInt_t e) const { return sideCombo(striggers[EAST][e],striggers[WEST][e]); }
	/// whether this event triggered the backing veto on the specified side
	bool hitBacking(UInt_t e, Side s) const;
	
	Trigger* TG;		//< event triggers

private:
	/// generate ouput plots
	void genHistograms();
	/// run-specific configuration
	void specialize();
	
	unsigned int nchan;	//< number of ADC channels
	bool* striggers[2];	//< muon triggers on each side
};

#endif

