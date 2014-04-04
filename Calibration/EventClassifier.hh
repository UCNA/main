#ifndef EVENTCLASSIFIER_HH
#define EVENTCLASSIFIER_HH

#include "Enums.hh"
#include "Types.hh"

/// base class for event PID/type/side identification logic; does not know how to fill in initial inputs
class EventClassifier {
public:
	/// constructor
	EventClassifier();
	/// destructor
	virtual ~EventClassifier() {}

	/// classify event based on inputs
	virtual void classifyEvent();

	/// return event description string
	std::string eventClassDescription() const;

	/// overall wirechamber cut
	virtual inline bool passedMWPC(Side s) const { return fPassedCathMaxSum[s]; }
	/// combined wirechamber+scintillator cut
	virtual inline bool Is2fold(Side s) const { return passedMWPC(s) && fPassedScint[s]; }
	/// overall muon tag check
	virtual inline bool taggedMuon() const { return fTaggedBack[EAST] || fTaggedBack[WEST] || fTaggedDrift[EAST] || fTaggedDrift[WEST] || fTaggedTop[EAST]; }
	/// sis-tagged LED events
	virtual inline bool isLED() const {return SIS00 & (1<<7); }
	/// sis-tagged UCN Mon events
	virtual inline bool isUCNMon() const { return SIS00 & ((1<<2) | (1<<8) | (1<<9) | (1<<10) | (1<<11)); }
	/// specific UCN monitors
	virtual inline bool isUCNMon(unsigned int n) const { return (SIS00 & (1<<2)) && (SIS00 & (1<<(8+n))); }
	/// sis-tagged scintillator trigger events
	inline bool isScintTrigger() const { return SIS00 & 3; }
	
	/// Type I initial hit side determination --- needed in subclass
	virtual Side getFirstScint() const = 0;
	/// Type II/III separation probability --- needed in subclass
	virtual float getProbIII() const = 0;
	
			
	// classification input variables
	Int_t fPassedScint[BOTH];		//< whether passed scintillator 2-of-4 on each side
	Int_t fPassedAnode[BOTH];		//< whether passed anode cut on each side
	Int_t fPassedCath[BOTH];		//< whether passed cathode sum cut on each side
	Int_t fPassedCathMax[BOTH];		//< whether passed cathode max cut on each side
	Int_t fPassedCathMaxSum[BOTH];	//< whether passed cathode max sum on each side
	
	Int_t fTaggedBack[BOTH];		//< whether event was tagged by the muon backing veto on each side
	Int_t fTaggedDrift[BOTH];		//< whether event was tagged by muon veto drift tubes on each side
	Int_t fTaggedTop[BOTH];			//< whether event was tagged by top veto on each side (only meaningful on East)

	Int_t fEvnbGood;				//< Event number header good
	Int_t fBkhfGood;				//< Block header good
	Int_t SIS00;					//< DAQ trigger bits
	
	// classification output variables
	Side fSide;						//< event primary scintillator side
	EventType fType;				//< event backscatter type
	PID fPID;						//< event particle ID
	Float_t fProbIII;				//< event estimated probability of being a Type III backscatter
};

#endif
