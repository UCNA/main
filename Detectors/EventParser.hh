/// \file EventParser.hh \brief Utility class for parsing event classification flags
#ifndef EVENTPARSER_HH
/// make sure this file is only included once
#define EVENTPARSER_HH 1

#include "Types.hh"
#include "OutputManager.hh"
#include <vector>
#include <TH1F.h>

/// parse event id flags into useful information
class EventParser {
public:
	/// constructor
	EventParser(int* fE, int* fW) { setFlagsReadPoint(fE,fW); }
	
	/// which side the scintillators triggered on
	inline Side whichScintSide(UInt_t e = 0) const { return sideCombo( evtF[EAST][e] & HIT_SCINT, evtF[WEST][e] & HIT_SCINT ); }
	/// whether the event was tracked on the given wire side
	inline bool hitWireSide(Side s, UInt_t e = 0) const { return evtF[s][e] & HIT_WIRES; }
	/// which side the MWPCs triggered on
	inline Side whichWiresSide(UInt_t e = 0) const { return sideCombo(hitWireSide(EAST,e), hitWireSide(WEST,e)); }
	/// whether the event was fully tracked on the specified side
	inline bool hitSide(Side s, UInt_t e = 0) const { return (evtF[s][e] & HIT_SCINT) && (evtF[s][e] & HIT_WIRES); }
	/// which side(s) the event was fully tracked on
	inline Side whichSide(UInt_t e = 0) const { return sideCombo(hitSide(EAST,e),hitSide(WEST,e)); }
	/// which side was hit first (e.g. TDC in zero-time peak)
	inline Side whichSideFirst(UInt_t e = 0) const { return sideCombo(evtF[EAST][e] & HIT_SIDE_FIRST, evtF[WEST][e] & HIT_SIDE_FIRST); }
	
	/// whether this is a GMS event
	inline bool isGMS(UInt_t e = 0) const { return (evtF[EAST][e] | evtF[WEST][e]) & IS_GMS; }
	/// whether the event is a muon
	inline bool isMuon(UInt_t e = 0) const { return (evtF[EAST][e] | evtF[WEST][e]) & IS_MUON; }
	/// whether the event is any non-beta-decay class
	inline bool isNotBeta(UInt_t e = 0) const { return (evtF[EAST][e] | evtF[WEST][e]) & IS_NOT_BETADATA; }
	
	/// whether the event is a "correct" (usually Beta decay) event
	inline bool isType0Beta(Side s, UInt_t e = 0) const { return (!isNotBeta(e)) && (evtF[s][e] & IS_CORRECT); }
	/// whether the event is a Type I backscatter
	inline bool isTypeIBeta(Side s, UInt_t e = 0) const { return (!isNotBeta(e)) && (evtF[s][e] & IS_TYPE_I); }
	/// whether the event is a Type II backscatter
	inline bool isTypeIIBeta(Side s, UInt_t e = 0) const { return (!isNotBeta(e)) && (evtF[s][e] & (IS_TYPE_II)); }
	/// whether the event is a Type III backscatter
	inline bool isTypeIIIBeta(Side s, UInt_t e = 0) const { return (!isNotBeta(e)) && (evtF[s][e] & (IS_TYPE_III)); }
	/// whether the event is a Type II or III backscatter
	inline bool isTypeII_IIIBeta(Side s, UInt_t e = 0) const { return (!isNotBeta(e)) && (evtF[s][e] & (IS_TYPE_II | IS_TYPE_III)); }
	/// whether the event is any type of beta/backscatter)
	inline bool isAnyBeta(Side s, UInt_t e = 0) const { return (!isNotBeta(e)) && (evtF[s][e] & IS_CLASSIFIED_HIT); }
	/// whether the event is hit-type classified
	inline bool isClassified(UInt_t e = 0) const { return (evtF[EAST][e] | evtF[WEST][e]) & IS_CLASSIFIED_HIT; }
	
	/// whether the original (SIMULATION) particle was going in this direction
	inline bool isPrimarySide(Side s, UInt_t e = 0) const { return evtF[s][e] & PRIMARY_DIRECTION; }
	
	/// classify event type based on hits information
	void classifyEventType(UInt_t e = 0) const;
	
	/// get flags for side
	int getFlags(Side s, UInt_t e = 0) const { return evtF[s][e]; }
	
	/// set event data readout point
	void setFlagsReadPoint(int* fE, int* fW) { evtF[EAST] = fE; evtF[WEST] = fW; }
	
private:
	int* evtF[2];	//< pointer to event classification flags
};

#endif
