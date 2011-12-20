#include "EventParser.hh"

void EventParser::classifyEventType(UInt_t e) const {
	for( Side s = EAST; s <= WEST; s = nextSide(s) ) {
		// Correct or Type-IV (unknown backscatter) event
		if( whichWiresSide(e) == s && whichScintSide(e) == s ) {
			evtF[s][e] |= IS_CORRECT;
			break;
		}
		// Type-I backscatter event, fires both scintillators and WCs
		if( whichSide(e) == BOTH && whichSideFirst(e) == s ) {
			evtF[s][e] |= IS_TYPE_I;
			break;
		}
		// Type-II/III backscatter event; final classification method for type III elsewhere
		if( whichScintSide(e) == s && whichWiresSide(e) == BOTH) {
			evtF[s][e] |= IS_TYPE_II;
			break;
		}
	}
}
