#ifndef __MCEvent_hh__
#define __MCEvent_hh__

#include <TObject.h>
#include <TClonesArray.h>

#include "TrackInfo.hh"
#include "PrimaryInfo.hh"

/// one MC event, holding lists of processed tracks and primary events data
/// this is what we save to the MC output ROOT TTree
class MCEvent : public TObject {
public:
	/// constructor
	MCEvent();
	/// destructor
	~MCEvent();
	
	/// append a TrackInfo
	TrackInfo *AddTrackInfo(const TrackInfo& aTrackInfo);
	/// append a PrimaryInfo
	PrimaryInfo *AddPrimaryInfo(const PrimaryInfo& aPrimaryInfo);
	
	/// clear event data
	void  ClearEvent();
	
	Int_t eventID;						//< ID number for event
	Int_t trapped;						//< "trapped" event flag
	Double_t compTime;					//< computation time for event
	
	static const unsigned int reserveArraySize;	//< size of array to reserve in memory
	
	TClonesArray* trackInfo;			//< event tracks, with space allocated for 10000 tracks
	TClonesArray* primaryInfo;			//< primary info array, with space allocated for 10000 primaries
	
	ClassDef(MCEvent, 1);
};

#endif
