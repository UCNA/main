#ifndef __bmMCEvent_hh__
#define __bmMCEvent_hh__

#include "TObject.h"
#include "TClonesArray.h"
#include "bmTrackInfo.hh"
#include "bmPrimaryInfo.hh"

/// one MC event, holding lists of processed tracks and primary events data
/// this is what we save to the MC output ROOT TTree
class bmMCEvent : public TObject {
public:
	/// constructor
	bmMCEvent();
	/// destructor
	~bmMCEvent();
	
	/// append a bmTrackInfo
	bmTrackInfo *AddTrackInfo(const bmTrackInfo& aTrackInfo);
	/// append a bmPrimaryInfo
	bmPrimaryInfo *AddPrimaryInfo(const bmPrimaryInfo& aPrimaryInfo);
	
	/// clear event data
	void  ClearEvent();
	
	Int_t eventID;						//< ID number for event
	Int_t trapped;						//< "trapped" event flag
	Double_t compTime;					//< computation time for event
	
	static const unsigned int reserveArraySize;	//< size of array to reserve in memory
	
	TClonesArray* trackInfo;			//< event tracks, with space allocated for 10000 tracks
	TClonesArray* primaryInfo;			//< primary info array, with space allocated for 10000 primaries
	
	ClassDef(bmMCEvent, 1);
};

#endif
