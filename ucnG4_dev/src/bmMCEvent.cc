#include "bmMCEvent.hh"
#include <iostream>
#include <cassert>
#include "TClonesArray.h"

ClassImp(bmMCEvent)

const unsigned int bmMCEvent::reserveArraySize = 10000;

bmMCEvent::bmMCEvent() {
	primaryInfo = new TClonesArray("bmPrimaryInfo",reserveArraySize);	
	trackInfo = new TClonesArray("bmTrackInfo",reserveArraySize);
	trapped = compTime = 0;
}


bmMCEvent::~bmMCEvent() {
	ClearEvent();
	delete trackInfo; trackInfo = NULL;
	delete primaryInfo; primaryInfo = NULL;
}

void bmMCEvent::ClearEvent() {
	trackInfo->Clear("C");
	primaryInfo->Clear("C");
	trapped = compTime = 0;
}

bmTrackInfo* bmMCEvent::AddTrackInfo(const bmTrackInfo& aTrackInfo) {
	unsigned int nentries = trackInfo->GetEntriesFast();
	assert(nentries<reserveArraySize);
	return new((*trackInfo)[nentries])bmTrackInfo(aTrackInfo); // construct at pre-reserved spot at trackInfo end
}

bmPrimaryInfo *bmMCEvent::AddPrimaryInfo(const bmPrimaryInfo& aPrimaryInfo) {
	unsigned int nentries = primaryInfo->GetEntriesFast();
	assert(nentries<reserveArraySize);
    return new((*primaryInfo)[nentries])bmPrimaryInfo(aPrimaryInfo); // construct at pre-reserved spot at primaryInfo end
}
