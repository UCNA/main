#include "MCEvent.hh"
#include <iostream>
#include <cassert>
#include "TClonesArray.h"

ClassImp(MCEvent)

const unsigned int MCEvent::reserveArraySize = 10000;

MCEvent::MCEvent() {
	primaryInfo = new TClonesArray("PrimaryInfo",reserveArraySize);	
	trackInfo = new TClonesArray("TrackInfo",reserveArraySize);
	trapped = compTime = 0;
}


MCEvent::~MCEvent() {
	ClearEvent();
	delete trackInfo; trackInfo = NULL;
	delete primaryInfo; primaryInfo = NULL;
}

void MCEvent::ClearEvent() {
	trackInfo->Clear("C");
	primaryInfo->Clear("C");
	trapped = compTime = 0;
}

TrackInfo* MCEvent::AddTrackInfo(const TrackInfo& aTrackInfo) {
	unsigned int nentries = trackInfo->GetEntriesFast();
	assert(nentries<reserveArraySize);
	return new((*trackInfo)[nentries])TrackInfo(aTrackInfo); // construct at pre-reserved spot at trackInfo end
}

PrimaryInfo *MCEvent::AddPrimaryInfo(const PrimaryInfo& aPrimaryInfo) {
	unsigned int nentries = primaryInfo->GetEntriesFast();
	assert(nentries<reserveArraySize);
    return new((*primaryInfo)[nentries])PrimaryInfo(aPrimaryInfo); // construct at pre-reserved spot at primaryInfo end
}
