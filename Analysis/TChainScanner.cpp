#include "TChainScanner.hh"
#include <cassert>
#include <stdlib.h>
#include <time.h>

TChainScanner::TChainScanner(const std::string& treeName): nEvents(0), nFiles(0), Tch(new TChain(treeName.c_str())),
currentEvent(0), noffset(0), nLocalEvents(0) {
	Tch->SetMaxVirtualSize(100000000);
}


int TChainScanner::addFile(const std::string& filename) {
	unsigned int oldEvents = nEvents;
	int nfAdded = Tch->Add(filename.c_str());
	if(!nfAdded) {
		printf("*** No such file: '%s'\n",filename.c_str());
		return 0;
	}
	nEvents = Tch->GetEntries();
	nnEvents.push_back(nEvents-oldEvents);
	if(!nFiles)
		setReadpoints();
	nFiles+=nfAdded;
	return nfAdded;
}

void TChainScanner::gotoEvent(unsigned int e) {
	currentEvent = e;
	Tch->GetEvent(currentEvent);
	assert(Tch->GetTree());
	nLocalEvents = noffset = 0;
}

void TChainScanner::startScan(bool startRandom) { 
	if(!nEvents) {
		printf("Starting scan with no data... ");
		fflush(stdout);
		return;
	}
	if(startRandom) {
		if(!currentEvent) {
			srand(time(NULL));	// random random seed
			gotoEvent(rand()%Tch->GetEntries());
			printf("Scan Starting at offset %i/%i: ",currentEvent,nEvents);
		} else {
			printf("Scan Continuing at offset %i/%i: ",currentEvent,nEvents);
		}
	} else {
		gotoEvent(0);
		currentEvent = -1;
		printf(">%i< ",nEvents);
	}
	fflush(stdout);
}

void TChainScanner::SetBranchAddress(const std::string& bname, void* bdata) {
	assert(bdata);
	Tch->SetBranchAddress(bname.c_str(),bdata);
}

void TChainScanner::speedload(unsigned int e) {
	if(e < noffset || e-noffset >= nLocalEvents) {
		Tch->LoadTree(e);
		Tch->GetTree()->SetMaxVirtualSize(10000000);
		//Tch->GetTree()->LoadBaskets();
		nLocalEvents = Tch->GetTree()->GetEntries();
		noffset = Tch->GetChainOffset();
	}
	Tch->GetTree()->GetEvent(e-noffset);
}

bool TChainScanner::nextPoint() {
	if(!nEvents) return false;
	++currentEvent;
	if(currentEvent >= nEvents) {
		printf("\n");
		startScan();
		return false;
	}
	if(!(currentEvent%(nEvents/20))) {
		printf("*"); fflush(stdout);
	}
	speedload(currentEvent);
	return true;
}
