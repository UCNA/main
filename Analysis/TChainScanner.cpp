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

void TChainScanner::startScan(unsigned int startRandom) { 
	if(startRandom) {
		if(!currentEvent) {
			srand(time(NULL));	// random random seed
			currentEvent = rand()%Tch->GetEntries();
			Tch->GetEvent(currentEvent);
			assert(Tch->GetTree());
			printf("Scan Starting at offset %i/%i: ",currentEvent,nEvents);
		} else {
			printf("Scan Continuing at offset %i/%i: ",currentEvent,nEvents);
		}
	} else {
		Tch->GetEvent(0);
		assert(Tch->GetTree());
		currentEvent = -1;
		printf(">%i< ",nEvents);
	}
	fflush(stdout);
	nLocalEvents = noffset = 0;
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
