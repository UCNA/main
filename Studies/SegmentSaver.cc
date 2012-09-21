#include "SegmentSaver.hh"
#include "Types.hh"
#include "PathUtils.hh"
#include "SMExcept.hh"
#include <TString.h>

TH1* SegmentSaver::tryLoad(const std::string& hname) {
	if(!fIn) return NULL;
	TH1* h = NULL;
	fIn->GetObject(hname.c_str(),h);
	if(!h) {
		if(ignoreMissingHistos) {
			printf("Warning: missing histogram '%s' in '%s'\n",hname.c_str(),inflname.c_str());
		} else {
			SMExcept e("fileStructureMismatch");
			e.insert("fileName",inflname);
			e.insert("objectName",hname);
			throw(e);
		}
	} else {
		addObject(h);
	}
	return h;
}

TH1* SegmentSaver::registerSavedHist(const std::string& hname, const std::string& title,unsigned int nbins, float xmin, float xmax) {
	assert(saveHists.find(hname)==saveHists.end());	// don't duplicate names!
	TH1* h = tryLoad(hname);
	if(!h)
		h = registeredTH1F(hname,title,nbins,xmin,xmax);
	saveHists.insert(std::make_pair(hname,h));
	return h;
}

TH1* SegmentSaver::registerSavedHist(const std::string& hname, const TH1& hTemplate) {
	assert(saveHists.find(hname)==saveHists.end());	// don't duplicate names!
	TH1* h = tryLoad(hname);
	if(!h) {
		h = (TH1*)addObject(hTemplate.Clone(hname.c_str()));
		h->Reset();
	}
	saveHists.insert(std::make_pair(hname,h));
	return h;
}

SegmentSaver::SegmentSaver(OutputManager* pnt, const std::string& nm, const std::string& inflName):
OutputManager(nm,pnt), ignoreMissingHistos(false), inflname(inflName), inflAge(0) {		
	// open file to load existing data
	fIn = (inflname.size())?(new TFile((inflname+".root").c_str(),"READ")):NULL;
	assert(!fIn || !fIn->IsZombie());
	if(fIn) {
		inflAge = fileAge(inflname+".root");
		printf("Loading data from %s [%.1f hours old]...\n",inflname.c_str(),inflAge/3600.);
	}
}

SegmentSaver::~SegmentSaver() {
	if(fIn) {
		fIn->Close();
		delete(fIn);
	}
}

bool SegmentSaver::inflExists(const std::string& inflName) {
	return fileExists(inflName+".root") && fileExists(inflName+".txt");
}

TH1* SegmentSaver::getSavedHist(const std::string& hname) {
	std::map<std::string,TH1*>::iterator it = saveHists.find(hname);
	assert(it != saveHists.end());
	return it->second;
}

const TH1* SegmentSaver::getSavedHist(const std::string& hname) const {
	std::map<std::string,TH1*>::const_iterator it = saveHists.find(hname);
	assert(it != saveHists.end());
	return it->second;
}

void SegmentSaver::zeroSavedHists() {
	for(std::map<std::string,TH1*>::iterator it = saveHists.begin(); it != saveHists.end(); it++)
		it->second->Reset();
}

void SegmentSaver::scaleData(double s) {
	if(s==1.) return;
	for(std::map<std::string,TH1*>::iterator it = saveHists.begin(); it != saveHists.end(); it++)
		if(it->second->ClassName() != TString("TProfile"))
			it->second->Scale(s);
}

bool SegmentSaver::isEquivalent(const SegmentSaver& S) const {
	if(saveHists.size() != S.saveHists.size()) return false;
	for(std::map<std::string,TH1*>::const_iterator it = saveHists.begin(); it != saveHists.end(); it++) {
		std::map<std::string,TH1*>::const_iterator otherit = S.saveHists.find(it->first);
		if(otherit == S.saveHists.end()) return false;
		// TODO other checks?
	}
	return true;
}

void SegmentSaver::addSegment(const SegmentSaver& S) {
	assert(isEquivalent(S));
	// add histograms
	for(std::map<std::string,TH1*>::const_iterator it = saveHists.begin(); it != saveHists.end(); it++)
		it->second->Add(S.getSavedHist(it->first));
}
