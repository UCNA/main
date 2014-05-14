#ifndef TCHAINSCANNER_HH
#define TCHAINSCANNER_HH

#include <TChain.h>
#include <string>
#include <vector>

/// class for assembling and scanning a TChain
class TChainScanner {
public:
	/// constructor
	TChainScanner(const std::string& treeName);
	/// destructor
	virtual ~TChainScanner() { delete(Tch); }
	
	/// add a file to the TChain
	virtual int addFile(const std::string& filename);
	
	/// start a "speed scan," possibly at a random entry number
	virtual void startScan(bool startRandom = false);
	/// jump scanner to specified event
	virtual void gotoEvent(unsigned int e);
	/// load identified "speed scan" point
	virtual void speedload(unsigned int e);	
	/// load next "speed scan" point
	virtual bool nextPoint();	
	/// get current speed scan point
	unsigned int getCurrentEvent() const { return currentEvent; }
	/// load data for given event number
	virtual void getEvent(unsigned int e) { Tch->GetEvent(e); }
	/// get named branch address
	TChain* getChain() { return Tch; }
	/// get branch
	TBranch* getBranch(const char* bname) { return Tch->GetBranch(bname); }
	/// get local event number
	unsigned int getLocal(unsigned int e) { return Tch->LoadTree(e); }
	/// get number of files
	virtual unsigned int getnFiles() const { return nFiles; }
		
	UInt_t nEvents;						///< number of events in current TChain
	
protected:
	
	/// over-write this in subclass to automaticlly set readout points on first loaded file
	virtual void setReadpoints() {}
	/// "string friendly" SetBranchAddress
	void SetBranchAddress(const std::string& bname, void* bdata);
	
	std::vector<unsigned int> nnEvents;	///< number of events in each loaded TChain;
	unsigned int nFiles;				///< get number of loaded files
	
	TChain* Tch;						///< TChain of relevant runs
	unsigned int currentEvent;			///< event number of current event in chain
	unsigned int noffset;				///< offset of current event relative to currently loaded tree
	unsigned int nLocalEvents;			///< number of events in currently loaded tree
};

#endif
