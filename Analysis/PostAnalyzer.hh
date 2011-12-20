#ifndef POSTANALYZER_HH
#define POSTANALYZER_HH 1

#include "ProcessedDataScanner.hh"

/// Class for collecting pre-processed run data (potentially across multiple runs) for further analysis
class PostAnalyzer: public ProcessedDataScanner {
public:
	/// constructor
	PostAnalyzer(bool withCalibrators = false): ProcessedDataScanner("OutTree",withCalibrators) {}
	
	/// add run to chain
	virtual unsigned int addRun(RunNum r);
	
	TrigInfo trig;
	
protected:
	
	/// set branch data read points for TChain
	virtual void setReadpoints();
};

#endif
