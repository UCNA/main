#ifndef POSTOFFICIALANALYZER_HH
#define POSTOFFICIALANALYZER_HH 1

#include "ProcessedDataScanner.hh"

/// class for assembling and scanning TChain of official replay data
class PostOfficialAnalyzer: public ProcessedDataScanner {
public:
	/// constructor
	PostOfficialAnalyzer(bool withCalibrators = false): ProcessedDataScanner("phys",withCalibrators) {}
	
	/// find path to processed run .root file
	virtual std::string locateRun(RunNum r);	
	
	float Etrue;	//< reconstructed "true" energy
	
protected:
	/// set TChain branch data readpoints
	virtual void setReadpoints();
};

#endif
