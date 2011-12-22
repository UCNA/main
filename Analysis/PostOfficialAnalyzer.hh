#ifndef POSTOFFICIALANALYZER_HH
#define POSTOFFICIALANALYZER_HH 1

#include "ProcessedDataScanner.hh"

/// class for assembling and scanning TChain of official replay data
class PostOfficialAnalyzer: public ProcessedDataScanner {
public:
	/// constructor
	PostOfficialAnalyzer(bool withCalibrators = false): ProcessedDataScanner("phys",withCalibrators) {}
	
	/// add run file
	virtual unsigned int addRun(RunNum r);
	
	/// find path to processed run .root file
	static std::string locateRun(RunNum r);	
	
	float Etrue;				//< reconstructed "true" energy
	float cathodes[2][2][16];	//< pedestal-subtracted cathode values
	
protected:
	/// set TChain branch data readpoints
	virtual void setReadpoints();
};

#endif
