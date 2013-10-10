#ifndef OCTETSIMUCLONEMANAGER_HH
#define OCTETSIMUCLONEMANAGER_HH 1

#include "RunAccumulator.hh"
#include "PathUtils.hh"
#include <string>

/// class for handling octet-by-octet matching data and simulation analyses
class OctetSimuCloneManager {
public:
	/// constructor
	OctetSimuCloneManager(const std::string& dname, const std::string& bdir = getEnvSafe("UCNA_ANA_PLOTS"));
	/// destructor
	virtual ~OctetSimuCloneManager() {}
	
	/// scan one octet of data
	void scanOct(RunAccumulator& RA, unsigned int octn);
	/// combine all data octets
	void combineOcts(RunAccumulator& RA);
	/// simulate one octet
	void simOct(RunAccumulator& SimRA, unsigned int octn);
	/// combine simulated octets
	void combineSims(RunAccumulator& SimRA);
	
	/// load simulation data for octet
	virtual Sim2PMT* getSimdata(unsigned int octn);
	
	std::string outputDir;	//< output directory naming
	std::string baseDir;	//< base directory for output
	bool doPlots;			//< whether to make individual octet plots
	std::string simFile;	//< location of Geant4 simulation data to use
	float simFactor;		//< ratio of simulation to data events to produce
	unsigned int nTot;		//< total number of individual sim files
	unsigned int stride;	//< number of sim files to load in a chunk for each octet
};

#endif
