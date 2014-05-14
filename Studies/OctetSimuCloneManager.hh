#ifndef OCTETSIMUCLONEMANAGER_HH
#define OCTETSIMUCLONEMANAGER_HH

#include "RunAccumulator.hh"
#include "PathUtils.hh"
#include <string>

/// class for handling octet-by-octet matching data and simulation analyses
class OctetSimuCloneManager {
public:
	/// constructor
	OctetSimuCloneManager(const std::string& dname, const std::string& bdir = getEnvSafe("UCNA_ANA_PLOTS"));
	/// destructor
	virtual ~OctetSimuCloneManager() { setSimData(NULL); }
	
	/// scan runs by provided run list
	void scanOct(RunAccumulator& RA, const Octet& oct);
	/// scan one octet of data, by octet number
	void scanOct(RunAccumulator& RA, unsigned int octn);
	/// combine all data octets
	void combineOcts(RunAccumulator& RA);
	
	/// set simulation data source
	void setSimData(Sim2PMT* s2p);
	/// simulate one octet
	void simOct(RunAccumulator& SimRA, const Octet& oct);
	/// simulate one octet, by octet number
	void simOct(RunAccumulator& SimRA, unsigned int octn);
	/// combine simulated octets; optionally, compare to data
	void combineSims(RunAccumulator& SimRA, RunAccumulator* OrigRA = NULL);
	
	/// run re-calculation on all octets
	unsigned int recalcAllOctets(RunAccumulator& RA, bool doPlots);
	
	std::string outputDir;	///< output directory naming
	std::string baseDir;	///< base directory for output
	bool doPlots;			///< whether to make individual octet plots
	bool doCompare;			///< whether to run data/MC comparison on cloned octets
	unsigned int hoursOld;	///< do not re-scan data for runs less than this number of hours old
	std::string simFile;	///< location of Geant4 simulation data to use
	float simFactor;		///< ratio of simulation to data events to produce
	unsigned int nTot;		///< total number of individual sim files
	unsigned int stride;	///< number of sim files to load in a chunk for each octet

protected:

	/// set up simulation data for specified octet
	virtual void setOctetSimdata(unsigned int octn);

	bool ownSimData;		///< whether this class ``owns'' simulation data
	Sim2PMT* simData;		///< simulated data source to use
};

#endif
