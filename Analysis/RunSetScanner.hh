#ifndef RUNSETSCANNER_HH
#define RUNSETSCANNER_HH 1

#include "Enums.hh"
#include "TChainScanner.hh"
#include "EnergyCalibrator.hh"

#include "QFile.hh"
#include "CalDBSQL.hh"
#include "TagCounter.hh"

#include <map>

/// Class for scanning through files specified by run number
class RunSetScanner: public TChainScanner {
public:
	/// constructor
	RunSetScanner(const std::string& treeName, bool withCalibrators = false);
	/// destructor
	virtual ~RunSetScanner();
	
	/// add run to data (return whether successful)
	virtual bool addRun(RunNum rn);
	/// add list of runs to data; return number successfully added
	unsigned int addRuns(const std::vector<RunNum>& rns);
	/// return path to run .root file
	virtual std::string locateRun(RunNum r) { assert(false); return ""; }
	
	/// speedload, keeping track of currently loaded run number
	virtual void speedload(unsigned int e);
	/// subclass this for routines when new run is loaded
	virtual void loadNewRun(RunNum rn) {}
	/// get run number of current event
	virtual RunNum getRun() const { return evtRun; }
	/// check whether this is simulated data
	virtual bool isSimulated() const { return false; }
	
	PMTCalibrator* ActiveCal;	//< PMTCalibrator currently active for loaded run
	
	/// print info about this scanner
	virtual void display();
	/// write run calibrations info to QFile
	void writeCalInfo(QFile& qout, std::string tag);
	
	RunNum evtRun;					//< run number for current event
	BlindTime totalTime;			//< combined length of runs in seconds	
	unsigned int nAFP[2];			//< number of events in each AFP state
	TagCounter<RunNum>	runTimes;	//< times for each run loaded
	bool withCals;					//< whether to use energy recalibrators
	
protected:
	
	std::vector<RunNum> runlist;			//< list of loaded runs
	std::map<RunNum,PMTCalibrator*> PCals;	//< calibrators for each run
};

#endif
