#ifndef PROCESSEDDATASCANNER_HH
#define PROCESSEDDATASCANNER_HH 1

#include "Types.hh"
#include "Enums.hh"

#include "TChainScanner.hh"
#include "EnergyCalibrator.hh"
#include "WirechamberReconstruction.hh"

#include "QFile.hh"
#include "PMTGenerator.hh"
#include "CalDBSQL.hh"
#include "TagCounter.hh"

#include <map>

/// Generic class for processed data TChains
class ProcessedDataScanner: public TChainScanner {
public:
	/// constructor
	ProcessedDataScanner(const std::string& treeName, bool withCalibrators = false);
	/// destructor
	virtual ~ProcessedDataScanner();
	
	/// add run to data
	virtual unsigned int addRun(RunNum rn);
	/// add list of runs to data; return number successfully added
	virtual unsigned int addRuns(const std::vector<RunNum>& rns);
	
	/// radius squared of event
	virtual float radius2(Side s) const;
	/// radius of event
	virtual float radius(Side s) const { return sqrt(radius2(s)); }
	/// event energy
	virtual float getEnergy() const { return scints[EAST].energy.x + scints[WEST].energy.x; }
	/// get event true (reconstructed) energy
	virtual float getEtrue();
	/// re-calibrate tube energy of currently loaded event
	virtual void recalibrateEnergy();
	/// speedload, keeping track of currently loaded run number
	virtual void speedload(unsigned int e);
	/// get run number of current event
	virtual RunNum getRun() const { return evtRun; }
	/// whether event passes fiducia/position cut on side
	virtual bool passesPositionCut(Side s);
	
	
	/// check whether this is simulated data
	virtual bool isSimulated() const { return false; }
	
	PMTCalibrator* ActiveCal;	//< PMTCalibrator currently active for loaded run
	
	/// print info about this scanner
	void display();
	/// write run calibrations info to QFile
	void writeCalInfo(QFile& qout, std::string tag);
	
	ScintEvent scints[2];		//< readout point for scintillator data
	Float_t led_pd[2];			//< readout point for reference photodiode
	wireHit wires[2][2];		//< readout point for wirechamber data [side][plane]
	MWPCevent mwpcs[2];			//< readout point for mwpc data (anode & cathode sum)
	Float_t mwpcEnergy[2];		//< calibrated wirechamber energy deposition on each side
	BlindTime runClock;			//< time of current event since run start
	RunNum evtRun;				//< run number for current event
	
	PID fPID;					//< analysis particle ID
	EventType fType;			//< analysis event type
	Side fSide;					//< analysis event side
	
	BlindTime totalTime;		//< combined length of runs in seconds	
	unsigned int nAFP[2];		//< number of events in each AFP state
	AnalysisChoice anChoice;	//< which analysis choice to use in identifying event types
	float fiducialRadius;		//< radius for position cut
	
	TagCounter<RunNum>	runTimes;	//< times for each run loaded
	
protected:
	
	/// generate event classification flags
	virtual void calcEventFlags() {}
	
	bool withCals;							//< whether to use energy recalibrators
	CalDB* CDB;								//< calibrations DB
	std::vector<RunNum> runlist;			//< list of loaded runs
	std::map<RunNum,PMTCalibrator*> PCals;	//< calibrators for each run
};

#endif
