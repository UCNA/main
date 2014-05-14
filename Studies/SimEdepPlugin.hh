#ifndef SIMEDEPPLUGIN_HH
#define SIMEDEPPLUGIN_HH

#include "RunAccumulator.hh"
#include <TProfile.h>

const unsigned int NUM_EDEP_VOLS = 11;

/// plugin with extra plots for simulated runs
class SimEdepPlugin: public AnalyzerPlugin {
public:
	/// constructor
	SimEdepPlugin(RunAccumulator* RA);
	
	/// fill from scan data point
	virtual void fillCoreHists(ProcessedDataScanner& PDS, double weight);
	/// make combined eloss table
	virtual void calculateResults();
	/// plot results
	virtual void makePlots();
	
	/// make energy losses table
	void makeBigTable();
	
	std::string volNames[NUM_EDEP_VOLS];
	TProfile* hEdep[WEST+1][TYPE_IV_EVENT][NUM_EDEP_VOLS];
	TProfile* hEdepCombo[NUM_EDEP_VOLS];
};

class SimEdepAnalyzer: public RunAccumulator {
public:
	/// constructor
	SimEdepAnalyzer(OutputManager* pnt, const std::string& nm = "SimEdepAnalyzer", const std::string& inflName = "");
	/// create a new instance of this object (cloning self settings) for given directory
	virtual SegmentSaver* makeAnalyzer(const std::string& nm, const std::string& inflname) { return new SimEdepAnalyzer(this,nm,inflname); }

	SimEdepPlugin* myEdep;
};

#endif
