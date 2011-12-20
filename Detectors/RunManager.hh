#ifndef RUNMANAGER_HH
#define RUNMANAGER_HH 1

#include <string>
#include <map>
#include <stdlib.h>
#include <cassert>

#include "RunInfo.hh"
#include "SourcesInfo.hh"
#include "EnergyCalibrator.hh"
#include "OutputManager.hh"

#include <TBranch.h>
#include <TTree.h>
#include <TFile.h>


class Subsystem;

// analysis mode flags
enum RunMode {
	RUNMODE_MINIMAL		= 0,	//< do bare minimum to not crash
	
	RUNMODE_MAKEPLOTS	= 1<<1,	//< make plots
	RUNMODE_PLOTSOUT	= 1<<2,	//< write plots to output file
	RUNMODE_TREEOUT		= 1<<3,	//< write full tree to output file
	RUNMODE_PEDESTALS	= 1<<4,	//< skip all but pedestal production
	RUNMODE_LEDTRACK	= 1<<5, //< track GMS LED
	RUNMODE_POSTRACK	= 1<<6,	//< track wirechamber positions
	RUNMODE_FULLID		= 1<<7, //< full event identification
	
	RUNMODE_FULLPLOTS	= RUNMODE_POSTRACK | RUNMODE_FULLID | RUNMODE_MAKEPLOTS | RUNMODE_PLOTSOUT,	//< full plots sequence output 
	RUNMODE_FULLOUT		= RUNMODE_FULLPLOTS | RUNMODE_TREEOUT					 					//< full condensed tree output
};

/// manages data TTree and run info for analyzing a UCNA data run
class RunManager: public OutputManager {
public:
	/// constructor
	RunManager(RunNum r, RunMode RM = RUNMODE_FULLPLOTS, std::string bp = "../", bool noload = false);
	/// destructor
	virtual ~RunManager() {
		clearItems();
		if(myTree) delete(myTree);
		if(treeFile) { treeFile->Close(); delete(treeFile); }
		delete(PCal); 
	}
	/// check existence of named branch'
	bool checkBranch(std::string s) const;
	/// write rundata file
	virtual void write(std::string outName = "");
	
	RunMode runMode;				//< analysis running mode
	bool manualAnalyze;				//< true for runs not in Cal DB, etc. to force manual analysis
	UInt_t nEvents;					//< number of events in run
	RunInfo RI;						//< info about this run
	SourcesInfo SI;					//< info about sources in this run
	
	/// get path to data for run
	static std::string getDataPath(RunNum rn);
	
	/// read all data from named branch, return as array
	float* readBranch(std::string bname);
	
	/// attach a subsystem to this run analysis
	void attachSubsystem(Subsystem* S) {
		attach((OutputManager*)S);
		subsystems.push_back(S);
	}
	
	/// write output tree
	void makeRootOut(Subsystem* trigger = NULL);

	PMTCalibrator* PCal;
	CalDB* CDB;
	
private:
	TTree* myTree;									//< data TTree for this run
	TFile* treeFile;								//< file containing data TTree
	
	std::vector<Subsystem*> subsystems;				//< subsystems in detector, for adding ROOT output branches
};


#endif

