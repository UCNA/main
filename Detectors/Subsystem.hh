#ifndef SUBSYSTEM_HH
#define SUBSYSTEM_HH 1

#include "SpectrumPeak.hh"
#include "Types.hh"
#include "OutputManager.hh"
#include "RunManager.hh"

#include <TObject.h>
#include <TCanvas.h>
#include <cmath>
#include <vector>
#include <string>
#include <map>
#include "Enums.hh"

#include <TF1.h>
#include <TGraph.h>
#include <TProfile.h>

/// generic UCNA detector subsystem base class
class Subsystem: public OutputManager {
public:
	/// constructor
	Subsystem(RunManager* T, std::string nm, Side s);
	/// destructor
	virtual ~Subsystem() {
		if(triggers) delete[](triggers);
		for(std::vector<float*>::iterator it = edata.begin(); it != edata.end(); it++)
			delete(*it);
	}
	
	/// whether an event is a valid trigger for this subsystem
	virtual bool triggered(UInt_t e) const { return triggers[e]; }
	
	/// add appropriate branches to the output TTree
	virtual void addOutBranches(TTree*) { }
	/// load appropriate data for ouput TTree filling
	virtual void fillEvent(UInt_t) { }
	
	/// make histograms of bit distributions
	void bitsHisto(UInt_t n, bool* selector = NULL) const;
				
	Side mySide;		//< which side this detector lies on
		
	std::vector<std::string> sensorNames;	//< physical name of each sensor
	PedestalCorrector PC;					//< pedestal corrections for this run
	
	/// add the named branch to the cached event data
	void loadData(std::string branchName, std::string sensorName = "", std::string dataName = "");
		
	/// get named data for given event
	float* getData(std::string dataName, int offset = 0);
	/// get named data for given event, const version
	const float* getData(std::string dataName, int offset = 0) const;
	
	std::vector<float*> edata;					//< cached event data, by sensor
	std::map<std::string,unsigned int> nameMap;	//< map data names to data columns
	RunManager* thisRun;						//< RunManager for run being analyzed
	bool* triggers;								//< event triggers recognized by this detector
	
	const UInt_t nEvents;		//< number of events in TTree
		
	/// check whether pedestals need be calculated for ADC n
	bool verifyPedestal(int n);

	virtual void setName(std::string nm) { OutputManager::setName(nm); printf("\n--- %i: %s ---\n",thisRun->RI.runNum,name.c_str()); }
	
	/// continuous peak monitoring function
	void monitorPeak(std::vector< std::pair<float,float> > (*selector)(Subsystem*,void*),
						  float tmin, unsigned int cmin, void* selector_params, std::string mon_name, bool savePed = true);
};

#endif

