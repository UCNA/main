#ifndef TH1TOPMT_HH
#define TH1TOPMT_HH 1

#include "ProcessedDataScanner.hh"
#include "PMTGenerator.hh"
#include <cassert>

// supply event data from an input energy spectrum
class TH1toPMT: public ProcessedDataScanner {
public:
	/// constructor
	TH1toPMT(TH1* h);
	
	// ----- inapplicable disabled functions ----- //
	/// add run to data -- NO
	virtual unsigned int addRun(RunNum rn) { assert(false); }
	/// add list of runs to data; return number successfully added -- NO
	virtual unsigned int addRuns(const std::vector<RunNum>& rns) { assert(false); }
	/// speedload, keeping track of currently loaded run number
	virtual void speedload(unsigned int e) { assert(false); }
	/// get run number of current event
	virtual RunNum getRun() const { return 0; }
	
	// ------ hijack these for our use ------ //
	/// load next "speed scan" point
	virtual bool nextPoint();	
	/// start a "speed scan," possibly at a random entry number
	void startScan(unsigned int startRandom = 0);
	
	/// set calibrator to use for simulations
	void setCalibrator(PMTCalibrator& PCal);
	/// set event generation postion
	void setPosition(float x, float y);
	
	TH1* mySpectrum;			//< spectrum to throw events from
	bool stochasticEnergy;		//< whether to select energies randomly or deterministically from spectrum
	float randomPositionRadius;	//< random event positioning radius (set <0 for fixed position)
	
protected:
	unsigned int nToSim;	//< total number of events to simulate
	unsigned int nSimmed;	//< number of events simulated so far
	PMTGenerator PGen[2];	//< PMT simulator for each side
};


#endif
