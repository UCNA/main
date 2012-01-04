#ifndef TH1TOPMT_HH
#define TH1TOPMT_HH 1

#include "ProcessedDataScanner.hh"
#include "PMTGenerator.hh"
#include "G4toPMT.hh"
#include <cassert>

// supply event data from an input energy spectrum
class TH1toPMT: public Sim2PMT {
public:
	/// constructor
	TH1toPMT(TH1* h);
	
	// ----- inapplicable disabled functions ----- //
	/// add run to data -- NO
	virtual unsigned int addRun(RunNum rn) { assert(false); }
	/// add list of runs to data; return number successfully added -- NO
	virtual unsigned int addRuns(const std::vector<RunNum>& rns) { assert(false); }
	/// speedload: doesn't make sense for this class
	virtual void speedload(unsigned int e) { assert(false); }
	
	TH1* mySpectrum;			//< spectrum to throw events from
	float randomPositionRadius;	//< random event positioning radius (set <0 for fixed position)
	float genpos[2];			//< position of simulated events
	Side genside;				//< side to generate events on
	unsigned int nToSim;		//< total number of events to simulate (set to 0 for random energy selection)
	
protected:
	/// select energy for simulation
	virtual void doUnits();
};


#endif
