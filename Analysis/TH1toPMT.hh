#ifndef TH1TOPMT_HH
#define TH1TOPMT_HH 1

#include "ProcessedDataScanner.hh"
#include "PMTGenerator.hh"
#include "G4toPMT.hh"
#include "SectorCutter.hh"
#include <cassert>

/// base class for generating event positions
class PosGen {
public:
	/// constructor
	PosGen() { pos[X_DIRECTION]=pos[Y_DIRECTION]=pos[Z_DIRECTION]=0; }
	/// destructor
	virtual ~PosGen() {}
	/// generate next position
	virtual void next() {}
	
	float pos[3];	//< position
};

/// position generator from sector cutter
class SectPosGen: public PosGen {
public:
	/// constructor
	SectPosGen(const SectorCutter& S): PosGen(), sects(S) {}
	/// generate next position
	virtual void next() { sects.randPos(m,pos[X_DIRECTION],pos[Y_DIRECTION]); }
	
	SectorCutter sects;	//< position-generating sector cutter
	unsigned int m;		//< sector to generator positions for
};

/// supply event data from an input energy spectrum
class TH1toPMT: public Sim2PMT {
public:
	/// constructor
	TH1toPMT(TH1* h, PosGen* P);
	
	/// whether to count this event as successfully generated (all events count)
	virtual double simEvtCounts() const { return 1.0; }
	
	// ----- inapplicable disabled functions ----- //
	/// speedload: doesn't make sense for this class
	virtual void speedload(unsigned int e) { assert(false); }
	
	TH1* mySpectrum;		//< spectrum to throw events from
	PosGen* PG;				//< for generating event positions
	Side genside;			//< side to generate events on
	unsigned int nToSim;	//< total number of events to simulate (set to 0 for random energy selection)
	
protected:
	/// select energy for simulation
	virtual void doUnits();
};


#endif
