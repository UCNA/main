#ifndef SIMCALIBRATIONS_HH
#define SIMCALIBRATIONS_HH 1

#include "EnergyCalibrator.hh"
#include "Types.hh"
#include "EfficCurve.hh"
#include <vector>

/// Class for generating PMT signals with energy resolution, efficiency considerations
class PMTGenerator {
public:
	/// constructor
	PMTGenerator(Side s = EAST, float xx = 0, float yy = 0);
	/// destructor
	~PMTGenerator() {}
	
	/// load a PMTCalibrator for event generation
	void setCalibrator(PMTCalibrator* P);
	
	/// generate an event for a given quenched energy
	ScintEvent generate(float en);
	
	/// calculate and count number of PMT triggers for event
	unsigned int triggers();
	/// whether the 2-of-4 trigger fired
	bool triggered() { return triggers() >= 2; }
	
	/// set generator position (sets scint and wirechamber position)
	void setPosition(float xx, float yy, float dxw=0, float dyw=0);
	/// set event side
	void setSide(Side s) { mySide = s; }
	
	/// get current calibrator
	const PMTCalibrator* getCalibrator() const { return currentCal; }
		
	float x,y;						//< event hit position in scintillator (projected back to decay trap)
	float xw,yw;					//< wirechamber hit offset from source position
	float presmear;					//< nPE/keV already smeared in input spectrum
	float threshnoise;				//< electronic threshold noise (equivalent photoelectrons)
	float propnoise;				//< extra fraction of noise proportional to signal
	
	unsigned int nTrigs;			//< number of individual PMTs triggered
	bool pmtTriggered[nBetaTubes];	//< whether each PMT triggered above threshold
	
protected:

	PMTCalibrator* currentCal;		//< current PMT Calibrator in use
	
	Side mySide;					//< side to simulate
	float pmtRes[2][nBetaTubes];	//< individual PMT nPE per keV
	ScintEvent sevt;				//< current generated event
};

#endif
