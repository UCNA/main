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
	void setPosition(float xx, float yy);
	/// set scint/wirechamber offsets
	void setOffsets(float xxs, float yys, float xxw, float yyw);
	/// set event side
	void setSide(Side s) { mySide = s; }
	
	/// get current calibrator
	const PMTCalibrator* getCalibrator() const { return currentCal; }
	
	bool calcADC;					//< whether to calculate PMT ADCs and trigger efficiency
	
	float x,y;						//< source position in scintillator (projected back to decay trap)
	float dsx,dsy;					//< scintillator hit offset from source position
	float dwx,dwy;					//< wirechamber hit offset from source position
	float presmear;					//< nPE/keV already smeared in input spectrum
	float larmorField;				//< optional field spreading out hit positions (0 for no effect, and use ds*, dw* for offsets)
	
	unsigned int nTrigs;			//< number of individual PMTs triggered
	bool pmtTriggered[nBetaTubes];	//< whether each PMT triggered above threshold
	
protected:

	PMTCalibrator* currentCal;		//< current PMT Calibrator in use
	
	Side mySide;					//< side to simulate
	float pmtRes[2][nBetaTubes];	//< individual PMT nPE per keV
	float resTot[2];				//< total resolution nPE per keV for all 4 tubes
	ScintEvent sevt;				//< current generated event
};

#endif
