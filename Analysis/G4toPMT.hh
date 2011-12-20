#ifndef G4TOPMT_HH
#define G4TOPMT_HH 1

#include "PMTGenerator.hh"
#include "ProcessedDataScanner.hh"
#include <string>

/// generic class for converting simulation data to standard form
class Sim2PMT: public ProcessedDataScanner {
public:
	/// constructor
	Sim2PMT(const std::string& treeName);
	
	/// set calibrator to use for simulations
	void setCalibrator(PMTCalibrator& PCal);
	/// generate type 0 events with simulated detector response (with optional position cuts)
	void genType0(unsigned int nToSim = 0, double wx = 0, double wy = 0);
	
	/// this does nothing for processed data
	virtual void recalibrateEnergy() {}
	/// overrides ProcessedDataScanner::nextPoint to insert reverse-calibrations
	virtual bool nextPoint();
	
	/// get true energy
	virtual float getEtrue();
	
	/// check whether this is simulated data
	virtual bool isSimulated() const { return true; }
	
	/// return AFP state for data (note: may need to use physicsWeight for this to be meaningful)
	virtual AFPState getAFP() const { return afp; }
	/// set desired AFP state for simulation data
	virtual void setAFP(AFPState a) { afp=a; }
	
	PMTGenerator PGen[2];		//< PMT simulator for each side
	bool reSimulate;			//< whether to re-simulate energy or use "raw" values
	TH1F* inputEnergy[2];		//< histograms to fill with simulation input energy
	double eQ[2];				//< Scintillator quenched energy
	double eDep[2];				//< Scintillator deposited energy
	double eW[2];				//< Wirechamber deposited energy
	double scintPos[2][3];		//< hit position in scintillator
	double mwpcPos[2][3];		//< hit position in MWPC
	double primPos[4];			//< primary event vertex position (4=radius)
	double time[2];				//< hit time in each scintillator
	double costheta;			//< primary event cos pitch angle
	double ePrim;				//< primary event energy
	double physicsWeight;		//< event spectrum re-weighting factor
	
protected:
	/// perform unit conversions, etc.
	virtual void doUnits() { assert(false); }
	/// "reverse calibration" from simulated data
	virtual void reverseCalibrate();
	/// calculate spectrum re-weighting factor
	virtual void calcReweight();
	
	AFPState afp;				//< AFP state for data
};


/// converts Geant 4 simulation results to PMT spectra
class G4toPMT: public Sim2PMT {
public:
	/// constructor
	G4toPMT(): Sim2PMT("anaTree"), matchPenelope(false) { }
	/// unit conversions
	virtual void doUnits();
	
	bool matchPenelope;		//< whether to apply fudge factors to better match Penelope data (until this is fixed)
	
protected:
	virtual void setReadpoints();
};


/// converts Robbie's Penelope data to PMT spectra
class PenelopeToPMT: public Sim2PMT {
public:
	/// constructor
	PenelopeToPMT(): Sim2PMT("h34") { }	
	/// unit conversions
	virtual void doUnits();
	
	float fEprim;			//< float version for primary energy
	float fEdep[2];			//< float version for scintillator energy
	float fEW[2];			//< float version of wirechamber energy
	float fMWPCpos[2][2];	//< float version of MWPC position
	float fPrimPos[3];		//< float version of primary position
	float fTime[2];			//< float version of time
	float fCostheta;		//< float version of cos theta
	
protected:
	virtual void setReadpoints();
};

#endif
