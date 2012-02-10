#ifndef G4TOPMT_HH
#define G4TOPMT_HH 1

#include "PMTGenerator.hh"
#include "ProcessedDataScanner.hh"
#include <string>
#include <TRandom3.h>


/// generic class for converting simulation data to standard form
class Sim2PMT: public ProcessedDataScanner {
public:
	/// constructor
	Sim2PMT(const std::string& treeName);
	
	/// set calibrator to use for simulations
	virtual void setCalibrator(PMTCalibrator& PCal);
	
	/// this does nothing for processed data
	virtual void recalibrateEnergy() {}
	/// start scanning events
	virtual void startScan(unsigned int startRandom = 0) { nSimmed = 0; ProcessedDataScanner::startScan(startRandom); }
	/// overrides ProcessedDataScanner::nextPoint to insert reverse-calibrations
	virtual bool nextPoint();
	/// whether to count this event as successfully generated
	virtual double simEvtCounts() const { return fPID==PID_BETA && fType==TYPE_0_EVENT?physicsWeight:0; }
	/// get event info
	virtual Stringmap evtInfo();
	
	/// get true energy
	virtual float getEtrue();
	
	/// check whether this is simulated data
	virtual bool isSimulated() const { return true; }
	
	/// return AFP state for data (note: may need to use physicsWeight for this to be meaningful)
	virtual AFPState getAFP() const { return afp; }
	/// set desired AFP state for simulation data
	virtual void setAFP(AFPState a) { afp=a; }
	
	/// Determine event classification flags
	virtual void classifyEvent();
	/// calculate spectrum re-weighting factor
	virtual void calcReweight();
	
	PMTGenerator PGen[2];		//< PMT simulator for each side
	bool reSimulate;			//< whether to re-simulate energy or use "raw" values
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
	unsigned int nToSim;		//< total number of events to simulate (set to 0 for random energy selection)
	unsigned int nSimmed;		//< number of events simulated since scan start
	double mwpcThresh[2];		//< MWPC trigger threshold on each side
	
protected:
	/// perform unit conversions, etc.
	virtual void doUnits() { assert(false); }
	/// "reverse calibration" from simulated data
	virtual void reverseCalibrate();
	
	AFPState afp;				//< AFP state for data
	bool passesScint[2];		//< whether simulation passed scintillator cut
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


/// mixes several simulations
class MixSim: public Sim2PMT {
public:
	/// constructor
	MixSim(double tinit=0): Sim2PMT(""), currentSim(NULL), t0(tinit), t1(tinit) {}
	
	/// start scan
	virtual void startScan(unsigned int startRandom = 0);
	/// load sub-simulation
	virtual bool nextPoint();
	
	/// add sub-simulation
	void addSim(Sim2PMT* S, double r0, double thalf);
	/// set simulation time (determines different line strengths)
	void setTime(double t);
	
	/// set desired AFP state for simulation data
	virtual void setAFP(AFPState a);
	/// set calibrator to use for simulations
	virtual void setCalibrator(PMTCalibrator& PCal);
	
protected:
	
	virtual void doUnits() { }
	
	std::vector<Sim2PMT*> subSims;
	std::vector<double> initStrength;
	std::vector<double> halflife;
	std::vector<double> cumStrength;
	Sim2PMT* currentSim;
	double t0;	//< initial time
	double t1;	//< current time
	
};

#endif
