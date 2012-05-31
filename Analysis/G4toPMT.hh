#ifndef G4TOPMT_HH
#define G4TOPMT_HH 1

#include "PMTGenerator.hh"
#include "ProcessedDataScanner.hh"
#include <string>
#include <TRandom3.h>

class Sim2PMT;

/// base class for setting simulated data position
class SimPositioner {
public:
	/// constructor
	SimPositioner() { offPos[X_DIRECTION] = offPos[Y_DIRECTION] = 0; }
	/// destructor
	virtual ~SimPositioner() {}
	/// apply position offset to event
	virtual void applyOffset(Sim2PMT& S);
	
protected:	
	/// calculate offset based on primary position
	virtual void calcOffset(const Sim2PMT& S) {}
	
	double offPos[2];		//< position offset to apply
};

/// set positions for source droplet
class SourcedropPositioner: public SimPositioner {
public:
	/// constructor
	SourcedropPositioner(double x, double y, double r): SimPositioner(), x0(x), y0(y), r0(r) {}
	/// calculate offset based on primary position
	virtual void calcOffset(const Sim2PMT& S);
	
	double x0;	//< x position center
	double y0;	//< y position center
	double r0;	//< spot radius
};


/// generic class for converting simulation data to standard form
class Sim2PMT: public ProcessedDataScanner {
public:
	/// constructor
	Sim2PMT(const std::string& treeName);
	
	/// set calibrator to use for simulations
	virtual void setCalibrator(PMTCalibrator& PCal);
	
	/// this does nothing for processed data
	virtual void recalibrateEnergy() {}
	/// overrides ProcessedDataScanner::nextPoint to insert reverse-calibrations, offsets
	virtual bool nextPoint();
	/// whether to count this event as successfully generated
	virtual double simEvtCounts() const { return fPID==PID_BETA && fType==TYPE_0_EVENT?physicsWeight:0; }
	/// reset simulation counters
	virtual void resetSimCounters() { nSimmed = nCounted = 0; }
	/// get event info
	virtual Stringmap evtInfo();
	
	/// get true energy
	virtual float getEtrue();
	/// primary event radius
	virtual float primRadius() const { return sqrt(pow(primPos[X_DIRECTION],2)+pow(primPos[Y_DIRECTION],2)); }
	
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
	SimPositioner* SP;			//< optional postion modifier
	bool reSimulate;			//< whether to re-simulate energy or use "raw" values
	bool fakeClip;				//< whether to fake clipping on wirechamber entrance edge
	double eQ[2];				//< Scintillator quenched energy
	double eDep[2];				//< Scintillator deposited energy
	double eW[2];				//< Wirechamber deposited energy
	double scintPos[2][3];		//< hit position in scintillator
	double mwpcPos[2][3];		//< hit position in MWPC
	double primPos[4];			//< primary event vertex position (4=radius)
	double time[2];				//< hit time in each scintillator
	double costheta;			//< primary event cos pitch angle
	double ePrim;				//< primary event energy
	unsigned int nSimmed;		//< number of events simulated since scan start
	double nCounted;			//< physics-weighted number of counted events
	double mwpcThresh[2];		//< MWPC trigger threshold on each side
	double mwpcAccidentalProb;	//< probability of MWPC accidental triggers
	
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
	G4toPMT(): Sim2PMT("anaTree") { }
	/// unit conversions
	virtual void doUnits();
		
protected:
	/// set read points for input tree
	virtual void setReadpoints();
};

/// For consistency checks, swaps E/W sides on Geant4 sim data
class G4toPMT_SideSwap: public G4toPMT {
public:
	/// constructor
	G4toPMT_SideSwap(): G4toPMT() { }
	/// unit conversions
	virtual void doUnits();
};



/// converts Robby's Penelope data to PMT spectra
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

/// multiply Geant4 data to fill all segments of SectionCutter
class G4SegmentMultiplier: public G4toPMT, SimPositioner {
public:
	/// constructor
	G4SegmentMultiplier(const SectorCutter& S);
	/// start scan
	virtual void startScan(bool startRandom = false);
	/// overrides G4toPMT::nextPoint to re-use points
	virtual bool nextPoint();
	/// unit conversions, only done once per rotation sequence
	virtual void doUnits();
	/// calculate spectrum re-weighting factor
	virtual void calcReweight();
	/// apply rotational offset
	virtual void applyOffset(Sim2PMT& S);
protected:
	/// rotate a point
	void rotpt(double& x0, double& y0);
	
	SectorCutter SC;			//< SectorCutter to determine event multiplication
	unsigned int nrots;			//< number of remaining point rotations
	unsigned int rcurrent;		//< current radial ring
	std::vector<double> vc,vs;	//< pre-calculated rotation matrix cosines, sines for each ring
	bool morePts;				//< whether the data still has more points to come
};

/// mixes several simulations
class MixSim: public Sim2PMT {
public:
	/// constructor
	MixSim(double tinit=0): Sim2PMT(""), currentSim(NULL), t0(tinit), t1(tinit) {}
	
	/// start scan
	virtual void startScan(bool startRandom = false);
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
	/// get number of files loaded by sub simulations
	virtual unsigned int getnFiles() const;
	
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
