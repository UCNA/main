#ifndef SIM2PMT_HH
#define SIM2PMT_HH 1

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
	
	/// set whether to simulate cathodes response --- implement in subclass
	virtual void runCathodeSim(bool = true) {}
	/// this does nothing for processed data
	virtual void recalibrateEnergy() {}
	/// overrides ProcessedDataScanner::startScan to clear simulation counters
	virtual void startScan(bool startRandom = false);
	/// overrides ProcessedDataScanner::nextPoint to insert reverse-calibrations, offsets
	virtual bool nextPoint();
	/// whether to count this event as successfully generated
	virtual double simEvtCounts() const { return fPID==PID_BETA && fType==TYPE_0_EVENT?physicsWeight:0; }
	/// reset simulation counters
	virtual void resetSimCounters() { nSimmed = nCounted = 0; }
	/// get event info
	virtual Stringmap evtInfo();
	
	/// get true energy
	virtual float getErecon() const;
	/// primary event radius
	virtual float primRadius() const { return sqrt(pow(primPos[X_DIRECTION],2)+pow(primPos[Y_DIRECTION],2)); }
	
	/// check whether this is simulated data
	virtual bool isSimulated() const { return true; }
	/// whether event was simulated as triggering the given side
	virtual bool Sis00_2fold(Side s) { return PGen[s].triggered(); }
	
	/// return AFP state for data (note: may need to use physicsWeight for this to be meaningful)
	virtual AFPState getAFP() const { return afp; }
	/// set desired AFP state for simulation data
	virtual void setAFP(AFPState a) { afp=a; }
	
	/// Determine event classification flags
	virtual void classifyEvent();
	/// calculate spectrum re-weighting factor
	virtual void calcReweight();
	
	PMTGenerator PGen[BOTH];		//< PMT simulator for each side
	SimPositioner* SP;				//< optional postion modifier
	bool reSimulate;				//< whether to re-simulate energy or use "raw" values
	bool fakeClip;					//< whether to fake clipping on wirechamber entrance edge
	bool weightAsym;				//< whether to weight simulated events by beta asymmetry
	
	double eQ[BOTH];						//< Scintillator quenched energy [keV]
	double eDep[BOTH];						//< Scintillator deposited energy [keV]
	double eW[BOTH];						//< Wirechamber active volume deposited energy [keV]
	double scintPos[BOTH][Z_DIRECTION+1];	//< hit position in scintillator [mm in decay trap]
	double mwpcPos[BOTH][Z_DIRECTION+1];	//< hit position in MWPC
	double primPos[Z_DIRECTION+2];			//< primary event vertex position (4=radius)
	double time[BOTH];						//< hit time [s] in each scintillator
	double costheta;						//< primary event cos pitch angle
	Side primSide;							//< side primary event is heading towards
	double ePrim;							//< primary event energy
	
	double edepFoils[BOTH];			//< energy deposition [keV] in decay trap foils
	double edepWinOut[BOTH];		//< energy deposition [keV] in outer wirechamber window
	double edepWinIn[BOTH];			//< energy deposition [keV] in inner wirechamber window
	double edepDeadMWPC[BOTH];		//< energy deposition [keV] in MWPC dead volume
	double edepKevlar[BOTH];		//< energy deposition [keV] in kevlar strings
	double edepWires[BOTH];			//< energy deposition [keV] in wire planes
	double edepDeadScint[BOTH];		//< energy deposition [keV] in dead scintillator
	float cath_chg[BOTH][Y_DIRECTION+1][kMaxCathodes];	//< signal on each cathode segment (portion of eW)
	
	double cosThetaInFoils[BOTH];	//< entrance angle cosine to decay trap foils
	double cosThetaInWinOut[BOTH];	//< entrance angle cosine to outer wirechamber window
	double cosThetaInWinIn[BOTH];	//< entrance angle cosine to inner wirechamber window
	double cosThetaInScint[BOTH];	//< entrance angle cosine to scintillator
	
	double cosThetaOutFoils[BOTH];	//< exit angle cosine from decay trap foils
	double cosThetaOutWinOut[BOTH];	//< exit angle cosine from outer wirechamber windo
	double cosThetaOutWinIn[BOTH];	//< exit angle cosine from inner wirechamber window
	double cosThetaOutScint[BOTH];	//< exit angle cosine from scintillator
	
	unsigned int nSimmed;			//< count of number of events simulated
	unsigned int nToSim;			//< number of events to simulate
	double nCounted;				//< physics-weighted number of counted events
	double mwpcThresh[BOTH];		//< MWPC trigger 50% threshold on each side
	double mwpcWidth[BOTH];			//< MWPC threshold width
	
protected:
	/// perform unit conversions, etc.
	virtual void doUnits() { assert(false); }
	/// generate event time stamp
	virtual void updateClock();
	/// "reverse calibration" from simulated data
	virtual void reverseCalibrate();
	
	AFPState afp;				//< AFP state for data
	bool passesScint[BOTH];		//< whether simulation passed scintillator cut
	bool simCathodes;			//< whether to simulate cathode response
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
