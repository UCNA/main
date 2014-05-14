#ifndef SIMCALIBRATIONS_HH
#define SIMCALIBRATIONS_HH

#include "EnergyCalibrator.hh"
#include "Types.hh"
#include "EfficCurve.hh"
#include "SMExcept.hh"
#include <vector>
#include <TRandom3.h>
#include <TMultiLayerPerceptron.h>

/// Base class for estimating 2-fold trigger probability
class TriggerProb {
public:
	/// constructor
	TriggerProb() {}
	/// destructor
	virtual ~TriggerProb() {};
	/// calculate 2-fold trigger probability
	virtual double calcProb();
	
	Double_t tubeProbs[(nBetaTubes*(nBetaTubes+1))/2];	///< input location for individual PMT probabilities
	ScintEvent sevt;									///< scintillator event to calculate trigger probability from
};

/// Trigger probability using TMultiLayerPerceptron
class TriggerProbMLP: public TriggerProb {
public:
	/// constructor
	TriggerProbMLP(TMultiLayerPerceptron* M): TriggerProb(), TMLP(M) { smassert(TMLP); }
	/// calculate 2-fold trigger probability
	virtual double calcProb();
	
	/// pre-modify input array
	static void condition(Double_t* aIn);
	
protected:
	TMultiLayerPerceptron* TMLP;
	static double unfold(double x);
};

/// Class for simulating detector response
class PMTGenerator {
public:
	/// constructor
	PMTGenerator(Side s = EAST, float xx = 0, float yy = 0);
	/// destructor
	~PMTGenerator() {}
	
	/// load a PMTCalibrator for event generation
	void setCalibrator(PMTCalibrator* P);
	/// set trigger probability estimator
	void setTriggerProb(TriggerProb* TP);
	
	/// generate a scintillator event for a given quenched energy
	ScintEvent generate(float en);
	
	/// calculate and count number of PMT triggers for event
	unsigned int triggers();
	/// whether the 2-of-4 trigger fired
	bool triggered();
	
	/// set generator position (sets scint and wirechamber position)
	void setPosition(float xx, float yy, float dxw=0, float dyw=0);
	/// set event side
	void setSide(Side s);
	/// set custom light balance
	void setLightbal(Side s, float l1, float l2, float l3, float l4);
		
	/// get current calibrator
	const PMTCalibrator* getCalibrator() const { return currentCal; }
		
	float x,y;						///< event hit position in scintillator (projected back to decay trap)
	float xw,yw;					///< wirechamber hit position, possibly offset from source position
	float evtm;						///< event time during run
	float presmear;					///< nPE/keV already smeared in input spectrum
	float dgain;					///< gain at first photomultiplier stage (smears single-PE resolution)
	float pedcorr;					///< pedestal noise correlation
	float crosstalk;				///< inter-channel noise crosstalk
	float xscatter;					///< "extra" proportional random scatter
	float trigThreshScale;			///< ADC scaling factor for trigger efficiency (to test efficiency changes)
	
	unsigned int nTrigs;			///< number of individual PMTs triggered
	bool pmtTriggered[nBetaTubes];	///< whether each PMT triggered above threshold
	
	static TRandom3 sim_rnd_source;	///< PMTGenerator random number generator
	
	
	/// convert simulated cathode charge distribution to ADC, updating wireHit results
	void calcCathodeSignals(Side s, AxisDirection d, const float* cath_chg, float* cath_adc, wireHit& w) const;
	
protected:

	PMTCalibrator* currentCal;			///< current PMT Calibrator in use
	TriggerProb* TProb;
	Side mySide;						///< side to simulate
	float pmtRes[BOTH][nBetaTubes];		///< individual PMT nPE per keV
	float lightBal[BOTH][nBetaTubes];	///< custom light balancing factor between PMTs for LED events
	ScintEvent sevt;					///< current generated event
};

#endif
