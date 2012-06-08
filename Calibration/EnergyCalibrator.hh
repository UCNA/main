#ifndef ENERGYCALIBRATOR_HH
#define ENERGYCALIBRATOR_HH 1

#include "Types.hh"
#include "CalDB.hh"
#include "CalDBSQL.hh"
#include "GainStabilizer.hh"
#include "EvisConverter.hh"
#include "WirechamberCalibrator.hh"
#include "QFile.hh"
#include <map>
#include <string>
#include <vector>

/// pedestal subtraction class
class PedestalCorrector {
public:
	/// constructor
	PedestalCorrector(RunNum rn, CalDB* cdb): myRun(rn), totalTime(cdb->totalTime(rn)), pCDB(cdb) {}
	/// destructor
	virtual ~PedestalCorrector();
	/// check whether pedestal info exists
	bool checkPedestals(const std::string& sensorName);
	/// get sensor pedestal at given time from beginning of run
	float getPedestal(const std::string& sensorName, float time);
	/// get sensor pedestal at given time from beginning of run
	float getPedwidth(const std::string& sensorName, float time);
	
	RunNum myRun;			//< run number for this run
	float totalTime;		//< run time for this run
	
	/// insert a pedestal graph
	void insertPedestal(const std::string& sensorName, TGraph* g);

private:
	std::map<std::string,TGraph*> pedestals;	//< pedestals history for each sensor
	std::map<std::string,TGraph*> pedwidths;	//< pedestal width history for each sensor
	CalDB* pCDB;								//< pedestal-containing DB
};


/// GMS and Linearity Corrections class
class LinearityCorrector {
public:
	/// constructor
	LinearityCorrector(RunNum myRun, CalDB* cdb);
	/// destructor
	virtual ~LinearityCorrector();
	/// get GMS correction factor
	float gmsFactor(Side s, unsigned int t, float time) const;
	/// get energy resolution factor
	float getDeltaL(Side s, unsigned int t) const { return deltaL[s][t]; }
	/// get energy resolution factor
	float getDeltaADC(Side s, unsigned int t) const { return deltaADC[s][t]; }
	/// positioning intensity factor eta
	virtual float eta(Side s, unsigned int t, float x, float y) const { return P->eval(s,t,x,y,true); }
	/// linearize tube adc (plus GMS correction), ADC -> L = eta*Evis
	float linearityCorrector(Side s, unsigned int t, float adc, float time) const;
	/// invert linearity correction to raw ADC
	float invertLinearity(Side s, unsigned int t, float l, float time) const;
	/// invert eta*Evis -> gain stabilized ADC
	float invertLinearityStabilized(Side s, unsigned int t, float l) const;
	/// derivative of inverse linearity
	float dInverse(Side s, unsigned int t, float l, float time) const;	
	/// invert linearity with float_err
	float_err invertLinearity(Side s, unsigned int t, float_err l, float time) const;
	/// linearity corrector derivative at given adc value
	float dLinearity(Side s, unsigned int t, float adc, float time) const;	
	/// get ref run t0 GMS factor
	float getGMS0(Side s, unsigned int t) const { return gms0[s][t]; }
	
	/// whether this is a reference run
	bool isRefRun() const { return rGMS == rn || !rGMS; }
		
	bool scaleNoiseWithL;	//< whether to scale noise with sqrt(Light) or sqrt(ADC)
	
	PositioningCorrector* P;		//< endpoint map for this run
	GainStabilizer* GS;				//< PMT gain stabilizer
	
	std::string sensorNames[2][nBetaTubes];	//< names for the PMT sensors
	RunNum rn;						//< run number for this run
	RunNum rGMS;					//< run number this is GMS corrected to
	LinearityCorrector* LCRef;		//< linearity corrector for GMS reference run
	CalDB* CDB;						//< Calibration DB link
	
protected:

	float deltaL[2][nBetaTubes];		//< energy resolution for each tube delta(eta*E)/sqrt(eta*E)
	float deltaADC[2][nBetaTubes];		//< energy resolution for each tube delta(ADC)/sqrt(ADC)
	float gms0[2][nBetaTubes];			//< GMS factor at reference run t=0
	float expected_adc[2][nBetaTubes];	//< expected ADC value for calibration peak
	
	TGraph* linearityFunctions[2][nBetaTubes];	//< linearity correction for each side, tube (including ref. pmt)
	TGraph* linearityInverses[2][nBetaTubes];	//< inverse linearity correction for each side, tube
	
	static std::map<RunNum,LinearityCorrector*> cachedRuns;			//< cache of run correctors for faster access
	static LinearityCorrector* getCachedRun(RunNum r,CalDB* cdb);	//< retrieve a cached corrector, creating if necessary
};

/// Energy reconstruction class
class PMTCalibrator: public LinearityCorrector, public PedestalCorrector, public EvisConverter, public WirechamberCalibrator {
public:
	/// Constructor, for specified run number
	PMTCalibrator(RunNum rn, CalDB* cdb = CalDBSQL::getCDB());
	/// Destructor
	~PMTCalibrator();
	/// resolution deltaLight for given light seen in tube
	float lightResolution(Side s, unsigned int t, float l, float time) const;
	/// expected width in raw ADC counts at given adc value
	float adcResolution(Side s, unsigned int t, float adc, float time) const;	
	/// expected energy resolution at given position, visible energy
	float energyResolution(Side s, unsigned int t, float e0, float x, float y, float time) const;	
	/// combined energy resolution of all 4 tubes
	float combinedResolution(Side s, float e0, float x, float y, float time = 0) const;	
	/// expected number of photoelectrons for an event
	float nPE(Side s, unsigned int t, float e0, float x, float y, float time = 0) const;	
	/// expected number of photoelectrons for an ADC value
	float nPE(Side s, unsigned int t, float adc, float time = 0) const;
	/// effective nPE for PMT sum
	float pmtSumPE(Side s, float e0, float x, float y, float time = 0) const;
	/// effective combined 4-PMT `eta'
	float combEta(Side s, float x, float y, float e0 = 500, float time = 0) const { return nPE(s,nBetaTubes,e0,x,y,time)/nPE(s,nBetaTubes,e0,0,0,time); }
	/// positioning intensity factor eta
	virtual float eta(Side s, unsigned int t, float x, float y) const { return (t<nBetaTubes)?LinearityCorrector::eta(s,t,x,y):combEta(s,x,y); }
	/// trigger efficiency at given ADC value for side tube
	float trigEff(Side s, unsigned int t, float adc) const;
	/// convert a single PMT ADC readout to an energy with an error estimate
	float_err calibratedEnergy(Side s, unsigned int t, float x, float y, float adc, float time = 0) const;	
	/// invert corrections to find raw ADC channel corresponding to specified energy
	float invertCorrections(Side s, unsigned int t, float e0, float x, float y, float time) const;
	/// invert corrections with error
	float_err invertCorrections(Side s, unsigned int t, float_err e0, float x, float y, float time) const;
	/// subtract pedestals from 4-tube ADC array (do this before energy calibrations!)
	void pedSubtract(Side s, float* adc, float time);
	/// convert all 4 tubes ped-subtracted ADC to energy estimates
	void calibrateEnergy(Side s, float x, float y, ScintEvent& evt, float time, float poserr = 0) const;	
	/// convert all 4 tubes ped-subtracted ADC to energy estimates and return QADC-sum averaged energy (ignores individual tube resolutions)
	void summedEnergy(Side s, float x, float y, ScintEvent& evt, float time) const;		
	/// print summary of energy calibrations
	virtual void printSummary();
	/// stringmap energy calibrations summary
	Stringmap calSummary() const;
	/// get clipping threshold for a PMT
	float getClipThreshold(Side s, unsigned int t) { assert(s<=WEST && t<nBetaTubes); return clipThreshold[s][t]; }
	
protected:
	
	float clipThreshold[2][nBetaTubes];		//< threshold to de-weight ADC in tube combination due to "clipping"
	EfficCurve* pmtEffic[2][nBetaTubes];	//< efficiency curves for each PMT
};

#endif
