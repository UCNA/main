#ifndef GAINSTABILIZER_HH
#define GAINSTABILIZER_HH 1

#include "CalDB.hh"
class LinearityCorrector;

/// generic gain stabilization class for matching PMT gain to reference run start
class GainStabilizer {
public:
	/// constructor
	GainStabilizer(RunNum myRun, CalDB* cdb, LinearityCorrector* myCorrector): CDB(cdb), LCor(myCorrector), rn(myRun) {}
	/// destructor
	virtual ~GainStabilizer() {}
	/// compute gain stabilization factor for given side, PMT
	virtual float gmsFactor(Side s, unsigned int t, float time) const;
	/// print gain stabilization info
	virtual void printSummary() { printf("Null Gain Stabilization\n"); }
	/// get a summary of GMS calibration parameters
	virtual Stringmap gmsSummary() const;
	
	CalDB* CDB;					//< reference to calibration DB
	LinearityCorrector* LCor;	//< reference to linearity corrector
	RunNum rn;					//< run num of this run
};

/// Chris Pulser gain stabilizer
class ChrisGainStabilizer: public GainStabilizer {
public:
	/// constructor
	ChrisGainStabilizer(RunNum myRun, CalDB* cdb, LinearityCorrector* myCorrecter);
	/// compute gain stabilization factor for given side, PMT
	virtual float gmsFactor(Side s, unsigned int t, float time) const;
	/// print gain stabilization info
	virtual void printSummary();
protected:
	TGraph* pulserPeak[2][nBetaTubes];			//< Chris Pulser peak position
	float pulser0[2][nBetaTubes];				//< Chris Pulser peak at reference time
};

/// Gain stabilizer wrapper with final manual tweaks
class TweakedGainStabilizer: public GainStabilizer {
public:
	/// constructor
	TweakedGainStabilizer(GainStabilizer* BG);
	/// compute gain stabilization factor for given side, PMT
	virtual float gmsFactor(Side s, unsigned int t, float time) const;
	/// print gain stabilization info
	virtual void printSummary();
	/// get a summary of GMS calibration parameters
	virtual Stringmap gmsSummary() const;
protected:
	GainStabilizer* baseGain;		//< base gain stabilization before tweaks
	float eOrig[2][nBetaTubes+1];	//< starting energy for each PMT
	float eFinal[2][nBetaTubes+1];	//< where starting energy gets scaled to
};



/*
 /// use GMS to correct LED brightness
 class LEDCorrector {
 public:
 /// constructor
 LEDCorrector(RunNum myRun, CalDB* cdb);
 /// destructor
 virtual ~LEDCorrector();
 /// get brightness of LED to reference tube (rel. to Co60 source)
 float ledBrightness(Side s, float time);
 /// combine 2 Co60 peaks into one "combined" ADC
 float adc_Co(Side s, float time);
 
 CalDB* CDB;					//< reference to calibration DB
 RunNum rn;					//< run num of this run
 
 protected:
 TGraph* refLED[2];			//< reference PMT LED ADC
 TGraph* co60Peaks[2][2];	//< reference PMT Co60 peaks
 TGraph* refLinearity[2];	//< reference PMT linearity
 TGraph* refInverses[2];		//< reference PMT inverse linearity
 };
 */

/// get LED ADC
//float getLEDadc(Side s, unsigned int t, float time);
/// get muon peak ADC
//float getMuonPeak(Side s, unsigned int t, float time) { if(muonPeaks[s][t]) return muonPeaks[s][t]->Eval(time); return 1000.0; }

//float muon0[2][4];				//< muon peak ADC at reference run t=0
//float kgms[2][4];				//< GMS factor for Kurie-endpoint based GMS
//TGraph* ledPeaks[2][4];				//< GMS LED Peaks history
//TGraph* muonPeaks[2][4];			//< Muon peaks history
//float energy_peg[2][4];				//< light equivalent of stabilized LED peak to calibration source light

#endif
