#ifndef PMTCALIBRATOR_HH
#define PMTCALIBRATOR_HH

#include "EnergyCalibrator.hh"


/// PMT reconstruction class
class PMTCalibrator: public LinearityCorrector, public PedestalCorrector, public EvisConverter, public WirechamberCalibrator  {
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
	float getClipThreshold(Side s, unsigned int t) { smassert(s<=WEST && t<nBetaTubes); return clipThreshold[s][t]; }
	
protected:
	
	float clipThreshold[2][nBetaTubes];		///< threshold to de-weight ADC in tube combination due to "clipping"
	EfficCurve* pmtEffic[2][nBetaTubes];	///< efficiency curves for each PMT
	bool disablePMT[2][nBetaTubes];			///< flags to disable summing PMT into combined result
};

#endif
