#ifndef CALDB_HH
#define CALDB_HH 1

#include <TGraphErrors.h>
#include <vector>
#include "Types.hh"
#include "PositionResponse.hh"
#include "EfficCurve.hh"
#include "RunInfo.hh"
#include "CathSegCalibrator.hh"

/// abstract class for interface to database of calibration settings and peaks
class CalDB {
public:
	/// destructor
	virtual ~CalDB() {}
	
	/// check if valid data available for run
	virtual bool isValid(RunNum rn) = 0;
	
	/// get linearity correction data
	virtual TGraph* getLinearity(RunNum rn, Side s, unsigned int t) = 0;
	
	/// get noise estimate calibration point ADC width
	virtual float getNoiseWidth(RunNum rn, Side s, unsigned int t) = 0;
	/// get noise estimate calibration point raw ADC value
	virtual float getNoiseADC(RunNum rn, Side s, unsigned int t) = 0;
	
	/// get a run monitor graph
	virtual TGraphErrors* getRunMonitor(RunNum rn, const std::string& sensorName, const std::string& monType, bool centers = true) = 0;
	/// get initial value for run monitor graph
	virtual float getRunMonitorStart(RunNum rn, const std::string& sensorName, const std::string& monType) = 0;
	
	/// get pedestals for named sensor
	virtual TGraph* getPedestals(RunNum rn, const std::string& sensorName) = 0;
	/// get pedestal widths for named sensor
	virtual TGraph* getPedwidths(RunNum rn, const std::string& sensorName) = 0;
	
	/// get Calibrations DB name
	virtual std::string getName() const = 0;
	
	/*
	/// get LED peak data for named sensor
	virtual TGraph* getLED(RunNum rn, const std::string& sensorName) = 0;
	/// get Co60 reference tube data
	virtual TGraph* getCo60(RunNum rn, Side s, unsigned int peakNum) = 0;
	/// get muon peak history for given sensor
	virtual TGraph* getMuonPeak(RunNum rn, const std::string& sensorName) = 0;
	/// get kurie endpoint calibration energy for given run, side, PMT
	virtual float getKurieEnergy(RunNum rn, Side s, unsigned int t) = 0;
	/// get kurie endpoint calibration energy for given run, side, PMT
	virtual float getKurieADC(RunNum rn, Side s, unsigned int t) = 0;
	*/
	
	/// get GMS Cal run number
	virtual RunNum getGMSRun(RunNum rn) = 0;
	
	/// get energy calibration point ADC value
	virtual float getEcalADC(RunNum rn, Side s, unsigned int t) = 0;
	/// get energy calibration point visible energy
	virtual float getEcalEvis(RunNum rn, Side s, unsigned int t) = 0;
	/// get energy calibration point x position
	virtual float getEcalX(RunNum rn, Side s) = 0;
	/// get energy calibration point y position
	virtual float getEcalY(RunNum rn, Side s) = 0;
	
	/// get positioning corrector for given run
	virtual PositioningCorrector* getPositioningCorrector(RunNum rn) = 0;
	/// get anode positioning corrector for given run
	virtual PositioningCorrector* getAnodePositioningCorrector(RunNum rn) = 0;
	/// get anode gain correction factor for run
	virtual float getAnodeGain(RunNum rn, Side s) = 0;
	/// get GMS gain tweaking factors
	virtual void getGainTweak(RunNum rn, Side s, unsigned int t, float& orig, float& final) = 0;
	
	/// get trigger efficiency function
	virtual EfficCurve* getTrigeff(RunNum rn, Side s, unsigned int t) = 0;
	/// get E_vis -> E_true parametrization
	virtual TGraph* getEvisConversion(RunNum rn, Side s, EventType tp) = 0;
	/// get list of cathode segment calibrators (caller is responsible for deletion)
	virtual std::vector<CathSegCalibrator*> getCathSegCalibrators(RunNum rn, Side s, AxisDirection d) = 0;
	
	/// run start time
	virtual int startTime(RunNum rn, int t0 = 0) { return 0; }
	/// run end time
	virtual int endTime(RunNum rn, int t0 = 0) { return 0; }
	/// run time after cuts, blinded
	virtual BlindTime fiducialTime(RunNum rn) = 0;
	/// total run time
	virtual float totalTime(RunNum rn) { return 0; }
	/// get RunInfo for given run
	virtual RunInfo getRunInfo(RunNum r) = 0;
};

#endif
