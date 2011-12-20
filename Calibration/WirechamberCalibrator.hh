#ifndef WIRECHAMBERCALIBRATOR_HH
#define WIRECHAMBERCALIBRATOR_HH 1

#include "CalDB.hh"
#include "PositionResponse.hh"

/// class for calibrating wirechamber data
class WirechamberCalibrator {
public:
	/// constructor
	WirechamberCalibrator(RunNum rn, CalDB* cdb);
	/// destructor
	virtual ~WirechamberCalibrator() {}
	
	/// calibrate anode padc signal to deposited energy
	float calibrateAnode(float adc, Side s, float x, float y, float t) const;
	
	/// get anode gain correction factor
	float wirechamberGainCorr(Side s, float t) const;
	
	/// print calibrations summary
	void printWirecalSummary() const;
	
private:
	PositioningCorrector* anodeP;	//< anode calibration maps
	float anodeGainCorr[2];			//< anode correction factor for each side
};


#endif
