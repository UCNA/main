#ifndef CATHSEGCALIBRATOR_HH
#define CATHSEGCALIBRATOR_HH 1

#include <string>
#include <TGraph.h>
#include <vector>
#include "Types.hh"

/// class for cathode calibrations and positioning adjustement
class CathSegCalibrator: private NoCopy {
public:
	/// constructor
	CathSegCalibrator(): norm(1.0) {}
	/// destructor
	~CathSegCalibrator();
	/// adjust position for event with given energy; assumes x0 in local coordinates [-0.5,0.5]
	double adjustPos(double x0, double E) const;
	
	double pos;						//< position
	std::string channel;			//< readout channel name
	double norm;					//< signal normalization for positioning
	std::vector<TGraph*> pcoeffs;	//< positioning correction coefficients
};

#endif
