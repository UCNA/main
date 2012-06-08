#ifndef WIRECHAMBERCALIBRATOR_HH
#define WIRECHAMBERCALIBRATOR_HH 1

#include "Types.hh"
#include "CalDB.hh"
#include "CalDBSQL.hh"
#include "PositionResponse.hh"
#include "CathSegCalibrator.hh"
#include <vector>

/// enumeration for various types of wirechamber reconstruction errors
enum reconstructionErrors {
	WIRES_GOOD		= 0,	//< perfectly good wirechamber events
	WIRES_NONE		= 1<<0,	//< no wires fired
	WIRES_EDGE		= 1<<1,	//< reconstructed event on edge of wire chamber
	WIRES_SINGLET	= 1<<2,	//< single wire fired
	WIRES_DOUBLET	= 1<<3,	//< 2 wires fired
	WIRES_CLIPPED	= 1<<4,	//< clipping in wirechamber signal
	WIRES_MULTICLIP	= 1<<5,	//< multiple wire clipping
	WIRES_NONUNIF	= 1<<6	//< non-uniform wire spacing required in calculations
};

/// struct for reconstruction information about a wirechamber hit
struct wireHit {
	float center;				//< reconstructed center of charge cloud
	float width;				//< reconstructed 1-sigma width of charge cloud
	float maxValue;				//< maximum wire signal for event
	float cathodeSum;			//< sum of cathode signals for event
	unsigned int maxWire;		//< wire with largest signal
	unsigned int nClipped;		//< number of wires clipped at maximum value
	unsigned int multiplicity;	//< number of wires firing above threshold for this event
	int errflags;				//< reconstruction error warnings
	float rawCenter;			//< center position before energy-dependent tweaking
};

/// class for calibrating wirechamber data
class WirechamberCalibrator: private NoCopy {
public:
	/// constructor
	WirechamberCalibrator(RunNum rn, CalDB* cdb = CalDBSQL::getCDB());
	/// destructor
	virtual ~WirechamberCalibrator();
	
	/// calibrate anode padc signal to deposited energy
	float calibrateAnode(float adc, Side s, float x, float y, float t) const;
	
	/// get anode gain correction factor
	float wirechamberGainCorr(Side s, float t) const;
	
	/// calculate hit position from wire values array
	wireHit calcHitPos(Side s, AxisDirection d, std::vector<float>& wireValues, std::vector<float>& wirePeds);

	/// fine-tweak hit position once energy is known
	void tweakPosition(Side s, AxisDirection d, wireHit& h, double E);
	
	/// get list of cathode sensor names
	std::vector<std::string> getCathChans(Side s, AxisDirection d) const;
	
	/// display summary info
	virtual void printSummary();
	
protected:
	PositioningCorrector* anodeP;					//< anode calibration maps
	float anodeGainCorr[2];							//< anode correction factor for each side
	std::vector<CathSegCalibrator*> cathsegs[2][2];	//< cathode segments for each [side][plane]
	std::vector<double> wirePos[2][2];				//< cathode wire positions on each plane
};


#endif
