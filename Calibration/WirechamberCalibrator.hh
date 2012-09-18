#ifndef WIRECHAMBERCALIBRATOR_HH
#define WIRECHAMBERCALIBRATOR_HH 1

#include "Types.hh"
#include "CalDB.hh"
#include "CalDBSQL.hh"
#include "PositionResponse.hh"
#include "CathSegCalibrator.hh"
#include <vector>

/// maximum number of MWPC cathode wires per plane (actual number may be less if some dead)
#define kMaxCathodes 16

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
	wireHit calcHitPos(Side s, AxisDirection d, std::vector<float>& wireValues, std::vector<float>& wirePeds) const;

	/// fine-tweak hit position once energy is known
	void tweakPosition(Side s, AxisDirection d, wireHit& h, double E) const;
	
	/// get number of cathodes
	unsigned int getNCaths(Side s, AxisDirection d) const { return cathsegs[s][d].size(); }
	
	/// get list of cathode sensor names
	std::vector<std::string> getCathChans(Side s, AxisDirection d) const;
	
	/// get normalization for given cathode
	float getCathNorm(Side s, AxisDirection d, unsigned int c) const;
	
	/// display summary info
	virtual void printSummary();
	
	/// summary data about wirechamber calibration
	Stringmap wirecalSummary() const;
	
	/// convert to local normalized position
	void toLocal(Side s, AxisDirection d, float x, unsigned int& n, float& c) const;
	/// convert back from local normalized position
	float fromLocal(Side s, AxisDirection d, unsigned int n, float c) const;
	
	/// Type II/III separation MWPC energy cut as a function of scintillator energy
	static float sep23Cut(Side s, float Escint);
	/// normalized MWPC signal for Type II/III events
	static float normMWPC(Side s, float Escint, float Emwpc);
	/// Type II/III separation probability
	static float sep23Prob(Side s, float Escint, float Emwpc);
	
protected:
	PositioningCorrector* anodeP;					//< anode calibration maps
	float anodeGainCorr[2];							//< anode correction factor for each side
	std::vector<CathSegCalibrator*> cathsegs[2][2];	//< cathode segments for each [side][plane]
	std::vector<double> wirePos[2][2];				//< cathode wire positions on each plane
	std::vector<double> domains[2][2];				//< dividing lines between ``domains'' of each wire
};


#endif
