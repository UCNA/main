#ifndef WIRECHAMBERCALIBRATOR_HH
#define WIRECHAMBERCALIBRATOR_HH 1

#include "Types.hh"
#include "CalDB.hh"
#include "CalDBSQL.hh"
#include "PositionResponse.hh"
#include "CathSegCalibrator.hh"
#include <vector>

/// maximum number of MWPC cathode segments per plane (actual number may be less if some dead)
#define kMaxCathodes 16
/// number of individual cathode wires per segment
#define kWiresPerCathode 4

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
	float height;				//< reconstructed amplitude of charge cloud
};

/// struct for approximate wireplane hit reconstruction
struct hitRecon {
	float c;	//< center
	float w;	//< width
	float h;	//< height
	float a;	//< area
	int flags;	//< reconstruction flags
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
	
	/// calculate hit position from wire values array; optionally, calculate other centering variables
	wireHit calcHitPos(Side s, AxisDirection d, std::vector<float>& wireValues, std::vector<float>& wirePeds, hitRecon* hLorentzian = NULL) const;
	
	/// special case for calculating hit position from two points
	void calcDoubletHitPos(wireHit& h, float x0, float x1, float y0, float y1) const;
	
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
	
	double sigma;	//< expected charge cloud width for default reconstruction, in same units as cathode positions
	
	/// calculate ``Lorentzian center'' from three points
	static hitRecon calcLorentzCenter(double sm, double s0, double sp, double pm, double p0, double pp);
	/// draw wire positions on plot
	void drawWires(Side s, AxisDirection p, TVirtualPad* C, Int_t color = 4, AxisDirection onAxis=X_DIRECTION) const;
	
	/// turn on/off cathode shape calibrations (for baseline comparison)
	static bool calibrateCathodes;
	
protected:
	PositioningCorrector* anodeP;						//< anode calibration maps
	float anodeGainCorr[BOTH];							//< anode correction factor for each side
	std::vector<CathSegCalibrator*> cathsegs[BOTH][2];	//< cathode segments for each [side][plane]
	std::vector<double> wirePos[BOTH][2];				//< cathode wire positions on each plane
	std::vector<double> domains[BOTH][2];				//< dividing lines between ``domains'' of each wire
};


#endif
