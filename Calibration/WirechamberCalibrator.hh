#ifndef WIRECHAMBERCALIBRATOR_HH
#define WIRECHAMBERCALIBRATOR_HH

#include "Types.hh"
#include "CalDB.hh"
#include "CalDBSQL.hh"
#include "PositionResponse.hh"
#include "CathSegCalibrator.hh"
#include "ManualInfo.hh"
#include <vector>
#include <map>

/// maximum number of MWPC cathode segments per plane (actual number may be less if some dead)
#define kMaxCathodes 16
/// number of individual cathode wires per segment
#define kWiresPerCathode 4

/// enumeration for various types of wirechamber reconstruction errors
enum reconstructionErrors {
	WIRES_GOOD		= 0,	///< perfectly good wirechamber events
	WIRES_NONE		= 1<<0,	///< no wires fired
	WIRES_EDGE		= 1<<1,	///< reconstructed event on edge of wire chamber
	WIRES_SINGLET	= 1<<2,	///< single wire fired
	WIRES_DOUBLET	= 1<<3,	///< 2 wires fired
	WIRES_CLIPPED	= 1<<4,	///< clipping in wirechamber signal
	WIRES_MULTICLIP	= 1<<5,	///< multiple wire clipping
	WIRES_NONUNIF	= 1<<6	///< non-uniform wire spacing required in calculations
};

/// struct for reconstruction information about a wirechamber hit
class wireHit {
public:
	float center;				///< reconstructed center of charge cloud
	float width;				///< reconstructed 1-sigma width of charge cloud
	float maxValue;				///< maximum wire signal for event
	float cathodeSum;			///< sum of cathode signals for event
	unsigned int maxWire;		///< wire with largest signal
	unsigned int nClipped;		///< number of wires clipped at maximum value
	unsigned int multiplicity;	///< number of wires firing above threshold for this event
	int errflags;				///< reconstruction error warnings
	float rawCenter;			///< center position before energy-dependent tweaking
	float height;				///< reconstructed amplitude of charge cloud
	
	/// charge cloud size
	float ccloud_size() const { return height*width; }
};

/// class for calibrating wirechamber data
class WirechamberCalibrator: private NoCopy {
public:
	/// constructor
	WirechamberCalibrator(RunNum rn, CalDB* cdb = CalDBSQL::getCDB());
	/// destructor
	virtual ~WirechamberCalibrator();
	
	/// calibrate wire energy
	float wirechamberEnergy(Side s, const wireHit& x_wires, const wireHit& y_wires, const MWPCevent& mwpce) const;
	
	/// get gain correction in use
	float wirechamberGainCorr(Side s) const { return mwpcGainCorr[s]; }
	/// get charge by method
	float chargeProxy(Side s, ChargeProxyType c, const wireHit& x_wires, const wireHit& y_wires, const MWPCevent& mwpce) const;
	/// get alternate energy calibration method
	const MWPC_Ecal_Spec& getAltEcal(Side s, ChargeProxyType tp) const;
	
	/// calculate hit position from wire values array; optionally, calculate other centering variables
	wireHit calcHitPos(Side s, AxisDirection d, const float* cathADC, const float* cathPed = NULL) const;
	
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
	/// get charge cloud gain factor for given cathode
	float getCathCCloudGain(Side s, AxisDirection d, unsigned int c) const;
	/// get charge cloud gain factor for overall wirechamber
	float getCcloudGain(Side s) const { return ccloudGainCorr[s]; }
	
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
	
	double sigma;	///< expected charge cloud width for default reconstruction, in same units as cathode positions
	
	/// draw wire positions on plot
	void drawWires(Side s, AxisDirection p, TVirtualPad* C, Int_t color = 4, AxisDirection onAxis=X_DIRECTION) const;
	
	/// turn on/off cathode shape calibrations (for baseline comparison)
	static bool calibrateCathodes;
	
	// wirechamber trigger cut information
	CutVariable fCathSum[BOTH];		///< combined x+y cathode sum on each side
	CutVariable fCathMax[BOTH];		///< min(max cathode each plane) for each side
	CutVariable fCathMaxSum[BOTH];	///< sum of max cathode from each plane for each side

	
	// some specialized data for detector response simulation
	unsigned int nWires(Side s, AxisDirection d) const { if(s>WEST||d>Y_DIRECTION) return 0; return nCaths[s][d]; }
	std::vector<std::string> cathNames[BOTH][2];		///< cathode sensor names
	std::vector<float> cathPeds0[BOTH][2];				///< t=0 cathode pedestal values [side][plane]
	std::vector<float> cathPedW0[BOTH][2];				///< t=0 cathode pedesatl widths [side][plane]
	std::vector<double> cath_ccloud_gains[BOTH][2];		///< gain conversion factors between cathode signal and charge cloud size
	std::vector<float> cathseg_energy_norm[BOTH][2];	///< cathode ADC to energy normalization ADC/(Ew*f*eta) = g_i/cnorm
	PositioningCorrector* ccloud_eta[BOTH];				///< position map for charge cloud method
	
protected:
	unsigned int nCaths[BOTH][2];									///< number of cathodes on each [side][plane]
	std::map<ChargeProxyType, MWPC_Ecal_Spec> ecalMethods[BOTH];	///< alternate energy calibration methods
	ChargeProxyType myChargeProxy[BOTH];							///< what signal to use for energy calibration
	PositioningCorrector* chargeP[BOTH];							///< charge signal position calibration maps
	float mwpcGainCorr[BOTH];										///< gain correction factor for each side
	float ccloudGainCorr[BOTH];										///< charge cloud method gain correction factor for each side
	std::vector<CathSegCalibrator*> cathsegs[BOTH][2];				///< cathode segments for each [side][plane]
	std::vector<double> wirePos[BOTH][2];							///< cathode wire positions on each plane
	std::vector<double> domains[BOTH][2];							///< dividing lines between ``domains'' of each wire
};


#endif
