#ifndef PEAK_HH
#define PEAK_HH 1

#include <vector>
#include <TF1.h>
#include "QFile.hh"
#include "EnergyCalibrator.hh"

/// spectrum peak types
enum PeakType {
	UNKNOWN_PEAK= 0,		//< Unclassified peak type
	REF_CO60_1	= 1,		//< 1st peak of GMS reference Co60 source
	REF_CO60_2	= 2,		//< 2nd peak of GMS reference Co60 source
	TUBE1_LED	= 3,		//< GMS PMT1 LED
	TUBE2_LED	= 4,		//< GMS PMT2 LED
	TUBE3_LED	= 5,		//< GMS PMT3 LED
	TUBE4_LED	= 6,		//< GMS PMT4 LED
	REF_LED		= 7,		//< GMS reference LED
	BI_PEAK_1	= 8,		//< 1st peak of Bi source ~481keV
	BI_PEAK_2	= 9,		//< 2nd peak of Bi source ~974keV + 1047keV
	BI_PEAK_3	= 10,		//< low energy auger peak of Bi source
	SN_PEAK		= 11,		//< Sn source ~364keV
	SR85_PEAK	= 12,		//< Sr source ~500.1keV
	CD109_PEAK	= 13,		//< Cd109 source ~75.2KeV
	IN114_PEAK	= 14,		//< Indium 114 peak
	CE139_PEAK	= 15,		//< Cerium 139 peak, 126.93 (17.1%) + 159.59 (2.3%)
	M000_PEAK	= 16,		//< imaginary Monosmuth peak
	M500_PEAK	= 17,		//< imaginary 500keV monoenergetic peak
	LINE_PEAK	= 18,		//< monoenergetic line peak
	BI_COINC	= 19		//< Bi207 coincidence peak between 500 and 1000
};

/// location of a peak in the spectrum
class SpectrumPeak {
public:	
	
	/// constructor for given peak type
	SpectrumPeak(PeakType t = UNKNOWN_PEAK, unsigned int sid = 0, float x = 0, float y = 0):
	center(0), width(0), xpos(x), ypos(y), type(t), sID(sid), simulated(false) {};
	/// print to stdout
	void print() const;
		
	/// convert to Stringmap
	Stringmap toStringmap() const;
	/// read from Stringmap
	static SpectrumPeak fromStringmap(const Stringmap& M);
	/// load from a Gaussian fit (or multi-gaussian peak n)
	void fromGaussian(const TF1* f, unsigned int n = 0);
	
	/// known peak energies in keV
	Float_t energy() const;
	
	/// printable name for a given peak
	const char* name() const;
	
	float_err center;		//< peak location, raw pedestal-subtracted ADC
	float_err width;		//< peak width, raw ADC counts
	float_err energyCenter;	//< energy-calibrated center position
	float_err energyWidth;	//< energy-calibrated peak width
	float nPE;				//< actual nPE at energyCenter
	float xpos;				//< wirechamber position x center
	float ypos;				//< wirechamber position y center
	float_err h;			//< height of peak
	float integral;			//< 1-sigma integral of peak
	PeakType type;			//< type of peak
	unsigned int sID;		//< ID number of corresponding source
	bool simulated;			//< whether this peak is from simulated data
	
};

/// quick gaussian peak fitter
SpectrumPeak gausFitter(std::vector<float> dat, Float_t mu0, Float_t sigma0, int nSigma = 3);

#endif
