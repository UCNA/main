#ifndef WIRECHAMBERSTUDY_HH
#define WIRECHAMBERSTUDY_HH 1

#include "SectorCutter.hh"
#include "RunAccumulator.hh"
#include "WirechamberCalibrator.hh"

/// Struct for anode calibration data
struct AnodeSeg {
	Side s;					//< anode side
	float energy;			//< energy bin center
	unsigned int energybin;	//< energy bin number
	unsigned int segment;	//< position segment number
	int fiterr;				//< fitter error code
	float_err mpv;			//< most probable value
	float_err sigma;		//< width
	float_err constant;		//< height
	float npts;				//< number of data points
};

/// Struct for cathode calibration data
struct CathodeSeg {
	Side s;					//< side
	AxisDirection d;		//< plane measuring direction
	unsigned int i;			//< cathode number
	float_err height;		//< normalized height
	float_err width;		//< width
	float pos;				//< cathode position
};

/// class for position-segmented analysis
class WirechamberAnalyzer: public RunAccumulator {
public:
	/// constructor
	WirechamberAnalyzer(OutputManager* pnt, const std::string& nm, float r, unsigned int nr, const std::string& infl="");
	
	/// cloning generator: just return another of the same subclass (with any settings you want to preserve)
	virtual SegmentSaver* makeAnalyzer(const std::string& nm,
									   const std::string& inflname) { return new WirechamberAnalyzer(this,nm,sects.r,sects.n,inflname); }
	
	/// process a data point into position histograms
	virtual void fillCoreHists(ProcessedDataScanner& PDS, double weight);
	
	/// fit enpoints in each sector
	virtual void calculateResults();
	/// make output plots
	virtual void makePlots();
	/// virtual routine for MC/Data comparison plots/calculations
	virtual void compareMCtoData(RunAccumulator& OAdata);
	
	double eMax;											//< maximum energy bin
	unsigned int nEnergyBins;								//< number of energy bins
	fgbgPair masterEnergySpectrum;							//< master energy spectrum (for binning points)
	std::vector< std::vector<fgbgPair> > anodeSpectra[2];	//< anode spectra in each energy bin for each position
	std::vector< std::vector<AnodeSeg> > anodeDat[2];		//< processed anode data
	std::vector<float> cathPos[2][2];						//< cathode positions array
	std::vector<fgbgPair> cathHists[2][2];					//< normalized cathode histograms for each side,plane,cathode
	std::vector<TH1D*> cathSlices[2][2][2];					//< cathode histogram slice [side][direction][center/width]
	std::vector<CathodeSeg> cathDat[2][2];					//< processed cathode data
	fgbgPair cathHitDist[2][2][kMaxCathodes];				//< position distribution around each cathode
	SectorCutter sects;										//< sector cutter
};

/// process runs with wirechamber analyzer
void runWirechamberAnalyzer(RunNum r0, RunNum r1, unsigned int nrings);


#endif
