#ifndef CATHODETWEAKANALYZER_HH
#define CATHODETWEAKANALYZER_HH 1

#include "RunAccumulator.hh"
#include "WirechamberCalibrator.hh"
#include "WirechamberAnodePlugins.hh"

/// Struct for cathode calibration data
struct CathodeSeg {
	Side s;					//< side
	AxisDirection d;		//< plane measuring direction
	unsigned int i;			//< cathode number
	float n_exp;			//< expected number from data
	float n_obs;			//< observed number from simulation
	float dndx_lo;			//< dn/dx at low edge, observed, normalized corrected position
	float dndx_hi;			//< dn/dx at high edge, observed, normalized corrected position
	float pos;				//< cathode position
};
/// convert cathode segment to Stringmap
Stringmap cathseg2sm(const CathodeSeg& c);

/// analyzer plugin for wirechamber position/calibration analysis
class CathodeGainPlugin: public AnalyzerPlugin {
public:
	/// constructor
	CathodeGainPlugin(RunAccumulator* RA);
	/// fill histograms
	virtual void fillCoreHists(ProcessedDataScanner& PDS, double weight);
	/// fit cathode shapes
	virtual void calculateResults();
	/// make output plots
	virtual void makePlots();
	
	fgbgPair* cathNorm[BOTH][2][kMaxCathodes];				//< cathode normalization histograms by [side][plane][cathode]
	std::vector<TH1D*> slicefits[BOTH][2][kMaxCathodes];	//< Gaussian fit parameters for cathode response shape at each position
};

/// analyzer plugin for wirechamber position/calibration analysis
class CathodeTweakPlugin: public AnalyzerPlugin {
public:
	/// constructor
	CathodeTweakPlugin(RunAccumulator* RA);
	/// process a data point into position histograms
	virtual void fillCoreHists(ProcessedDataScanner& PDS, double weight);
	/// make output plots
	virtual void makePlots();
	
	fgbgPair* hitPos[2];						//< corrected hit positions on each side
	fgbgPair* hitPosRaw[2];						//< uncorrected hit positions on each side
	fgbgPair* cathHitpos[2][2][kMaxCathodes];	//< raw position distribution around each cathode by energy
};

/// analyzer for wirechamber gain and position tweaking
class CathodeTuningAnalyzer: public RunAccumulator {
public:
	/// constructor
	CathodeTuningAnalyzer(OutputManager* pnt, const std::string& nm, const std::string& inflName = "");
	/// create a new instance of this analyzer
	virtual SegmentSaver* makeAnalyzer(const std::string& nm, const std::string& inflname) { return new CathodeTuningAnalyzer(this,nm,inflname); }
	
	CathodeGainPlugin* myCG;		//< cathode gain
	CathodeTweakPlugin* myCT;		//< cathode position tweaking
};

/// comparison analyzer for simulation
class SimCathodeTuningAnalyzer: public RunAccumulator {
public:
	/// constructor
	SimCathodeTuningAnalyzer(OutputManager* pnt, const std::string& nm, const std::string& inflName = "");
	/// create a new instance of this analyzer
	virtual SegmentSaver* makeAnalyzer(const std::string& nm, const std::string& inflname) { return new SimCathodeTuningAnalyzer(this,nm,inflname); }
	
	CathodeTweakPlugin* myCT;		//< cathode position tweaking
};

/// fit/plot CathodeTweakPlugin cathode shape correction data
void processCathTweak(CathodeTweakPlugin& CTDat, CathodeTweakPlugin& CTSim);

#endif
