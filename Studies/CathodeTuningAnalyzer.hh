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
	float_err height;		//< normalized height
	float_err width;		//< width
	float_err center;		//< measured center
	float_err max;			//< maximum normalized value
	float fill_frac;		//< proportion of expected events
	float pos;				//< cathode position
};
/// convert cathode segment to Stringmap
Stringmap cathseg2sm(const CathodeSeg& c);
/// convert Stringmap to CathodeSeg
CathodeSeg sm2cathseg(const Stringmap& m);

/// analyzer plugin for wirechamber position/calibration analysis
class CathodeGainPlugin: public AnalyzerPlugin {
public:
	/// constructor
	CathodeGainPlugin(RunAccumulator* RA);
	/// fill histograms
	virtual void fillCoreHists(ProcessedDataScanner& PDS, double weight);
	
	fgbgPair* cathNorm[BOTH][2][kMaxCathodes];		//< cathode normalization histograms by [side][plane][cathode]
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

/// fit/plot CathodeGainPlugin histograms
void processCathNorm(CathodeGainPlugin& CGA);

/// process wirechamber calibration data
void processWirechamberCal(CathodeTuningAnalyzer& WCdat, CathodeTuningAnalyzer& WCsim);
/// process wirechamber calibration data from Xenon posmap specified by run range, nrings
void processWirechamberCal(RunNum r0, RunNum r1, unsigned int nrings);

#endif
