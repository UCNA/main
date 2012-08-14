#ifndef CATHODETWEAKANALYZER_HH
#define CATHODETWEAKANALYZER_HH 1

#include "RunAccumulator.hh"
#include "WirechamberCalibrator.hh"
#include "WirechamberGainAnalyzer.hh"

/// analyzer plugin for wirechamber position/calibration analysis
class CathodeTweakAnalyzer: public AnalyzerPlugin {
public:
	/// constructor
	CathodeTweakAnalyzer(RunAccumulator* RA);
	/// process a data point into position histograms
	virtual void fillCoreHists(ProcessedDataScanner& PDS, double weight);
	/// make output plots
	virtual void makePlots();
	
	fgbgPair* hitPos[2];						//< corrected hit positions on each side
	fgbgPair* hitPosRaw[2];						//< uncorrected hit positions on each side
	fgbgPair* cathHitpos[2][2][kMaxCathodes];	//< raw position distribution around each cathode by energy
};

/// analyzer for wirechamber gain and position tweaking
class MWPCTuningAnalyzer: public RunAccumulator {
public:
	/// constructor
	MWPCTuningAnalyzer(OutputManager* pnt, const std::string& nm, const std::string& inflName = "");
	/// create a new instance of this analyzer
	virtual SegmentSaver* makeAnalyzer(const std::string& nm, const std::string& inflname) { return new MWPCTuningAnalyzer(this,nm,inflname); }
	
	CathodeGainAnalyzer* myCG;		//< cathode gain
	CathodeTweakAnalyzer* myCT;		//< cathode position tweaking
};

/// process wirechamber calibration data
void processWirechamberCal(MWPCTuningAnalyzer& WCdat, MWPCTuningAnalyzer& WCsim);
/// process wirechamber calibration data from Xenon posmap specified by run range, nrings
void processWirechamberCal(RunNum r0, RunNum r1, unsigned int nrings);

#endif
