#ifndef WIRECHAMBERGAINANALYZER
#define WIRECHAMBERGAINANALYZER 1

#include "RunAccumulator.hh"
#include "WirechamberCalibrator.hh"

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
class WirechamberGainAnalyzer: public AnalyzerPlugin {
public:
	/// constructor
	WirechamberGainAnalyzer(RunAccumulator* RA);
	/// process a data point into position histograms
	virtual void fillCoreHists(ProcessedDataScanner& PDS, double weight);
	/// calculations on filled histograms
	virtual void calculateResults();
	/// make output plots
	virtual void makePlots();
	/// MC/data comparison
	virtual void compareMCtoData(AnalyzerPlugin* AP);
	
	fgbgPair* cathNorm[2][2][kMaxCathodes];	//< cathode normalization histograms
	fgbgPair* anodeCal[2];					//< anode calibration spectrum (Type 0, Erecon>225)
};



#endif
