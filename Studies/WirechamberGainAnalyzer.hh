#ifndef WIRECHAMBERGAINANALYZER
#define WIRECHAMBERGAINANALYZER 1

#include "RunAccumulator.hh"
#include "WirechamberCalibrator.hh"
#include <TGraphErrors.h>

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
class CathodeGainAnalyzer: public AnalyzerPlugin {
public:
	/// constructor
	CathodeGainAnalyzer(RunAccumulator* RA);
	/// fill histograms
	virtual void fillCoreHists(ProcessedDataScanner& PDS, double weight);
	
	fgbgPair* cathNorm[BOTH][2][kMaxCathodes];		//< cathode normalization histograms by [side][plane][cathode]
};

/// analyzer plugin for monitoring overall wirechamber gain
class AnodeGainAnalyzer: public AnalyzerPlugin {
public:
	/// constructor
	AnodeGainAnalyzer(RunAccumulator* RA);
	/// destructor
	~AnodeGainAnalyzer();
	/// process a data point into position histograms
	virtual void fillCoreHists(ProcessedDataScanner& PDS, double weight);
	/// calculations on filled histograms
	virtual void calculateResults();
	/// make output plots
	virtual void makePlots();
	/// MC/data comparison
	virtual void compareMCtoData(AnalyzerPlugin* AP);
	
	fgbgPair* anodeCal[BOTH][TYPE_III_EVENT+1];	//< anode energy spectrum binned by event energy by [side][type]
	fgbgPair* anodeGaincorr;					//< average anode gain correction factor in use
	fgbgPair* norm23[BOTH];						//< "normalized" Type II, III anode signals by [side]
	
	std::vector<TH1F*> hSlices[BOTH][TYPE_III_EVENT+1];	//< fit anode spectrum "slices"
	TGraphErrors* gAnode[BOTH][TYPE_III_EVENT+1];		//< anode energy distribution fit
};

/// analyzer plugin for evaluating simulated Type II/III MWPC energy deposition split
class WirechamberSimTypeID: public AnalyzerPlugin {
public:
	/// constructor
	WirechamberSimTypeID(RunAccumulator* RA);
	/// fill histograms
	virtual void fillCoreHists(ProcessedDataScanner& PDS, double weight);
	/// generate Type II/III separation data
	void make23SepInfo(OutputManager& OM);
	
	fgbgPair* anodeTypeID[BOTH][TYPE_III_EVENT+1];		//< anode energy spectrum binned by event energy by [side][type]
	fgbgPair* anodeNormCoords[BOTH][TYPE_III_EVENT+1];	//< anode energy spectrum in cut-normalized coordinates by [side][type]
};

/// extract properties of wirechamber charge proxies
class WirechamberChargeProxyAnalyzer: public AnalyzerPlugin {
public:
	/// constructor
	WirechamberChargeProxyAnalyzer(RunAccumulator* RA);
	/// fill histograms
	virtual void fillCoreHists(ProcessedDataScanner& PDS, double weight);
	
	fgbgPair* cathsumVanode[BOTH];		//< cathode sum versus anode
	fgbgPair* chargesizeVanode[BOTH];	//< fit charge cloud size vs. anode
};

#endif
