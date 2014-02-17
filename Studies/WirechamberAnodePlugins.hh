#ifndef WIRECHAMBERGAINANALYZER
#define WIRECHAMBERGAINANALYZER 1

#include "RunAccumulator.hh"
#include "WirechamberCalibrator.hh"
#include <TGraphErrors.h>

/// analyzer plugin for monitoring overall wirechamber gain
class AnodeGainPlugin: public AnalyzerPlugin {
public:
	/// constructor
	AnodeGainPlugin(RunAccumulator* RA);
	/// destructor
	~AnodeGainPlugin();
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
class WirechamberSimBackscattersPlugin: public AnalyzerPlugin {
public:
	/// constructor
	WirechamberSimBackscattersPlugin(RunAccumulator* RA);
	/// fill histograms
	virtual void fillCoreHists(ProcessedDataScanner& PDS, double weight);
	/// generate Type II/III separation data
	void make23SepInfo(OutputManager& OM);
	
	fgbgPair* anodeTypeID[BOTH][TYPE_III_EVENT+1];		//< anode energy spectrum binned by event energy by [side][type]
	fgbgPair* anodeNormCoords[BOTH][TYPE_III_EVENT+1];	//< anode energy spectrum in cut-normalized coordinates by [side][type]
};

/// extract properties of wirechamber charge proxies
class WirechamberChargeProxyPlugin: public AnalyzerPlugin {
public:
	/// constructor
	WirechamberChargeProxyPlugin(RunAccumulator* RA);
	/// fill histograms
	virtual void fillCoreHists(ProcessedDataScanner& PDS, double weight);
	
	fgbgPair* cathsumVanode[BOTH];		//< cathode sum versus anode
	fgbgPair* chargesizeVanode[BOTH];	//< fit charge cloud size vs. anode
};

#endif
