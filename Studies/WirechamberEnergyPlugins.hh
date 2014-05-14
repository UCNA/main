#ifndef WIRECHAMBERGAINANALYZER
#define WIRECHAMBERGAINANALYZER 1

#include "RunAccumulator.hh"
#include "WirechamberCalibrator.hh"
#include <TGraphErrors.h>

/// analyzer plugin for monitoring overall wirechamber gain
class MWPCGainPlugin: public AnalyzerPlugin {
public:
	/// constructor
	MWPCGainPlugin(RunAccumulator* RA);
	/// destructor
	~MWPCGainPlugin();
	/// process a data point into position histograms
	virtual void fillCoreHists(ProcessedDataScanner& PDS, double weight);
	/// calculations on filled histograms
	virtual void calculateResults();
	/// make output plots
	virtual void makePlots();
	/// MC/data comparison
	virtual void compareMCtoData(AnalyzerPlugin* AP);
	
	fgbgPair* mwpcE[BOTH][TYPE_III_EVENT+1];	///< MWPC energy spectrum binned by event energy by [side][type]
	fgbgPair* mwpcGaincorr;						///< average MWPC gain correction factor already use
	fgbgPair* norm23[BOTH];						///< "normalized" Type II, III energy signals by [side]
	
	std::vector<TH1F*> hSlices[BOTH][TYPE_III_EVENT+1];	///< fit energy spectrum "slices"
	std::vector< std::vector<double> > fitParams[BOTH][TYPE_III_EVENT+1];	///< Landau fit parameters for each scintillator energy bin
	std::vector< std::vector<double> > fitErrs[BOTH][TYPE_III_EVENT+1];		///< Landau fit errors for each scintillator energy bin
	
	TGraphErrors* gEw[BOTH][TYPE_III_EVENT+1];			///< wirechamber energy deposition vs scintillator energy
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
	
	fgbgPair* EwTypeID[BOTH][TYPE_III_EVENT+1];		///< anode energy spectrum binned by event energy by [side][type]
	fgbgPair* EwNormCoords[BOTH][TYPE_III_EVENT+1];	///< anode energy spectrum in cut-normalized coordinates by [side][type]
};

/// analyzer plugin for showing CathMaxSum calibration
class WirechamberCathMaxSumPlugin: public AnalyzerPlugin {
public:
	/// constructor
	WirechamberCathMaxSumPlugin(RunAccumulator* RA);
	/// fill histograms
	virtual void fillCoreHists(ProcessedDataScanner& PDS, double weight);
	/// MC/data comparison
	virtual void compareMCtoData(AnalyzerPlugin* AP);

	fgbgPair* mwpcCathMax[BOTH];		///< cathode max sum spectra on each side
};

/// analyzer plugin for evaluting wirechamber trigger efficiency curves
class WirechamberSimTrigEfficPlugin: public AnalyzerPlugin {
public:
	/// constructor
	WirechamberSimTrigEfficPlugin(RunAccumulator* RA);
	/// fill histograms
	virtual void fillCoreHists(ProcessedDataScanner& PDS, double weight);
	/// make output plots
	virtual void makePlots();
	
	TH1* mwpcHitEffic[BOTH][2];			///< MWPC event spectra for [side][triggered]
	TH1* mwpcCathMax[BOTH][2];			///< max cathode signal for events with and without actual energy deposition
};


#endif
