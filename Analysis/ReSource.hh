#ifndef RESOURCE_HH
#define RESOURCE_HH

#include "Source.hh"
#include "Types.hh"
#include "EnergyCalibrator.hh"
#include "OutputManager.hh"
#include "RunAccumulator.hh"
#include "WirechamberEnergyPlugins.hh"
#include <vector>
#include <TH1F.h>
#include <TProfile.h>

/// analyzer plugin for calibration source events
class SourceHitsPlugin: public AnalyzerPlugin {
public:
	/// constructor
	SourceHitsPlugin(RunAccumulator* RA, const Source& s, PMTCalibrator* P);

	Source mySource;		///< source being fit for
	PMTCalibrator* PCal;	///< PMTCalibrator for estimating energy scale
	bool dbgplots;			///< whether to print extra "debugging" plots
	
	unsigned int nBins;		///< histogram binning
	float eMin;				///< histogram lower energy
	float eMax;				///< histogram upper energy
	float pkMin;			///< minimum peak value (avoid Bi auger peak false alarms)
	float nSigma;			///< number of sigma to fit peaks
	
	/// process a data point into position histograms
	virtual void fillCoreHists(ProcessedDataScanner& PDS, double weight);
	/// calculations on filled histograms
	virtual void calculateResults();
	/// make output plots
	virtual void makePlots();
	/// MC/data comparison
	virtual void compareMCtoData(AnalyzerPlugin* AP);
	
	// raw event counts histograms
	fgbgPair* hTubes[nBetaTubes+1][TYPE_III_EVENT+1];	///< calibrated PMT energy spectra by event type
	fgbgPair* hTubesRaw[nBetaTubes+1];					///< raw ADC energy spectra
	fgbgPair* hErec;									///< total reconstructed energy histogram
	fgbgPair* hitPos[Y_DIRECTION+1];					///< x,y hit profiles
	
	// calculated rates histograms for fitting/plotting
	TH1* hTubesR[nBetaTubes+1][TYPE_III_EVENT+1];	///< calibrated PMT energy spectra by event type, rate
	TH1* hTubesRawR[nBetaTubes+1];					///< raw ADC energy spectra, rate
	TH1* hErecR;									///< total reconstructed energy histogram, rate
	TH1* hitPosR[Y_DIRECTION+1];					///< x,y hit profiles, rate
	std::vector<SpectrumPeak> tubePeaks[nBetaTubes+1];	///< found peaks for each PMT
	
	/// get number of counts
	unsigned int counts(EventType tp = TYPE_0_EVENT) const { return (unsigned int)(hTubes[nBetaTubes][tp]->h[GV_OPEN]->Integral()); }
};

/// simple source hit positions plugin
class SourcePositionsPlugin: public AnalyzerPlugin {
public:
	/// constructor
	SourcePositionsPlugin(RunAccumulator* RA);
	
	/// process a data point into position histograms
	virtual void fillCoreHists(ProcessedDataScanner& PDS, double weight);
	/// make output plots
	virtual void makePlots();
	
	TH2F* hitPos[BOTH];		///< hit positions on each side
};

/// Wirechamber energy calibrations analyzer
class SourceHitsAnalyzer: public RunAccumulator {
public:
	/// constructor
	SourceHitsAnalyzer(OutputManager* pnt, const std::string& nm, const std::string& inflname = "");
	
	/// cloning generator: just return another of the same subclass (with any settings you want to preserve)
	virtual SegmentSaver* makeAnalyzer(const std::string& nm, const std::string& inflname);

	/// add a source to analyze
	void addSource(const Source& s);
	
	PMTCalibrator* PCal;						///< PMTCalibrator for estimating energy scale
	SourcePositionsPlugin* spos_plgn;			///< source hit positions
	std::vector<SourceHitsPlugin*> srcPlugins;	///< plugins for each source to analyze
	MWPCGainPlugin* mwpcgain_plgn;				///< wirechamber gain measurement
};


/// re-generate source fits / plots
void reSource(RunNum rn);

/// upload sources from run log
void uploadRunSources(const std::string& rlogname = "UCNA Run Log.txt");










#endif
