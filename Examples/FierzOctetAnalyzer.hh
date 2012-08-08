#ifndef FIERZOCTETANALYZER_HH
#define FIERZOCTETANALYZER_HH

#include "OctetAnalyzer.hh"
#include "RunAccumulator.hh"
#include "PMTGenerator.hh"
#include "TGraphAsymmErrors.h"

using namespace std;

/// Fierz ratio extractor plugin
class FierzAnalyzerPlugin: public OctetAnalyzerPlugin {
public:
	/// constructor
	FierzAnalyzerPlugin(OctetAnalyzer* OA);
	
	/// fill from scan data point
	virtual void fillCoreHists(ProcessedDataScanner& PDS, double weight);
	/// calculate offset info
	virtual void calculateResults();
	/// output plot generation
	virtual void makePlots();
	/// MC/Data comparison routine
	/// NOTE: this MUST NOT change the contents of saved histograms (calculated ones are OK)
	virtual void compareMCtoData(AnalyzerPlugin* AP);
	
	quadHists* qFullEnergySpectrum[2];		//< set of histograms for extracting anode spectrum on each side
	TH1F* hFullEnergySR;					//< super-ratio asymmetry of anode data
	TH1F* hFullEnergySS;					//< super-sum of anode data
	
	PMTGenerator PGen;						//< event simulator for threshold extraction
	fgbgPair* pTriggerThreshold[2][2];		//< trigger threshold extraction histograms for each [side][all/trig]
	TGraphAsymmErrors* gTrigCurve[2];		//< extracted trigger efficiency curve
};

/// Octet analyzer using Fierz plugin
class FierzOctetAnalyzer: public OctetAnalyzer {
public:
	/// constructor
	FierzOctetAnalyzer(OutputManager* pnt, const string& nm, const string& inflname = "");
	
	/// cloning generator: just return another of the same subclass (with any settings you want to preserve)
	virtual SegmentSaver* makeAnalyzer(const string& nm,
									   const string& inflname) { 
		return new FierzOctetAnalyzer(this,nm,inflname); 
	}
	
	/// location of already-processed data (after first run) for errorbar estimation
	virtual string estimatorHistoLocation() const { return FierzOctetAnalyzer::processedLocation; }
	static string processedLocation;	//< set location here for already-processed files
	
	OctetAnalyzerPlugin* myFierz;
};

#endif
