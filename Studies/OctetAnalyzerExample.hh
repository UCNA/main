#ifndef OCTETANALYZEREXAMPLE_HH
#define OCTETANALYZEREXAMPLE_HH 1

#include "OctetAnalyzer.hh"
#include "RunAccumulator.hh"
#include "PathUtils.hh"

/// minimalist example subclass of OctetAnalyzer: generates super-ratio and super-sum of wirechamber energy spectra
class OctetAnalyzerExample: public OctetAnalyzer {
public:
	/// constructor
	OctetAnalyzerExample(OutputManager* pnt, const std::string& nm, const std::string& inflname = "");
	/// destructor... you probably won't need to destruct anything
	~OctetAnalyzerExample() {}
	
	/// cloning generator: just return another of the same subclass (with any settings you want to preserve)
	virtual SegmentSaver* makeAnalyzer(const std::string& nm,
									   const std::string& inflname) { return new OctetAnalyzerExample(this,nm,inflname); }
	
	/// location of already-processed data (after first run) for errorbar estimation
	virtual std::string estimatorHistoLocation() const { return OctetAnalyzerExample::processedLocation; }
	static std::string processedLocation;	//< set location here for already-processed files
	
	
	/// fill from scan data point
	virtual void fillCoreHists(ProcessedDataScanner& PDS, double weight);
	/// calculate super-ratio asymmetry from anode spectra
	virtual void calculateResults();
	/// output plot generation
	virtual void makePlots();
	/// MC/Data comparison routine
	/// NOTE: this MUST NOT change the contents of saved histograms (calculated ones are OK)
	virtual void compareMCtoData(RunAccumulator& OAdata);
	
	quadHists qAnodeSpectrum[2];	//< set of histograms for extracting anode spectrum on each side
	TH1F* hAnodeSpectrum[2];		//< convenient pointer for currently active histogram
	TH1F* hAnodeSR;					//< super-ratio asymmetry of anode data
	TH1F* hAnodeSS;					//< super-sum of anode data
};

#endif
