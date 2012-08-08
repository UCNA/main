#ifndef OCTETANALYZEREXAMPLE_HH
#define OCTETANALYZEREXAMPLE_HH 1

#include "OctetAnalyzer.hh"
#include "RunAccumulator.hh"
#include "PathUtils.hh"

/// example of an analyzer plugin
class ExampleAnalyzerPlugin: public OctetAnalyzerPlugin {
public:
	/// constructor
	ExampleAnalyzerPlugin(OctetAnalyzer* OA);
	
	/// fill from scan data point
	virtual void fillCoreHists(ProcessedDataScanner& PDS, double weight);
	/// calculate super-ratio asymmetry from anode spectra
	virtual void calculateResults();
	/// output plot generation
	virtual void makePlots();
	/// MC/Data comparison routine
	/// NOTE: this MUST NOT change the contents of saved histograms (calculated ones are OK)
	virtual void compareMCtoData(AnalyzerPlugin* AP);
	
	quadHists* qAnodeSpectrum[2];	//< set of histograms for extracting anode spectrum on each side
	TH1F* hAnodeSR;					//< super-ratio asymmetry of anode data
	TH1F* hAnodeSS;					//< super-sum of anode data
};


/// Analyzer class using example analyzer plugin
class OctetAnalyzerExample: public OctetAnalyzer {
public:
	/// constructor
	OctetAnalyzerExample(OutputManager* pnt, const std::string& nm, const std::string& inflname = ""):
	OctetAnalyzer(pnt,nm,inflname) { addPlugin(myPlugin = new ExampleAnalyzerPlugin(this)); }
	
	/// cloning generator: just return another of the same subclass (with any settings you want to preserve)
	virtual SegmentSaver* makeAnalyzer(const std::string& nm,
									   const std::string& inflname) { return new OctetAnalyzerExample(this,nm,inflname); }
	
	/// location of already-processed data (after first run) for errorbar estimation
	virtual std::string estimatorHistoLocation() const { return OctetAnalyzerExample::processedLocation; }
	static std::string processedLocation;	//< set location here for already-processed files

	ExampleAnalyzerPlugin* myPlugin;	//< pointer to analyzer plugin
};

#endif
