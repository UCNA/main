#ifndef SIMASYMMETRYANALYZER_HH
#define SIMASYMMETRYANALYZER_HH 1

#include "AsymmetryAnalyzer.hh"

/// Octet data analysis with extra plots for simulated runs
class SimAsymmetryAnalyzer: public AsymmetryAnalyzer {
public:
	/// constructor
	SimAsymmetryAnalyzer(OutputManager* pnt, const std::string& nm, const std::string& inflname = "");
	
	/// cloning generator
	virtual SegmentSaver* makeAnalyzer(const std::string& nm,const std::string& inflname);
	
	/// MC/data comparison
	void compareMCtoData(RunAccumulator& OAdata);
	
	/// output plot generation
	virtual void makePlots();
	
	quadHists* qBCT[2][TYPE_IV_EVENT+1];		//< average beta cos theta TProfiles by [side][event type]
	quadHists* qWrongSide[2][TYPE_IV_EVENT+1];	//< Energy spectra of events ID'd on wrong side by [side][type]
	quadHists* qMissedSpectrum;					//< energy spectrum of missed events
	
protected:
	
	/// fill from scan data point
	virtual void fillCoreHists(ProcessedDataScanner& PDS, double weight);
};

#endif
