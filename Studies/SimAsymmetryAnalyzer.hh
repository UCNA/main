#ifndef SIMASYMMETRYANALYZER_HH
#define SIMASYMMETRYANALYZER_HH 1

#include "AsymmetryAnalyzer.hh"

/// primary octet data analysis class
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
	
	quadHists qBCT[2][TYPE_IV_EVENT+1];	//< average beta cos theta TProfiles by [side][event type]
	
protected:
	
	/// fill from scan data point
	virtual void fillCoreHists(ProcessedDataScanner& PDS, double weight);

	TProfile* pBCT[2][TYPE_IV_EVENT+1];	//< fill point for beta cos theta profile
};

#endif
