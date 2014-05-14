#ifndef POSITIONANALYZER_HH
#define POSITIONANALYZER_HH

#include "OctetAnalyzer.hh"
#include "SectorCutter.hh"

/// event positions by type octet analyzer
class PositionsPlugin: public OctetAnalyzerPlugin {
public:
	/// constructor
	PositionsPlugin(OctetAnalyzer* OA);
	
	/// fill from scan data point
	virtual void fillCoreHists(ProcessedDataScanner& PDS, double weight);
	/// calculate offset info
	virtual void calculateResults();
	/// output plot generation
	virtual void makePlots();
	/// MC/data comparison plots
	virtual void compareMCtoData(AnalyzerPlugin* AP);
	
	SectorCutter offSects;						///< sector cutter for position binning data
	std::vector<fgbgPair*> poff[2];				///< East-West position offsets for [x/y direction] in each sector
	quadHists* qPositions[2][TYPE_III_EVENT+1];	///< event positions quad hists for [side][type]
	quadHists* qRadius2[2][TYPE_III_EVENT+1];	///< event radius^2 by [side][type]
	
	TH1* hSuperSumPos[TYPE_III_EVENT+1];		///< super-sum positions by event type
};

#endif
