#ifndef ANODEPOSITIONANALYZER_HH
#define ANODEPOSITIONANALYZER_HH 1

#include "PositionBinnedAnalyzer.hh"

class AnodePositionAnalyzer: public PositionBinnedAnalyzer {
public:
	/// constructor
	AnodePositionAnalyzer(RunAccumulator* RA, unsigned int nr = 0);
	
	/// process a data point into position histograms
	virtual void fillCoreHists(ProcessedDataScanner& PDS, double weight);
	/// overall spectrum info
	virtual void calculateResults();
	
	std::vector<fgbgPair*> sectAnode[2];	//< anode by position on each side
};

#endif
