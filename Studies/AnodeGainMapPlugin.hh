#ifndef ANODEPOSITIONANALYZER_HH
#define ANODEPOSITIONANALYZER_HH 1

#include "PositionBinnedPlugin.hh"

class AnodeGainMapPlugin: public PositionBinnedPlugin {
public:
	/// constructor
	AnodeGainMapPlugin(RunAccumulator* RA, unsigned int nr);
	
	/// process a data point into position histograms
	virtual void fillCoreHists(ProcessedDataScanner& PDS, double weight);
	/// generate and upload anode positioning map
	void genAnodePosmap();
	
	std::vector<fgbgPair*> sectAnode[2];	//< anode by position on each side
};

#endif
