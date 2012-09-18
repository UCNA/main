#ifndef BGDECAYANALYZER_HH
#define BGDECAYANALYZER_HH 1

#include "OctetAnalyzer.hh"

class BGDecayAnalyzer: public OctetAnalyzerPlugin {
public:
	/// constructor
	BGDecayAnalyzer(OctetAnalyzer* OA);
	/// fill from scan data point
	virtual void fillCoreHists(ProcessedDataScanner& PDS, double weight);
	
	quadHists* qBGDecay[2];		//< 5min E vs. time plot to see decaying BG components
};




#endif
