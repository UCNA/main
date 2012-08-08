#ifndef SIMASYMMETRYANALYZER_HH
#define SIMASYMMETRYANALYZER_HH 1

#include "OctetAnalyzer.hh"

/// plugin with extra plots for simulated runs
class SimAsymmetryAnalyzer: public OctetAnalyzerPlugin {
public:
	/// constructor
	SimAsymmetryAnalyzer(OctetAnalyzer* OA);
	
	/// fill from scan data point
	virtual void fillCoreHists(ProcessedDataScanner& PDS, double weight);
	/// output plot generation
	virtual void makePlots();
	
	quadHists* qBCT[2][TYPE_IV_EVENT+1];		//< average beta cos theta TProfiles by [side][event type]
	quadHists* qWrongSide[2][TYPE_III_EVENT+1];	//< Energy spectra of events ID'd on wrong side by [side][type]
	quadHists* qMissedSpectrum;					//< energy spectrum of missed events
};

#endif
