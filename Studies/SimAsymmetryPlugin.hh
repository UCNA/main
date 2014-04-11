#ifndef SIMASYMMETRYANALYZER_HH
#define SIMASYMMETRYANALYZER_HH

#include "OctetAnalyzer.hh"
#include "AsymmetryPlugin.hh"

/// plugin with extra plots for simulated runs
class SimAsymmetryPlugin: public OctetAnalyzerPlugin {
public:
	/// constructor
	SimAsymmetryPlugin(OctetAnalyzer* OA);
	
	/// fill from scan data point
	virtual void fillCoreHists(ProcessedDataScanner& PDS, double weight);
	/// output plot generation
	virtual void makePlots();
	/// calculate corrections given actual data; return asymmetry correction stages
	std::vector<TH1*> calculateCorrections(AsymmetryPlugin& Adat, AsymmetryPlugin& Asim);
	/// "old style" corrections calculation
	std::vector<TH1*> calculateCorrectionsOld(AsymmetryPlugin& Asim);
	
	quadHists* qBCT[2][TYPE_IV_EVENT+1];		///< average beta cos theta TProfiles by [side][event type]
	quadHists* qCosth[2][TYPE_IV_EVENT+1];		///< average cos theta TProfiles by [side][type]
	quadHists* qBeta[2][TYPE_IV_EVENT+1];		///< average energy TProfiles by [side][type]
	quadHists* qAsymAcc[2][TYPE_IV_EVENT+1];	///< average asymmetry acceptance
	quadHists* qWrongSide[2][TYPE_III_EVENT+1];	///< Energy spectra of events ID'd on wrong side by [side][type]
	quadHists* qMissedSpectrum;					///< energy spectrum of missed events
};

#endif
