#ifndef HIGHENERGYEXCESS_HH
#define HIGHENERGYEXCESS_HH 1

#include "OctetAnalyzer.hh"
#include "MuonAnalyzer.hh"
#include "AnalysisDB.hh"

/// fit beyond-endpoint background subtraction
void fitHighEnergyExcess(QFile& qOut, quadHists* qh, double e0, double e1);

/// energy spectra and asymmetry analysis class
class HighEnergyExcessAnalyzer: public OctetAnalyzerPlugin {
public:
	/// constructor
	HighEnergyExcessAnalyzer(OctetAnalyzer* OA);
	
	/// fill from scan data point
	virtual void fillCoreHists(ProcessedDataScanner& PDS, double weight);
	/// calculate super-ratio asymmetry from anode spectra
	virtual void calculateResults();
	/// output plot generation
	virtual void makePlots();
	
	quadHists* qExcessSpectra[2];	//< excess high-energy event spectrum
	quadHists* qExcessGamma[2];		//< excess high-energy gamma event spectrum
	quadHists* qExcessr2[2];		//< radius^2 distribution of >1keV "excess" Type 0,I events
	quadHists* qExcessTheta[2];		//< angular distribution of >1keV excess events
};

#endif
