#ifndef HIGHENERGYEXCESS_HH
#define HIGHENERGYEXCESS_HH 1

#include "OctetAnalyzer.hh"
#include "MuonPlugin.hh"
#include "AnalysisDB.hh"
#include "AsymmetryPlugin.hh"

/// fit beyond-endpoint background subtraction
void fitHighEnergyExcess(QFile& qOut, quadHists* qh, double e0, double e1);

/// energy spectra and asymmetry analysis class
class HighEnergyExcessPlugin: public OctetAnalyzerPlugin {
public:
	/// constructor
	HighEnergyExcessPlugin(OctetAnalyzer* OA);
	
	/// fill from scan data point
	virtual void fillCoreHists(ProcessedDataScanner& PDS, double weight);
	/// calculate super-ratio asymmetry from anode spectra
	virtual void calculateResults();
	/// output plot generation
	virtual void makePlots();
	/// MC output plot generation
	virtual void compareMCtoData(AnalyzerPlugin* AP);
	
	quadHists* qExcessSpectra[2];	//< excess high-energy event spectrum
	quadHists* qExcessGamma[2];		//< excess high-energy gamma event spectrum
	quadHists* qExcessr2[2];		//< radius^2 distribution of >1keV "excess" Type 0,I events
	quadHists* qExcessTheta[2];		//< angular distribution of >1keV excess events
};

class NGBGAnalyzer: public OctetAnalyzer {
public:
	/// constructor
	NGBGAnalyzer(OutputManager* pnt, const std::string& nm = "BetaDecayAnalyzer", const std::string& inflName = "");
	/// create a new instance of this object (cloning self settings) for given directory
	virtual SegmentSaver* makeAnalyzer(const std::string& nm, const std::string& inflname) { return new NGBGAnalyzer(this,nm,inflname); }
	
	HighEnergyExcessPlugin* myHEE;	//< high energy excess events
	AsymmetryPlugin* myAsym;		//< asymmetry plugin
};

/// process simulation for neutron generated background
void NGBGSpectra(std::string simName);

#endif
