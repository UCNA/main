#ifndef ASYMMETRYCORRECTIONS_HH
#define ASYMMETRYCORRECTIONS_HH 1

#include "SimAsymmetryAnalyzer.hh"
#include "BetaDecayAnalyzer.hh"
#include "OutputManager.hh"
#include <TGraph.h>
#include <TGraphErrors.h>
#include <string>

/// asymmetry correction/uncertainty read from file
class AsymCorr {
public:
	/// constructor with name to look up in Aux/Corrections
	AsymCorr(const std::string& nm);
	
	std::string name;	//< correction name
	TGraph gCor;		//< amount of correction
	TGraph gUnc;		//< uncertainty
};

/// do corrections to asymmetry
void doFullCorrections(AsymmetryAnalyzer& AA, OutputManager& OM);

/// systematic errors table generator
class ErrTables {
public:
	/// constructor
	ErrTables(const std::string& datset = "OctetAsym_Offic");
	/// destructor
	~ErrTables();
	/// experimental super-ratio R
	double getRexp(double e) const;
	/// experimental asymmetry at given energy
	double getAexp(double e) const;
	/// asymmetry from superratio
	static double AofR(double R);
	/// gain fluctuations errors
	void gainfluctsTable(double delta);
	/// pedestal fluctuations errors
	void pedShiftsTable(double delta);
	/// muon veto efficiency change errors
	void muonVetoEfficTable(double delta);
	/// uniform efficiency shifts tables (e.g. deadtime, veto accidentals)
	void efficShiftTable(double delta);
	///constant neutron generated background (in Hz/keV)
	void NGBGTable(double EScale, double dEScale, double WScale, double dWScale, double dAFPfrac);
	
protected:
	OutputManager OM;		//< unused OutputManager
	BetaDecayAnalyzer Adat;	//< data for error estimation
	TGraphErrors* S[2][2];	//< observed energy spectra for [side][afp] as TGraphs
};

/// calculate MC-based corrections given data, MC filenames
void calcMCCorrs(OutputManager& OM, const std::string& datin, const std::string& simin, bool writeAux = false, bool oldCorr = false);

/// comparison between two MCs
void compareMCs(OutputManager& OM, const std::string& sim0, const std::string& sim1, const std::string& fOut = "");

#endif
