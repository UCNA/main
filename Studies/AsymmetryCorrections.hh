#ifndef ASYMMETRYCORRECTIONS_HH
#define ASYMMETRYCORRECTIONS_HH

#include "SimAsymmetryPlugin.hh"
#include "BetaDecayAnalyzer.hh"
#include "OutputManager.hh"
#include "PathUtils.hh"
#include <TGraph.h>
#include <TGraphErrors.h>
#include <string>

/// asymmetry correction/uncertainty read from file
class AsymCorr {
public:
	/// constructor with name to look up in Aux/Corrections
	AsymCorr(const std::string& nm): name(nm) {}
	/// destructor
	virtual ~AsymCorr() {}
	
	/// get correction at energy
	virtual double getCor(double e) = 0;
	/// get uncertainty at energy
	virtual double getUnc(double e) = 0;
	
	std::string name;	//< correction name
};

class AsymCorrFile: public AsymCorr {
public:
	/// constructor
	AsymCorrFile(const std::string& nm, const std::string& basePath = getEnvSafe("UCNA_AUX")+"/Corrections/");
	/// get correction at energy
	virtual double getCor(double KE) { return gCor.Eval(KE); }
	/// get uncertainty at energy
	virtual double getUnc(double KE) { return gUnc.Eval(KE); }
protected:
	TGraph gCor;		//< amount of correction
	TGraph gUnc;		//< uncertainty
};

/// energy independent asymmetry correction
class ConstAsymCorr: public AsymCorr {
public:
	/// constructor
	ConstAsymCorr(const std::string& nm, double c, double e): AsymCorr(nm), corr(c), uncert(e) {}
	/// get correction at energy
	virtual double getCor(double) { return corr; }
	/// get uncertainty at energy
	virtual double getUnc(double) { return uncert; }
protected:
	double corr;
	double uncert;
};

/// do corrections to asymmetry
void doFullCorrections(AsymmetryPlugin& AA, OutputManager& OM, std::string mcBase = "");

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
	/// energy reconstruction errors
	void eLinearityTable(unsigned int yr = 2010);
	/// muon veto efficiency change errors
	void muonVetoEfficTable(double delta);
	/// uniform efficiency shifts tables (e.g. deadtime, veto accidentals)
	void efficShiftTable(double delta);
	///constant neutron generated background (in Hz/keV)
	void NGBGTable(double EScale, double dEScale, double WScale, double dWScale, double dAFPfrac);
	
	/// energy reconstruction uncertainty envelope
	double energyErrorEnvelope(double e, unsigned int year = 2010) const;
	
protected:
	OutputManager OM;		//< unused OutputManager
	BetaDecayAnalyzer Adat;	//< data for error estimation
	TGraphErrors* S[2][2];	//< observed energy spectra for [side][afp] as TGraphs
};

/// calculate MC-based corrections given data, MC filenames
void calcMCCorrs(OutputManager& OM, const std::string& datin, const std::string& simin, const std::string& outDir = "", bool oldCorr = false);

/// comparison between two MCs
void compareMCs(OutputManager& OM, const std::string& sim0, const std::string& sim1, const std::string& fOut = "");

#endif
