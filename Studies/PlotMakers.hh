#ifndef PLOTMAKERS_HH
#define PLOTMAKERS_HH 1

#include <TCanvas.h>
#include "EnergyCalibrator.hh"
#include "Source.hh"
#include "G4toPMT.hh"
#include "NuclEvtGen.hh"
#include "OutputManager.hh"
#include "BetaSpectrum.hh"
#include "BetaDecayAnalyzer.hh"
#include <vector>
#include <string>

/// Class for generating position plots
class PosPlotter {
public:
	/// constructor
	PosPlotter(OutputManager* O);
	/// plot number of PE
	void npePlot(PMTCalibrator* PCal);
	/// plot normalized gradient of nPE
	void npeGradPlot(PMTCalibrator* PCal);
	/// plot light transport
	void etaPlot(PositioningCorrector* P, double axisRange = 2.0);
	
	float rscale; 		//< extra radius to plot beyond edge of measured area
	unsigned int nbin;	//< number of position bins
	
protected:
	/// generate new position histogram
	TH2F* makeHisto(const std::string& nm, const std::string& title);
	
	OutputManager* OM;
	float x,y;
	int nx,ny;
	float r0;
	TH2* hCurrent;
	void startScan(TH2* h);
	bool nextPoint();
};


/// plot run-by-run GMS corrections, output data
void plotGMScorrections(const std::vector<RunNum>& runs, const std::string& foutPath = "../PostPlots/");

/// dump position map data to a file
void dumpPosmap(std::string basepath, unsigned int pnum);

/// extract Evis<->Etrue info from contiunuum spectrum simulation
void SimSpectrumInfo(Sim2PMT& S, OutputManager& OM);

/// generate a file with spectrum correction factors
void makeCorrectionsFile(int A = 1, int Z = 1, double Endpt = neutronBetaEp);

/// make plots for simulated source spectrum
void showSimSpectrum(const std::string& nm, OutputManager& OM, NucDecayLibrary& NDL, PMTCalibrator& PCal);
	
/// Various Xenon spectra
void compareXenonSpectra();

/// Decompose xenon spectrum for given run number
void decomposeXenon(RunNum rn, bool includeFast = false);

/// process simulation for neutron generated background
void NGBGSpectra(std::string datname);

/// determine Type II/III anode cuts
void separate23(std::string datname);

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


/// test for BG subtraction statistics errors
void lowStatsTest();

/// redo position map from Xe anode data
void refitXeAnode(std::string datname);

/// generate asymmetry spectra for each analysis choice
void calcAnalysisChoices(OutputManager& OM, const std::string& inflname);

#endif
