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
void NGBGSpectra(std::string simName);

/// determine Type II/III anode cuts
void separate23(std::string datname);

/// test for BG subtraction statistics errors
void lowStatsTest();

/// redo position map from Xe anode data
void refitXeAnode(std::string datname);

/// generate asymmetry spectra for each analysis choice
void calcAnalysisChoices(OutputManager& OM, const std::string& inflname);

#endif
