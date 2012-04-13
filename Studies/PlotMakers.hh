#ifndef PLOTMAKERS_HH
#define PLOTMAKERS_HH 1

#include <TCanvas.h>
#include "EnergyCalibrator.hh"
#include "Source.hh"
#include "G4toPMT.hh"
#include <vector>
#include <string>

/// plot run-by-run GMS corrections, output data
void plotGMScorrections(const std::vector<RunNum>& runs, const std::string& foutPath = "../PostPlots/");

/// generate position map plots from a PositioningCorrector
void etaPlot(OutputManager& OM, PositioningCorrector* P, bool normalize = true, float axisRange = 2.5);

/// generate position map gradient plots from a PositioningCorrector
void etaGradPlot(OutputManager& OM, PositioningCorrector* P);


/// dump position map data to a file
void dumpPosmap(std::string basepath, unsigned int pnum);

/// generate nPE plots from a PMTCalibrator
void npePlot(OutputManager& OM, PMTCalibrator* PCal, float e0 = 1000, float s0 = 0.5, bool dumbsum = false);

/// extract Evis<->Etrue info from contiunuum spectrum simulation
void SimSpectrumInfo(Sim2PMT& S, OutputManager& OM);

/// generate a file with spectrum correction factors
void makeCorrectionsFile(const std::string& fout);

/// Class for generating position plots
class PosPlotter {
public:
	/// constructor
	PosPlotter(OutputManager* O): rscale(1.2), OM(O) {}
	/// plot number of PE
	void npePlot(PMTCalibrator* PCal);
	
	float rscale; //< extra radius to plot beyond edge of measured area
	
protected:
	OutputManager* OM;
	float x,y;
	int nx,ny;
	float r0;
	TH2* hCurrent;
	void startScan(TH2* h);
	bool nextPoint();
};


#endif
