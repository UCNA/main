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

/// extract data from line simulations for eQuenched->eTrue calibration curves
void ProcessLineSims(const std::string& baseDir, unsigned int nLines, float eMin, float eMax);

/// Fit Xe spectrum 915keV endpoint
void XeEndpointStudy();

/// generate a file with spectrum correction factors
void makeCorrectionsFile(const std::string& fout);

#endif
