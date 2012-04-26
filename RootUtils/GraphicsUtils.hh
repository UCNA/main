#ifndef GRAPHICSUTILS_HH
#define GRAPHICSUTILS_HH 1

#include <math.h>
#include <TVirtualPad.h>
#include "Types.hh"
#include "Source.hh"
#include "SectorCutter.hh"
#include <vector>
#include <TH1.h>
#include <cfloat>

/// draw several histograms simultaneously; return max histogram height
double drawSimulHistos(std::vector<TH1*>& hists, const std::string& opt = "", const std::string& newTitle = "DEFAULT");
/// draw a pair of histograms (red and blue by default)
void drawHistoPair(TH1* hRed, TH1* hBlue, const std::string& opt = "", Int_t c1 = 2, Int_t c2 = 4);

/// 6" diameter, projected back
#define wirechamber_window_radius (3*25.4*sqrt(0.6))	
/// 4.590" diameter
#define collimator_radius (4.950*25.4/2)			
/// 5" diameter
#define decay_tube_radius (5*25.4/2)		

/// draw detector fiducial cuts at 50mm (decay trap radius) and 42.5mm (useful events radius)
void drawFiducialCuts(Int_t color = 6);
/// draw red ellipse around, e.g., a source's position
void drawEllipseCut(Source E, Float_t nSigma, std::string label = "");
/// draw vertical line marker
void drawVLine(Float_t x, TVirtualPad* C, Int_t color = 4);
/// draw shaded rectangle marker
void drawExcludedRegion(Float_t x0, Float_t x1, TCanvas* C, Int_t color = 4, Int_t fill = 1001);
/// draw sectors to current canvas
void drawSectors(const SectorCutter& S, int color = 2);
/// label SectorCutter sectors on current canvas
void labelSectors(const SectorCutter& S, int color = 2);

/// set up grayscale figures color palette
void makeGrayscalepalette();

#endif

