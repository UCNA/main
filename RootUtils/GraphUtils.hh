#ifndef GRAPHUTILS_HH
#define GRAPHUTILS_HH 1

#include "QFile.hh"
#include "Enums.hh"
#include <TF1.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TH3.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <vector>

/// convert TH1* to Stringmap
Stringmap histoToStringmap(const TH1* h);
/// convert Stringmap to TH1F*
TH1F* stringmapToTH1F(const Stringmap& m);
/// convert a TGraph to a Stringmap
Stringmap graphToStringmap(const TGraph& g);

/// convert a histogram to a TGraph
TGraphErrors* TH1toTGraph(const TH1& h);

/// make cumulative histogram
TH1F* cumulativeHist(const TH1F& h, bool normalize = false);

/// invert a TGraph
TGraph* invertGraph(const TGraph* g);

/// merge a list of TGrapherrors into a larger TGraphError, offsetting x values by given offsets
TGraphErrors* merge_plots(const std::vector<TGraphErrors*>& pin, const std::vector<int>& toffset);

/// draw several graphs on the same canvas
void drawTogether(std::vector<TGraphErrors*>& gs, float ymin, float ymax, TCanvas* C, const char* outname, const char* graphTitle);

/// find transform curve to match two histograms
TGraph* matchHistoShapes(const TH1F& h1, const TH1F& h2);

/// scale a TGraphErrors
void scale(TGraphErrors& tg, float s);

/// generate derivative graph
TGraph* derivative(TGraph& g);

/// transform axis on a TGraph, optionally using jacobean to preserve integral
void transformAxis(TGraph& g, TGraph& T, bool useJacobean = true);

/// interpolate a TGraphErrors to a new grid spacing ~dx
TGraphErrors* interpolate(TGraphErrors& tg, float dx);

/// get inverse CDF x: CDF(x)=p for histogram
double invCDF(TH1* h, double p);

/// check histogram for NaNs
void fixNaNs(TH1* h);

/// replacement for missing TH2::FitSlicesY in older version of ROOT
std::vector<TH1D*> replaceFitSlicesY(TH2* h, TF1* f1=NULL);
	
/// slice a TH3 into a stack of TH2Fs
std::vector<TH2F*> sliceTH3(const TH3& h3);

/// slice a TH2 into a stack of TH1Fs
std::vector<TH1F*> sliceTH2(const TH2& h2, AxisDirection d, bool includeOverflow = false);

/// split a list of counts into approximately equal segments
std::vector<unsigned int> equipartition(const std::vector<float>& elems, unsigned int n);

#endif
