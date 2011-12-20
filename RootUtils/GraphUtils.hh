#ifndef GRAPHUTILS_HH
#define GRAPHUTILS_HH 1

#include "QFile.hh"
#include <TH1.h>
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

#endif
