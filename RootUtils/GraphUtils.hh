#ifndef GRAPHUTILS_HH
#define GRAPHUTILS_HH

#include "QFile.hh"
#include "Enums.hh"
#include <TF1.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TH3.h>
#include <TGraph.h>
#include <TProfile.h>
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

/// convert a TProfile to a TGraph
TGraphErrors* TProf2TGraph(const TProfile& P, unsigned int minpts = 0);

/// make cumulative histogram
TH1F* cumulativeHist(const TH1F& h, bool normalize = false);

/// inplace cumulative sum of TGraph, with switch to accumulate from top end
void makeCumulative(TGraph& g, bool fromFirst = true);


/// invert a TGraph
TGraph* invertGraph(const TGraph* g);

/// combine a list of TGraphs
TGraph* combine_graphs(const std::vector<TGraph*> gs);

/// merge a list of TGrapherrors into a larger TGraphError, offsetting x values by given offsets
TGraphErrors* merge_plots(const std::vector<TGraphErrors*>& pin, const std::vector<int>& toffset);

/// draw several graphs on the same canvas
void drawTogether(std::vector<TGraphErrors*>& gs, float ymin, float ymax, TCanvas* C, const char* outname, const char* graphTitle);

/// find transform curve to match two histograms
TGraph* matchHistoShapes(const TH1F& h1, const TH1F& h2);

/// scale a TGraphErrors
void scale(TGraphErrors& tg, float s);

/// accumulate TGraphErrors
void accumPoints(TGraphErrors& a, const TGraphErrors& b, bool errorWeight = true, bool yonly = false);

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
	
/// slice a TH3 into a stack of TH2Fs
std::vector<TH2F*> sliceTH3(const TH3& h3, AxisDirection d = Z_DIRECTION);

/// slice a TH2 into a stack of TH1Fs
std::vector<TH1F*> sliceTH2(const TH2& h2, AxisDirection d, bool includeOverflow = false);

/// split a list of counts into approximately equal segments
std::vector<unsigned int> equipartition(const std::vector<float>& elems, unsigned int n);

/// project a TH2 onto an arbitrary line
TH1* projectTH2(const TH2& h, double nb, double cx, double cy);

/// get plot indicating optimal separation point between two histograms
TH1* histsep(const TH1& h1, const TH1& h2);

/// get plot indicating optimal separation point between two TGraphs; note, converts graphs to cumulative form
TGraph* graphsep(TGraph& g1, TGraph& g2, unsigned int npts = 400);


/// calculate optimum dividing point and overlap error between two histograms
void histoverlap(const TH1& h1, const TH1& h2, double& xdiv, double& o);


/// Quantiles (inverse CDF) distribution from a TF1
/// based on ROOT's TF1::GetQuantiles(...) function
class TF1_Quantiles {
public:
	/// constructor
	TF1_Quantiles(TF1& f);
	/// return quantile for 0 <= p <= 1
	double eval(double p) const;
	
protected:

	const unsigned int npx;
	const Double_t xMin;
	const Double_t xMax;
	const Double_t dx;
	TArrayD integral;
	TArrayD alpha;
	TArrayD beta;
	TArrayD gamma;
};

	
#endif
