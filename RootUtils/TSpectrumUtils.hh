#ifndef TSPECTRUMUTILS_HH
#define TSPECTRUMUTILS_HH 1

#include <vector>
#include <TH1F.h>
#include <TPad.h>
#include <TSpectrum.h>
#include "MultiGaus.hh"
#include "SpectrumPeak.hh"
#include "Types.hh"
#include <float.h>
#include <algorithm>

/// simple struct for peaks identified by TSpectrum
struct TSpectrumPeak {
	float x;	//< peak position
	float y;	//< peak height
};

/// more convenient wrapper for TSpectrumSearch
std::vector<TSpectrumPeak> tspectrumSearch(TH1F* hin, float sigma = 2.0, float thresh = 0.01);

/// perform a TSpectrum peak search on the input histogram for peaks of width sigma bins; report results up to specified threshold below top peak
std::vector<float> tspectrumSearch(TH1F* hin, TH1F** hout, float sigma, float thresh);

/// pre-fit using a TSpectrum to estimate peak positions
std::vector<SpectrumPeak> tspectrumPrefit(TH1F* indat, float searchsigma, const std::vector<SpectrumPeak>& expectedPeaks,
										  TH1F*& htout, float pkMin = -FLT_MAX, float pkMax = FLT_MAX);

/// generate fitter for fitting multi-peak spectum
MultiGaus multiPeakFitter(TH1F* indat, const std::vector<SpectrumPeak>& expectedPeaks, float nSigma = 1.5);

/// locate peaks in multi-peak spectrum based on initial guesses
std::vector<SpectrumPeak> fancyMultiFit(TH1F* indat, float searchsigma, const std::vector<SpectrumPeak>& expectedPeaks,
										bool bgsubtract = false, const std::string& drawName = "", float nSigma = 1.5, float pkMin = -FLT_MAX, float pkMax = FLT_MAX);

#endif
