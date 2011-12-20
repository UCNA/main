#ifndef KURIEFITTER_HH
#define KURIEFITTER_HH 1

#include "PMTGenerator.hh"
#include <TGraphErrors.h>
#include <TH1F.h>
#include "GraphUtils.hh"
#include "BetaSpectrum.hh"

/// use Kurie plot to estimate endpoint given initial endpoint estimate (and optional linearity correction)
/// input is linearized but uncalibrated spectrum and endpoint guess (optional output graph, normalization at 500, and error estimate)
float_err kuriePlotter(TH1* spectrum, float endpoint, TGraphErrors** tgout = NULL, float targetEP = neutronBetaEp, float fitStart = 250, float fitEnd = 700);

/// iterate Kurie plot until convergence is reached
float_err kurieIterator(TH1* spectrum, float iguess, TGraphErrors** tgout = NULL, float targetEP = neutronBetaEp,
						float fitStart = 250, float fitEnd = 700, unsigned int nmax = 50, float etol = 0.02 );


#endif
