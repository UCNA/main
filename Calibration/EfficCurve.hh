#ifndef EFFICCURVE_HH
#define EFFICCURVE_HH 1

#include "OutputManager.hh"
#include <TH1F.h>
#include <TGraphAsymmErrors.h>

class EfficCurve: public OutputManager {
public:
	/// constructor
	EfficCurve(std::string nm="", OutputManager* prnt=NULL): OutputManager(nm,prnt), gEffic(NULL) {}
	/// calculate efficiency curves from input hitograms
	virtual void genEffic(TH1F* hAll, TH1F* hTrig, bool adcChan = false);
	/// return efficiency at given point
	virtual double effic(double x) const;
	/// invert efficiency effects on  spectrum histogram
	virtual void invertEffic(TH1F* hIn, float th=0.25);
	/// get 50% trigger threshold
	virtual float getThreshold() const { return params[0]; }
	
	double params[4];			//< efficiency curve fit parameters
	TGraphAsymmErrors* gEffic;	//< full curve as TGraph
};

/// CDF for poisson function (for trigger efficiency fits)
Double_t poiscdf(const Double_t *x, const Double_t *par);
/// Fancier trigger efficiency model
Double_t fancyfish(const Double_t *x, const Double_t *par);

#endif
