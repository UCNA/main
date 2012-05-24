#ifndef MULTIGAUS_HH
#define MULTIGAUS_HH 1

#include <TF1.h>
#include <TH1F.h>
#include "Types.hh"

/// class for fitting multi-peak gaussians
class MultiGaus {
public:	
	
	/// correlated subpeaks specification
	struct corrPeak {
		unsigned int mainPeak;
		double relCenter;
		double relHeight;
		double relWidth;
	};
	
	/// constructor
	MultiGaus(unsigned int n, const char* name, float ns = 1.5): nSigma(ns), npks(n), iguess(new double[3*n]), myTF1(new TF1(name,this,0,0,3*n)) { }
	
	/// destructor
	~MultiGaus();	
	
	/// add correlated peak
	void addCorrelated(unsigned int n, double relCenter, double relHeight, double relWidth = 0);
	/// fill initial values array
	void setParameter(unsigned int n, double p);	
	/// get fit parameter
	double getParameter(unsigned int n) const;	
	/// get fit parameter error
	double getParError(unsigned int n) const;
	/// get parameter+error as float_err
	float_err getPar(unsigned int n) const;
	
	/// get TF1 with appropriate pre-set values
	TF1* getFitter();	
	
	/// fit a TH1F after initial peak centers/widths have been guessed; update inital guess
	void fit(TH1F* h, bool draw = true);
	
	/// gaussian evaluation function
	double operator() (double* x, double* par);
	
	float nSigma;				//< number of sigma peak width to fit
	
protected:
	const unsigned int npks;			//< number of peaks being fitted
	double* iguess;						//< inital guess at peak positions
	TF1* myTF1;							//< TF1 using this class as its fit function
	std::vector<corrPeak> corrPeaks;	//< correlated subpeaks
};

int iterGaus(TH1* h0, TF1* gf, unsigned int nit, float mu, float sigma, float nsigma = 1.5, float asym = 0);

#endif
