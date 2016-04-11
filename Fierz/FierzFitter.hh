#ifndef FIERZ_FITTER
#define FIERZ_FITTER

// UCNA includes
#include "G4toPMT.hh"
#include "PenelopeToPMT.hh"
#include "CalDBSQL.hh"

// ROOT includes
#include <TH1.h>
#include <TLegend.h>
#include <TF1.h>
#include <TVirtualFitter.h>
#include <TList.h>
#include <TStyle.h>
#include <TApplication.h>

// c++ includes
#include <iostream>
#include <fstream>
#include <string>

// c includes
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

///
/// physical constants
///
const double    pi          = 3.1415926535;         /// pi
const double    alpha       = 1/137.035999;         /// fine-structure constant
const double    a2pi        = alpha/(2*pi);         /// alpha over 2pi
const double    m_e         = 510.9988484;          /// electron mass       (14)
const double    m_p         = 938272.046;           /// proton mass         (21)
const double    m_n         = 939565.378;           /// neutron mass        (21)
const double    mu_p        = 2.792846;             /// proton moment       (7)
const double    mu_n        = -1.9130427;           /// neutron moment      (5)
const double    muV         = mu_p - mu_n;          /// nucleon moment      (9)
const double    Q           = 782.344;              /// end point KE        (30)
const double    E0          = m_e + Q;              /// end point E         (30)
const double    lambda      = -1.27590;             /// gA/gV               (450)
const double    M           = 1 + 3*lambda*lambda;  /// matrix element
const double    I_0         = 1.63632;              /// 0th moment          (25)
const double    I_1         = 1.07017;              /// 1st moment          (15)
const double    x_1         = I_1/I_0;              /// first m/E moment    (9)


using namespace std;
#if 0
int output_histogram(string filename, TH1D* h, double ax, double ay)
{
	ofstream ofs;
	ofs.open(filename.c_str());
	for (int i = 1; i < h->GetNbinsX(); i++)
	{
		double x = ax * h->GetBinCenter(i);
		//double sx = h->GetBinWidth(i);
		double r = ay * h->GetBinContent(i);
		double sr = ay * h->GetBinError(i);
		ofs << x << '\t' << r << '\t' << sr << endl;
	}
	ofs.close();
	return 0;
}
#endif

struct STLhistogram {
    int bins;
    int detector;
    int spin;
    double min, max;
    TH1D* histogram;
    vector<double> energy;        
    vector<double> values;        
    vector<double> errors;        

    void init(int bins, double min, double max) {
        this->bins = bins;
        this->min = min;
        this->max = max;
        energy = vector<double>(bins);
	    values = vector<double>(bins);
	    errors = vector<double>(bins);
    }

    void fill(TH1D *histogram) {
        TAxis *axis  = histogram->GetXaxis();
        bins = histogram->GetNbinsX(); 
        for (int ix = 1; ix <= bins + 1; ix++)
        {
            double E = axis->GetBinCenter(ix);
            if (min < E and E < max)
            {
                double Y = histogram->GetBinContent(ix);
                double eY = histogram->GetBinError(ix);
                energy.push_back(E);
                values.push_back(Y);
                errors.push_back(eY);
            }
        }
    }
};

struct UCNAModel {
    int    bins;
    double min;
    double max;
    TH1D *raw[2][2];
    //TH1D *super_sum_histogram;
    //TH1D *super_ratio_histogram;

    //STLhistogram raw[2][2];
    STLhistogram super_ratio;
    STLhistogram super_sum;
    STLhistogram asymmetry;

    void init(int bins, double min, double max) {
        this->bins = bins;
        this->min = min;
        this->max = max;
        // TODO init TH1Ds
        for (int side = 0; side < 2; side++)
            for (int spin = 0; spin < 2; spin++) {
                //raw[side][spin].init(bins, min, max);
                raw[side][spin] = 0;
            }
        super_ratio.init(bins, min, max);
        super_sum.init(bins, min, max);
        asymmetry.init(bins, min, max);
   }
};


struct UCNAEvent {
    double EdepQ;
    double Edep;
    double MWPCEnergy;
    double ScintPos;
    double MWPCPos;
    double time;
    double primMomentum;
};

class UCNAFierzFitter {
  public: 
    double minBin;
    double maxBin;
    unsigned int nBins;
    /*
    TH1D *fierz_super_sum_histogram;
    TH1D *sm_super_sum_histogram;
    TH1D *fierz_histogram[2][2];
    TH1D *sm_histogram[2][2];
    */
    UCNAModel sm;
    UCNAModel fierz;
    UCNAModel data;

    UCNAFierzFitter(unsigned int _nBins, double _minBin, double _maxBin) {
        minBin = _minBin;
        maxBin = _maxBin;
        nBins = _nBins;
        fierz.super_sum.histogram = new TH1D("fierz_histogram", "", nBins, minBin, maxBin);
        sm.super_sum.histogram = new TH1D("standard_model_histogram", "", nBins, minBin, maxBin);
        for (int side = 0; side < 2; side++)
            for (int spin = 0; spin < 2; spin++) {
                fierz.histogram[side][spin] = new TH1D("fierz_super_sum", "", nBins, minBin, maxBin);
                sm.histogram[side][spin] = new TH1D("standard_model_super_sum", "", nBins, minBin, maxBin);
            }
    }

    /*
    double evaluate(double *x, double*p) {
        double rv = 0;
        rv += p[0] * sm_histogram->GetBinContent(sm_histogram->FindBin(x[0]));        
        rv += p[1] * fierz_histogram->GetBinContent(fierz_histogram->FindBin(x[0]));
        return rv;
    }
    */

    void normalizeHistogram(TH1D* hist) {
        hist->Scale(1/(hist->GetBinWidth(2)*hist->Integral()));
    }

    void normalizeHistogram(TH1D* hist, double min, double max) {
		int _min = hist->FindBin(min);
		int _max = hist->FindBin(max);
        hist->Scale(1/(hist->GetBinWidth(2)*hist->Integral(_min, _max)));
    }
};


/// beta spectrum with little b term
double fierz_beta_spectrum(const double *val, const double *par) 
{
	const double K = val[0];                    /// kinetic energy
	if (K <= 0 or K >= Q)
		return 0;                               /// zero outside range

	const double b = par[0];                    /// Fierz parameter
	const int n = par[1];                    	/// Fierz exponent
	const double E = K + m_e;                   /// electron energy
	const double e = Q - K;                     /// neutrino energy
	const double p = sqrt(E*E - m_e*m_e);       /// electron momentum
	const double x = pow(m_e/E,n);              /// Fierz term
	const double f = (1 + b*x)/(1 + b*x_1);     /// Fierz factor
	const double k = 1.3723803E-11/Q;           /// normalization factor
	const double P = k*p*e*e*E*f*x;             /// the output PDF value

	return P;
}



/// beta spectrum with expected x^-n and beta^m
double beta_spectrum(const double *val, const double *par) 
{
	const double K = val[0];                    	///< kinetic energy
	if (K <= 0 or K >= Q)
		return 0;                               	///< zero beyond endpoint

	const double m = par[0];                    	///< beta exponent
	const double n = par[1];                    	///< Fierz exponent
	const double E = K + m_e;                   	///< electron energy
	const double B = pow(1-m_e*m_e/E/E,(1+m)/2);  	///< beta power factor
	const double x = E / m_e;                   	///< reduced electron energy
	const double y = (Q - K) / m_e;             	///< reduced neutrino energy
	const double z = pow(x,2-n);          			///< Fierz power term
	const double k = 1.3723803E-11/Q;           	///< normalization factor
	const double P = k*B*z*y*y;             		///< the output PDF value

	return P;
}




double evaluate_expected_fierz(double m, double n, double min, double max, int integral_size = 1234) 
{
    TH1D *h1 = new TH1D("beta_spectrum_fierz", "Beta spectrum with Fierz term", integral_size, min, max);
    TH1D *h2 = new TH1D("beta_spectrum", "Beta Spectrum", integral_size, min, max);
	for (int i = 0; i < integral_size; i++)
	{
		double K = min + double(i)*(max-min)/integral_size;
		double par1[2] = {m, n};
		double par2[2] = {0, 0};
		double y1 = beta_spectrum(&K, par1);
		double y2 = beta_spectrum(&K, par2);
		h1->SetBinContent(K, y1);
		h2->SetBinContent(K, y2);
	}
	return h1->Integral(0, integral_size) / h2->Integral(0, integral_size);
}



double evaluate_expected_fierz(double min, double max, int integral_size = 1234) 
{
	return evaluate_expected_fierz(0, 1, min, max, integral_size);
}



void compute_fit(TH1D* histogram, TF1* fierz_fit) 
{
	// compute chi squared
    double chisq = fierz_fit->GetChisquare();
    double N = fierz_fit->GetNDF();
	char A_str[1024];
	sprintf(A_str, "A = %1.3f #pm %1.3f", fierz_fit->GetParameter(0), fierz_fit->GetParError(0));
	char b_str[1024];
	sprintf(b_str, "b = %1.3f #pm %1.3f", fierz_fit->GetParameter(1), fierz_fit->GetParError(1));
	char chisq_str[1024];
    printf("Chi^2 / ( N - 1) = %f / %f = %f\n", chisq, N-1, chisq/(N-1));
	sprintf(chisq_str, "#frac{#chi^{2}}{n-1} = %f", chisq/(N-1));

	// draw the ratio plot
	histogram->SetStats(0);
    histogram->Draw();

	// draw a legend on the plot
    TLegend* ratio_legend = new TLegend(0.3,0.85,0.6,0.65);
    ratio_legend->AddEntry(histogram, "Asymmetry data", "l");
    ratio_legend->AddEntry(fierz_fit, "fit", "l");
    ratio_legend->AddEntry((TObject*)0, A_str, "");
    ratio_legend->AddEntry((TObject*)0, b_str, "");
    ratio_legend->AddEntry((TObject*)0, chisq_str, "");
    ratio_legend->SetTextSize(0.03);
    ratio_legend->SetBorderSize(0);
    ratio_legend->Draw();
}





double asymmetry_fit_func(double *x, double *par)
{
	double A = par[0];
	double b = par[1];
	double E = x[0] + m_e;

	return A / (1 + b*m_e/E);
}



double fierzratio_fit_func(double *x, double *par)
{
	double b = par[0];
	double m_e_E = par[1];
	double E = x[0] + m_e;

	return (1 + b*m_e/E) / (1 + b*m_e_E);
}



double GetEntries(TH1D* histogram, double min, double max)
{
	double entries = histogram->GetEffectiveEntries();
	double part_int = histogram->Integral(
					  histogram->FindBin(min),
					  histogram->FindBin(max));
	double full_int = histogram->Integral();
	double N = entries * part_int / full_int;

	return N;
}


#endif // FIERZ_FITTER
