#ifndef FIERZ_FITTER
#define FIERZ_FITTER

/// UCNA includes
#include "G4toPMT.hh"
#include "PenelopeToPMT.hh"
#include "CalDBSQL.hh"

/// ROOT includes
#include <TH1.h>
#include <TLegend.h>
#include <TF1.h>
#include <TVirtualFitter.h>
#include <TList.h>
#include <TStyle.h>
#include <TApplication.h>
#include <TMatrixD.h>
#include <TNtuple.h>
#include <TLeaf.h>
#include <TString.h>
#include <TRandom2.h>
#include <TRegexp.h>

/// c++ includes
#include <iostream>
#include <fstream>
#include <string>
#include <limits>
#include <iomanip>

/// c includes
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

using TMath::Sqrt;
using TMath::IsNaN;


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
const double    M_F         = 1;                    /// Fermi matrix element
const double    M_GT        = 3*lambda*lambda;      /// Gamow-Teller matrix element
const double    M_n         = M_F + M_GT;           /// neutron matrix element (was M)
const double    I_0         = 1.63632;              /// 0th moment          (25)
const double    I_1         = 1.07017;              /// 1st moment          (15)
const double    x_1         = I_1/I_0;              /// first m/E moment    (9)


using namespace std;



struct UCNAhistogram : TH1D {
    //int side;
    //int spin;
    //vector<double> energy;        
    //vector<double> values;        
    //vector<double> errors;

    /*
    UCNAhistogram(int bins, double min, double max) 
      : TH1D(name, title, bins, min, max),
        name(""),
        title(""),
        bins(bins), min(min), max(max)
        //histogram(NULL),
        //energy(bins),        
        //values(bins),        
        //errors(bins)
    {}
    */
       
    UCNAhistogram(TString name, TString title, int bins, double min, double max) 
      : TH1D(name, title, bins, min, max)
        //name(name),
        //title(title),
        //bins(bins), min(min), max(max)
        //histogram(NULL),
        //energy(bins),        
        //values(bins),        
        //errors(bins),
    {
        assert(min <= max);
        //histogram = new TH1D(name, title, bins, min, max);
    }

    int fill(TString filename, TString name, TString title);
    int fill(TString filename);
    void save(TString filename, TString name, TString title);
    //void save(TString filename) {save(filename,1,1);}
    void draw(TString name, TString title,
              TCanvas* canvas, TLegend* legend, 
              TString draw, int color, int marker);
    void save(TString filename, double ax = 1, double ay = 1);
    void snapshot(int every = 10);

    bool test_min();
    bool test_min(double min);
    bool test_max();
    bool test_max(double max);
    bool test_range();
    bool test_range(double min, double max);
    bool test_match(const UCNAhistogram & other);
    bool test_compatable(UCNAhistogram & other);

    double normalize(double min, double max);
    double normalize();
    double GetEffectiveEntries(double min, double max);

    double chi2(const UCNAhistogram & other);
};



struct UCNAmodel {
    TString name;
    TString title;
    int     bins;
    double  min, max;
    double  Neff;
    double spin_ratio;

    TRandom2 rand;

    TNtuple* ntuple;     /// another way to store the raw data
    UCNAhistogram* counts[2][2]; // TODO make member not pointer
    UCNAhistogram super_ratio;
    UCNAhistogram super_sum;
    UCNAhistogram asymmetry;
    // Add Yup and Ydown ?


    //UCNAmodel():UCNAmodel("","",0,0,0) {exit(1);} /// default constructor needed by stl map.

    UCNAmodel(TString name, TString title, int bins, double min, double max) 
      : name(name), title(title),
        bins(bins), min(min), max(max),
        rand(0),
        super_ratio(name+"_super_ratio",title+" Super Ratio",bins,min,max),
        super_sum(name+"_super_sum",title+" Super Sum",bins,min,max),
        asymmetry(name+"_asymmetry",title+" Asymmetry",bins,min,max)
        /*
        KEbins(0), KEmin(min), KEmax(max),
        KEbins_A(0), KEmin_A(min), KEmax_A(max),
        KEbins_b(0), KEmin_b(min), KEmax_b(max),
        */
        //fedutial_cut(0), fidcut2(0)//, bin_resolution(0)
    {
        assert(min <= max);
        for (int side = 0; side < 2; side++) {
            TString sub_name = name;
            TString sub_title = title;
            if (not side) {
                sub_name += "_E";
                sub_title += " East";
            } else {
                sub_name += "_W";
                sub_title += " West";
            }
            for (int spin = 0; spin < 2; spin++) {
                if (not spin) {
                    sub_name += "_off";
                    sub_title += " AFP Off";
                } else {
                    sub_name += "_on";
                    sub_title += " AFP On";
                }
                counts[side][spin] = new UCNAhistogram(sub_name, sub_title, bins, min, max);
            }
        }
        ntuple = new TNtuple("mc_ntuple", "MC NTuple", "s:load:energy");
        /*
        KEbins = (KEmax - KEmin)/bin_resolution;
        KEbins_A = (KEmax_A - KEmin_A)/bin_resolution;
        KEbins_b = (KEmax_b - KEmin_b)/bin_resolution;
        */
    }

    UCNAmodel & operator=(const UCNAmodel & other);

    void SetAllBranches(TChain *chain);
    double asymmetry_chi2(double A, double b);
    //int fill(TString filename, TString name, TString title);
    int fill(TString pattern, int first, int last, 
             TString name, TString title, 
             int type, double flip);
             //int type = 0, double flip = -1);
    int fill(TString filename, TString name, TString title, 
             int type, double flip);
             //int type = 0, double flip = -1);
    int fill(TChain *chain, int type, double flip);
             //int type = 0, double flip = -1);
    void save(TString filename, TString name, TString title);
    void save(TString filename);

    /// accessing data
    bool test_counts();
    bool test_construction();
    void get_counts(int bin, double n[2][2]);
    void get_counts(int bin, double n[2][2], double e[2][2]);

    /// compute super sum
    double compute_super_sum(double n[2][2]);
    double compute_super_sum(double n[2][2], double e[2][2], 
                             double& S, double& error);
    double compute_super_sum(int bin, double& count, double& error);
    double compute_super_sum(int bin);
    TH1D& compute_super_sum();
    TH1D& compute_super_sum(double min, double max, 
                            int& min_bin, int& max_bin);
    TH1D& compute_super_sum(int min_bin, int max_bin);
    #if 0
    /// compute super ratio
    double compute_super_ratio(double n[2][2]);
    double compute_super_ratio(double n[2][2], double e[2][2], 
                               double& S, double& error);
    double compute_super_ratio(int bin, double& count, double& error);
    double compute_super_ratio(int bin);
    TH1D& compute_super_ratio();
    TH1D& compute_super_ratio(double min, double max, 
                             int& min_bin, int& max_bin);
    TH1D& compute_super_ratio(int min_bin, int max_bin);
    #endif

    /// compute asymmetry (and super ratio)
    double compute_asymmetry(double n[2][2]);
    double compute_asymmetry(double n[2][2], double e[2][2], 
                             double& S, double& error);
    double compute_asymmetry(int bin, double& count, double& error);
    double compute_asymmetry(int bin);
    TH1D& compute_asymmetry();
    TH1D& compute_asymmetry(double min, double max, 
                            int& min_bin, int& max_bin);
    TH1D& compute_asymmetry(int min_bin, int max_bin);
};


struct UCNAEvent {
    double EdepQ;
    double Edep;
    double MWPCEnergy;
    double ScintPos;
    double MWPCPos;
    double time;
    double primMomentum;
    double spin;
    int type;
};

struct UCNAFierzFitter {
    TString name;
    TString title;
    int bins;                       /// number of bins to use fit spectral plots
    double min;                     /// min kinetic energy for plots
    double max;                     /// max kinetic range for plots

    int fit_bins;                   /// number of bins to use fit spectral plots
    double fit_min;                 /// min kinetic energy for asymmetry fit
    double fit_max;                 /// max kinetic range for asymmetry fit

    double Nsim_data;               /// = data.super_sum.GetEntries();
    double Nall_data;               /// = data.super_sum.GetEffectiveEntries(KEmin, KEmax);
    double Nfit_data;               /// = data.super_sum.GetEffectiveEntries(fit_min, fit_max);
    double Nfit_vector;             /// = vector.super_sum.GetEffectiveEntries(fit_min, fit_max);
    double Nfit_axial;              /// = axial.super_sum.GetEffectiveEntries(fit_min, fit_max);
    double Nfit_fierz;              /// = fierz.super_sum.GetEffectiveEntries(fit_min, fit_max);
    double Neff;                    /// = Nfit_data*Nfit_vector/(Nfit_data + Nfit_vector);

    UCNAmodel data;                 /// Measured foreground data to fit.
    UCNAmodel back;                 /// Measured background data to remove.

    UCNAmodel vector;               /// Standard Model vector Monte Carlo spectrum.
    UCNAmodel axial;                /// Standard Model axial-vector Monte Carlo spectrum.
    /// TODO UCNAmodel axial[2];    /// Standard Model axial-vector Monte Carlo spectrum.
    UCNAmodel fierz;                /// Fierz (Scaler + tensor) Monte Carlo spectrum.
    UCNAmodel fit;                  /// Vector + axial + Fierz Monte Carlo best fit.

    /// cuts and settings
    unsigned nToSim = 5e7;			/// how many triggering events to simulate
    double spin_ratio = 1/1.68; 	/// afp off probability per neutron (0.68/1.68 for on)
    double fedutial_radius = 50;    /// radial cut in millimeters 
    int type = 0;
    //double fidrad2;               /// mm^2 radial cut.

    /// set up free fit parameters with best guess
    //static const int nPar = 3;
    //TString paramNames[3] = {"A", "b", "N"};
    //double paramInits[3] = {-0.12, 0, 1e1};

    UCNAFierzFitter(TString name, TString title, 
                    int bins, double min, double max)
      : name(name), title(title),
        bins(bins), min(min), max(max),
        fit_bins(bins), fit_min(min), fit_max(max),
        data(name+"_data", title+" data", bins, min, max),
        back(name+"_back", title+" background", bins, min, max),
        //sm(name+"_sm", title+" Standard Model Monte Carlo", bins, min, max),
        vector(name+"_vector", title+" Standard Model vector current", bins, min, max),
        axial(name+"_axial", title+" Standard Model axial-vector current", bins, min, max),
        fierz(name+"_fierz", title+"BSM Fierz current", bins, min, max),
        fit(name+"_fit", title+" fit", bins, min, max) 
    { assert(min < max); }

    UCNAFierzFitter(TString name, TString title,
                    int bins, double min, double max, 
                    int fit_bins, double fit_min, double fit_max)
      : name(name), title(title),
        bins(bins), min(min), max(max),
        fit_bins(fit_bins), fit_min(fit_min), fit_max(fit_max),
        data(name+"_data", title+" data", bins, min, max),
        back(name+"_back", title+" background", bins, min, max),
        vector(name+"_vector", title+" Standard Model vector current", bins, min, max),
        axial(name+"_axial", title+" Standard Model axial-vector vector current", bins, min, max),
        fierz(name+"_fierz", title+" BSM Fierz current", bins, min, max),
        fit(name+"_fit", title+" fit", fit_bins, fit_min, fit_max),
        fedutial_radius(50)
    { assert(min < max); }
/*
    UCNAFierzFitter(double min, double max, double fit_min, double fit_max)
      : bins((max-min)/resolution), min(min), max(max),
        bins((fit_max-fit_min)/resolution), fit_min(fit_min), fit_max(fit_max),
        data(name+"_data_", "UCNA data", bins, min, max),
        back(name+"_back_", "UCNA background", bins, min, max),
        sm(name+"_sm_", "Standard Model Monte Carlo", bins, min, max),
        axial(name+"_axial_", "Axial-vector Monte Carlo", bins, min, max),
        fierz(name+"_fierz_", "Fierz Monte Carlo", bins, min, max),
        fit(name+"_fit_", "Standard Model + Fierz best fit", bins, min, max) 
    { assert(min < max); }
    */
    void fill(TString vector_pattern, 
              TString axial_pattern, 
              // TODO TString axial_up_pattern, 
              // TODO TString axial_down_pattern, 
              TString fierz_pattern,
              int min, int max, /// TODO read filename pattern
              TString name, /// not sure if this is needed (or wanted) if the FF was constructed correctly
              int type, double flip);
              //int type = 0, double flip = -1);

    void save(TString filename);

    //void combined_chi2(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *p, Int_t /*iflag */);
    double asymmetry_chi2(double A, double b);
    double supersum_chi2(double b, double N);
    double combined_chi2(double A, double b, double N);
    void compute_asymmetry_fit(double A, double b);
    void compute_supersum_fit(double b, double N);
    void compute_fit(double A, double b, double N);
    void compute_data(double A, double b, double N);
    void compute_fit(TF1* func);
    void set_spin_ratio(double flip) { spin_ratio = flip; }
    double get_spin_ratio() { return spin_ratio; }
    void set_fedutial_radius(double r) { fedutial_radius = r; }
    double get_fedutial_radius() { return fedutial_radius; }

    //TF1* combined_fit(TH1D* asymmetry, TH1D* super_sum, TMatrixD &cov, TF1 *func);
    TF1* combined_fit(TMatrixD &cov, TF1 *func,  
        void (*)(Int_t&, Double_t*, Double_t&, Double_t*, Int_t));
        /*
        fierz.super_sum.histogram = 0;
        sm.super_sum.histogram = 0;
        for (int side = 0; side < 2; side++)
            for (int spin = 0; spin < 2; spin++) {
                fierz.raw[side][spin] = 0;
                sm.raw[side][spin] = 0;
            }
        fierz.super_sum.histogram = new TH1D("fierz_histogram", "", bins, min, max);
        sm.super_sum.histogram = new TH1D("standard_model_histogram", "", bins, min, max);
        for (int side = 0; side < 2; side++)
            for (int spin = 0; spin < 2; spin++) {
                fierz.raw[side][spin] = new TH1D("fierz_super_sum", "", bins, min, max);
                sm.raw[side][spin] = new TH1D("standard_model_super_sum", "", bins, min, max);
            }
    }
            */

    /*
    double evaluate(double *x, double*p) {
        double rv = 0;
        rv += p[0] * sm_histogram->GetBinContent(sm_histogram->FindBin(x[0]));        
        rv += p[1] * fierz_histogram->GetBinContent(fierz_histogram->FindBin(x[0]));
        return rv;
    }

    void normalize(TH1D* hist) {
        hist->Scale(1/(hist->GetBinWidth(2)*hist->Integral()));
    }

    void normalize(TH1D* hist, double min, double max) {
		int _min = hist->FindBin(min);
		int _max = hist->FindBin(max);
        hist->Scale(1/(hist->GetBinWidth(2)*hist->Integral(_min, _max)));
    }
    */
    void display(TString &plots_dir);
    double comupte_sizes();         /// compute effective size
    double print_sizes();           /// Output the number of events for data, mcs and fits and with cut info.
};


/*
double random(double min, double max) 
{
    double p = rand();
    return min + (max-min)*p/RAND_MAX;
}
*/


/// beta spectrum with little b term
double fierz_beta_spectrum(const double *val, const double *par) ;
/*
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
*/


/// beta spectrum with expected x^-n and beta^m
double beta_spectrum(const double *val, const double *par);
/*
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
*/


double evaluate_expected_fierz(double m, double n, double min, double max);
/*
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
	double rv = h1->Integral(0, integral_size) / h2->Integral(0, integral_size);
    delete h1;
    delete h2;
    return rv;
}
*/


double evaluate_expected_fierz(double min, double max);
/*
{
	return evaluate_expected_fierz(0, 1, min, max, integral_size);
}
*/



/*
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
*/


double asymmetry_fit_func(double *x, double *par);
double fierzratio_fit_func(double *x, double *par);
double GetEntries(TH1D* histogram, double min, double max);



/**
 * compute_corrected_asymmetry
 * 
 * @param rate_histogram[2][2]
 * @param correction
TH1D* compute_corrected_asymmetry(TH1D* rate_histogram[2][2], TH1D* correction) 
{
    TH1D *asymmetry_histogram = new TH1D(*(rate_histogram[0][0]));
    int bins = asymmetry_histogram->GetNbinsX();
    for (int bin = 1; bin < bins; bin++) 
	{
        double r[2][2];
        for (int side = 0; side < 2; side++)
            for (int spin = 0; spin < 2; spin++)
                r[side][spin] = rate_histogram[side][spin]->GetBinContent(bin);
        double sqrt_super_ratio = Sqrt((r[0][0] * r[1][1]) / (r[0][1] * r[1][0]));
        if ( IsNaN(sqrt_super_ratio) ) 
            sqrt_super_ratio = 0;
		double denom = 1 + sqrt_super_ratio;
        double asymmetry = (1 - sqrt_super_ratio) / denom;
		double sqrt_inverse_sum = Sqrt(1/r[0][0] + 1/r[1][1] + 1/r[0][1] + 1/r[1][0]);
		double asymmetry_error = sqrt_inverse_sum * sqrt_super_ratio / (denom * denom);  
		double K = asymmetry_histogram->GetBinCenter(bin);
		double E = K + m_e;                   /// electron energy
		double p = Sqrt(E*E - m_e*m_e);       /// electron momentum
		double beta = p/E;				      /// v/c
        asymmetry_histogram->SetBinContent(bin, -2*asymmetry/beta);
        asymmetry_histogram->SetBinError(bin, 2*asymmetry_error/beta);
        asymmetry_histogram->Multiply(correction);
        //printf("Setting bin content for corrected asymmetry bin %d, to %f\n", bin, asymmetry);
    }
    return asymmetry_histogram;
}
*/

#if 0    
/// this causes a linking error but I don't know why. I would like to use these at some point so don't erase.
TH1D* compute_rate_function(TH1D* rate_histogram[2][2], 
                            double (*rate_function)(double r[2][2]))
{
    TH1D *out_histogram = new TH1D(*(rate_histogram[0][0]));
    int bins = out_histogram->GetNbinsX();
    for (int bin = 1; bin < bins+1; bin++) {
        double r[2][2];
        for (int side = 0; side < 2; side++)
            for (int spin = 0; spin < 2; spin++)
                r[side][spin] = rate_histogram[side][spin]->GetBinContent(bin);

        double value = rate_function(r);
        out_histogram->SetBinContent(bin, value);
    }
    return out_histogram;
}


TH1D* compute_rate_function(TH1D* rate_histogram[2][2], 
                            double (*rate_function)(double r[2][2]),
                            double (*error_function)(double r[2][2])) 
{
    TH1D *out_histogram = new TH1D(*(rate_histogram[0][0]));
    int bins = out_histogram->GetNbinsX();
    for (int bin = 1; bin < bins+2; bin++) {
        double r[2][2];
        double e[2][2];
        for (int side = 0; side < 2; side++)
            for (int spin = 0; spin < 2; spin++) {
                r[side][spin] = rate_histogram[side][spin]->GetBinContent(bin);
                e[side][spin] = rate_histogram[side][spin]->GetBinError(bin);
            }

        double value = 0;
        if (rate_function)
            value = rate_function(r);

        double error = 0; 
        if (error_function)
            error = error_function(e);

        out_histogram->SetBinContent(bin, value);
        out_histogram->SetBinError(bin, error);
    }
    return out_histogram;
}


/*
TH1D* compute_rate_error_function(TH1D* rate_histogram[2][2], 
                                  double (*rate_error_function)(double r[2][2])) 
{
    TH1D *out_histogram = new TH1D(*(rate_histogram[0][0]));
    int bins = out_histogram->GetNbinsX();
    for (int bin = 1; bin < bins+2; bin++) {
        double sr[2][2];
        for (int side = 0; side < 2; side++)
            for (int spin = 0; spin < 2; spin++)
                sr[side][spin] = rate_histogram[side][spin]->GetBinError(bin);
    }
    return out_histogram;
}
*/


double bonehead_sum(double r[2][2]) {
    return r[0][0]+r[0][1]+r[1][0]+r[1][1];
}

double bonehead_asymmetry(double r[2][2]) {
    return (r[0][0]-r[0][1])/(r[1][0]+r[1][1]);
}

double super_ratio_asymmetry(double r[2][2]) {
    double super_ratio = (r[0][0]*r[1][1])/(r[0][1]*r[1][0]);
    double sqrt_super_ratio = Sqrt(super_ratio);
    if ( IsNaN(sqrt_super_ratio) ) 
        sqrt_super_ratio = 0;
    return (1-sqrt_super_ratio)/(1+sqrt_super_ratio);
}


/*
double super_sum_error(double r[2][2]) {
    double super_ratio = (r[0][0] * r[1][1]) / (r[0][1] * r[1][0]);
    double sqrt_super_ratio = Sqrt(super_ratio);
    if ( IsNaN(sqrt_super_ratio) ) 
        sqrt_super_ratio = 0;
    //return (1 - sqrt_super_ratio) / (1 + sqrt_super_ratio);
}
*/
#endif

#endif // FIERZ_FITTER
