/// UCNA includes
#include "G4toPMT.hh"
#include "PenelopeToPMT.hh"
#include "CalDBSQL.hh"
#include "FierzFitter.hh"

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

/// C++ includes
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

/// C includes
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

using std::setw;

/// settings
double expected_fierz = 0.6540;				/// full range (will get overwritten) 
//static double expected_fierz = 0.6111;	/// for range 150-600
//static double expected_gluck = 11.8498;   /// for range 150-600
static unsigned nToSim = 5e7;				/// how many triggering events to simulate
///static double afp_on_prob = 0.68/1.68; 		/// afp on probability per neutron
static double afp_off_prob = 1/1.68; 	    /// afp off probability per neutron
int large_prime = 179430331;                /// a large prime number
static int bins = 150;						/// replace with value from data or smoothing fit
double scale_x = 1.00000;

double min_E = 220;                         /// min energy from the 2013 paper
double max_E = 670;                         /// max range from the 2013 paper
double fedutial_cut = 50;                   /// radial cut in millimeters 
double fidcut2 = 50*50;                     /// mm^2 radial cut

/// set up free fit parameters with best guess
static const int nPar = 3;
TString iniParamNames[3] = {"A", "b", "N"};
double iniParams[3] = {-0.12, 0, 1e6};

/// Path to experiment data files.
TString data_dir = "/media/hickerson/boson/Data/OctetAsym_Offic_2010_FINAL/"; 

/// Path to Monte Carlo files.
TString mc_dir = "/home/xuansun/Documents/SimData_Beta/";

//double expected[3][3];                      /// expected values (based on the energy range)
//TMatrixD expected(3,3);                      /// expected values (based on the energy range)
//                                            /// needs to be visible to the chi^2 code

// ug. Needs to be static
//FierzHistogram mc(0,1500,bins);

/**
 * x[0] : kenetic energy
 * p[0] : b, fierz term
double theoretical_fierz_spectrum(double *x, double*p) 
{
    double rv = 0;
    //unsigned n = ff.mc.sm_histogram->FindBin(p[3]*x[0]*x[0] + p[2]*x[0] - p[1]);        
    unsigned n = ff.sm.super_sum.histogram->FindBin(p[2]*x[0] - p[1]);        
    //unsigned n = mc.sm_histogram->FindBin(x[0]);        
    double b = p[0];
    double norm = 1 + expected_fierz * b;
    rv += ff.sm_super_sum_histogram->GetBinContent(n) / norm;
    rv += b * expected_fierz * ucna.fierz.super_sum.histogram->GetBinContent(n) / norm;
    return rv;
}
 */

/*
/// beta spectrum with little b term
double fierz_beta_spectrum(const double *val, const double *par) 
{
	const double K = val[0];                    /// kinetic energy
	if (K <= 0 or K >= Q)
		return 0;                               /// zero outside range

	const double b = par[0];                    /// Fierz parameter
	const int    n = par[1];                    /// Fierz exponent
	const double E = K + m_e;                   /// electron energy
	const double e = Q - K;                     /// neutrino energy
	const double p = sqrt(E*E - m_e*m_e);       /// electron momentum
	const double x = pow(m_e/E,n);              /// Fierz term
	const double f = (1 + b*x)/(1 + b*x_1);     /// Fierz factor
	const double k = 1.3723803e-11/Q;           /// normalization factor
	const double P = k*p*e*e*E*f*x;             /// the output PDF value

	return P;
}


/// beta spectrum with little b term
double beta_spectrum(const double *val, const double *par) 
{
	const double K = val[0];                    /// kinetic energy
	const int n = par[0];                    	/// Fierz exponent
	if (K <= 0 or K >= Q)
		return 0;                               /// zero outside range

	const double E = K + m_e;                   /// electron energy
	const double p = sqrt(E*E - m_e*m_e);       /// electron momentum 
	const double e = (Q - K) / Q;               /// reduced neutrino energy
	const double x = pow(m_e/E,n);              /// Fierz term
	const double k = 1.3723803E-11/Q;           /// normalization factor
	const double P = k*p*e*e*E*x;               /// the output PDF value

	return P;
}
*/


unsigned deg = 4;
double mc_model(double *x, double*p) 
{
    double _exp = 0;
    double _x = x[0] / m_e;
    for (int i = deg; i >= 0; i--)
        _exp = p[i] + _x * _exp;
    return TMath::Exp(_exp);
}


void normalize(TH1D* hist) 
{
    hist->Scale(1/(hist->GetBinWidth(2)*hist->Integral()));
}


void normalize(TH1D* hist, double min, double max) 
{
	int _min = hist->FindBin(min);
	int _max = hist->FindBin(max);
	hist->Scale(1/(hist->GetBinWidth(2)*hist->Integral(_min, _max)));
}


/*
double evaluate_expected_fierz(double min, double max) 
{
    TH1D *h1 = new TH1D("beta_spectrum_fierz", "Beta spectrum with Fierz term", integral_size, min, max);
    TH1D *h2 = new TH1D("beta_spectrum", "Beta Spectrum", integral_size, min, max);
	for (int i = 0; i < integral_size; i++)
	{
		double K = min + double(i)*(max-min)/integral_size;
		double par1[2] = {0, 1};
		double par2[2] = {0, 0};
		double y1 = fierz_beta_spectrum(&K, par1);
		double y2 = fierz_beta_spectrum(&K, par2);
		h1->SetBinContent(K, y1);
		h2->SetBinContent(K, y2);
	}
	return h1->Integral(0, integral_size) / h2->Integral(0, integral_size);
}
*/


/**
 * Si := (r[0][0] r[1][1]) / (r[0][1] * r[1][0])
 */
TH1D* compute_super_ratio(TH1D* rate_histogram[2][2], TH1D* super_ratio_histogram = NULL) 
{
    for (int side=0; side<2; side++)
        for (int spin=0; spin<2; spin++)
            if (not rate_histogram[side][spin]) {
                cout<<"Error: rate histogram on the "
                    <<(side? "west":"east")<<" side with afp "
                    <<(spin? "on":"off")<<" is not constructed.\n";
                exit(1);
            }

    if (not super_ratio_histogram)
        super_ratio_histogram = new TH1D(*(rate_histogram[0][0]));

    int bins = super_ratio_histogram->GetNbinsX();
	cout<<"Number of bins "<<bins<<endl;
    for (int bin = 1; bin < bins; bin++) {
        double r[2][2];
        for (int side = 0; side < 2; side++)
            for (int spin = 0; spin < 2; spin++)
                r[side][spin] = rate_histogram[side][spin]->GetBinContent(bin);
        double super_ratio = r[0][0]*r[1][1]/r[0][1]/r[1][0];
        if (TMath::IsNaN(super_ratio)) {
            cout<<"Warning: super ratio in bin "<<bin<<" is not a number:\n"
                <<"Was "<<super_ratio<<". Setting to zero and continuing.\n";
            super_ratio = 0;
        }
        super_ratio_histogram->SetBinContent(bin, super_ratio);
        super_ratio_histogram->SetBinError(bin, 0.01);   // TODO compute correctly!!
        cout<<"Warning: super ratio is not computed correctly.\n";
    }
    return super_ratio_histogram;
}


/**
 * S := sqrt(r[0][0] r[1][1]) + sqrt(r[0][1] + r[1][0])
 */
TH1D* compute_super_sum(TH1D* rate_histogram[2][2], TH1D* super_sum_histogram = NULL) 
{
    for (int side=0; side<2; side++)
        for (int spin=0; spin<2; spin++)
            if (not rate_histogram[side][spin]) {
                cout<<"Error: rate histogram on side: "
                         <<(side? "west":"east")<<"and afp: "
                         <<(spin? "on":"off")<<"is not constructed.\n";
                exit(1);
            }

    if (not super_sum_histogram) {
        cout<<"Warning: super sum histogram is not constructed.\n";
        super_sum_histogram = new TH1D(*(rate_histogram[0][0]));
    }

    int bins = super_sum_histogram->GetNbinsX();
    for (int side=0; side < 2; side++)
        for (int spin = 0; spin < 2; spin++) {
            int ss_bins = rate_histogram[side][spin]->GetNbinsX();
            if (bins != ss_bins) {
                cout<<"Error: super sum and side spin histogram sizes don't match.\n";
                cout<<"super sum bins: "<<bins<<"\n";
                exit(1);
            }
            if (ss_bins <= 0) {
                cout<<"Error: bad bin number.\n";
                cout<<"super sum bins: "<<bins<<"\n";
                exit(1);
            }
        }

    for (int bin=1; bin < bins+2; bin++) {
        double r[2][2];
        for (int side = 0; side < 2; side++)
            for (int spin = 0; spin < 2; spin++)
                r[side][spin] = rate_histogram[side][spin]->GetBinContent(bin);

        double super_sum = TMath::Sqrt(r[0][0] * r[1][1]) + TMath::Sqrt(r[0][1] * r[1][0]);
        double error = TMath::Sqrt(1/r[0][0] + 1/r[1][0] + 1/r[1][1] + 1/r[0][1]) * super_sum;
        if (TMath::IsNaN(super_sum)) {
            super_sum = 0;
            cout<<"Warning: super sum is not a number: "<<super_sum<<".\n";
        } else if (super_sum == 0)
            cout<<"Warning: super sum is zero.\n";

        if (TMath::IsNaN(error)) {
			error = 0;
            cout<<"Warning: super sum error: division by zero.\n";
        } else if (error <= 0) 
            cout<<"Warning: super sum error: error is non positive.\n";

        if (bin % 10 == 0)
            printf("Setting bin content for super sum bin %d, to %f\n", bin, super_sum);

        super_sum_histogram->SetBinContent(bin, super_sum);
        super_sum_histogram->SetBinError(bin, error);
    }
    return super_sum_histogram;
}


TH1D* compute_asymmetry(TH1D* rate_histogram[2][2]) {
    TH1D *asymmetry_histogram = new TH1D(*(rate_histogram[0][0]));
    int bins = asymmetry_histogram->GetNbinsX();
    for (int bin=1; bin <= bins; bin++) 
	{
        double r[2][2];
        for (int side = 0; side < 2; side++)
            for (int spin = 0; spin < 2; spin++)
                r[side][spin] = rate_histogram[side][spin]->GetBinContent(bin);
        double super_ratio = TMath::Sqrt(r[0][0]*r[1][1]/r[0][1]/r[1][0]);
        if (TMath::IsNaN(super_ratio)) {
            cout<<"Warning: super ratio in bin "<<bin<<" is not a number:\n";
            cout<<"Was "<<super_ratio<<". Setting to zero and continuing.\n";
            super_ratio = 0;
        }

		double norm = 1 + super_ratio;
        double asymmetry = (1 - super_ratio) / norm;
        asymmetry_histogram->SetBinContent(bin, asymmetry);

		double inv_sum = TMath::Sqrt(1/r[0][0] + 1/r[1][1] + 1/r[0][1] + 1/r[1][0]) / norm;
		double asymmetry_error = super_ratio * inv_sum / norm;  
        if (TMath::IsNaN(asymmetry_error) or asymmetry_error <= 0) {
            cout<<"Warning: super ratio in bin "<<bin<<" is not a number:\n";
            cout<<"Was "<<asymmetry_error<<". Setting to 0.01 and continuing.\n";
            asymmetry_error = 0.01;
        }
        asymmetry_histogram->SetBinError(bin, asymmetry_error);
        //printf("Setting bin content for asymmetry bin %d, to %f\n", bin, asymmetry);
    }
    return asymmetry_histogram;
}


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
        double sqrt_super_ratio = TMath::Sqrt((r[0][0] * r[1][1]) / (r[0][1] * r[1][0]));
        if ( TMath::IsNaN(sqrt_super_ratio) ) 
            sqrt_super_ratio = 0;
		double denom = 1 + sqrt_super_ratio;
        double asymmetry = (1 - sqrt_super_ratio) / denom;
		double sqrt_inverse_sum = TMath::Sqrt(1/r[0][0] + 1/r[1][1] + 1/r[0][1] + 1/r[1][0]);
		double asymmetry_error = sqrt_inverse_sum * sqrt_super_ratio / (denom * denom);  
		double K = asymmetry_histogram->GetBinCenter(bin);
		double E = K + m_e;                   /// electron energy
		double p = sqrt(E*E - m_e*m_e);       /// electron momentum
		double beta = p / E;				  /// v/c
        asymmetry_histogram->SetBinContent(bin, -2*asymmetry/beta);
        asymmetry_histogram->SetBinError(bin, 2*asymmetry_error/beta);
        asymmetry_histogram->Multiply(correction);
        //printf("Setting bin content for corrected asymmetry bin %d, to %f\n", bin, asymmetry);
    }
    return asymmetry_histogram;
}


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
    double sqrt_super_ratio = TMath::Sqrt(super_ratio);
    if ( TMath::IsNaN(sqrt_super_ratio) ) 
        sqrt_super_ratio = 0;
    return (1-sqrt_super_ratio)/(1+sqrt_super_ratio);
}


/*
double super_sum_error(double r[2][2]) {
    double super_ratio = (r[0][0] * r[1][1]) / (r[0][1] * r[1][0]);
    double sqrt_super_ratio = TMath::Sqrt(super_ratio);
    if ( TMath::IsNaN(sqrt_super_ratio) ) 
        sqrt_super_ratio = 0;
    //return (1 - sqrt_super_ratio) / (1 + sqrt_super_ratio);
}
*/


using namespace std;
void output_histogram(TString filename, TH1D* h, double ax, double ay)
{
	using namespace std;
	ofstream ofs;
	//ofs.open(filename.c_str());
	ofs.open(filename);
	for (int i = 1; i < h->GetNbinsX(); i++)
	{
		double x = ax * h->GetBinCenter(i);
		//double sx = h->GetBinWidth(i);
		double r = ay * h->GetBinContent(i);
		double sr = ay * h->GetBinError(i);
		ofs<<x<<'\t'<<r<<'\t'<<sr<<endl;
	}
	ofs.close();
}


double random(double min, double max) 
{
    //double p = (large_prime*nSimmed % rand_bins)/rand_bins;
    double p = rand();
    return min + (max-min)*p/RAND_MAX;
}


#if 0
using namespace std;
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


// data need to be globals to be visible by func 
/*
vector<double> asymmetry_energy;        
vector<double> asymmetry_values;        
vector<double> asymmetry_errors;        

vector<double> super_ratio_energy;        
vector<double> super_ratio_values;        
vector<double> super_ratio_errors;        

vector<double> super_sum_energy;        
vector<double> super_sum_values;        
vector<double> super_sum_errors;        
*/

/*
UCNAModel ucna_data; // Need construction.
UCNAModel ucna_sm_mc; // Need construction.
UCNAModel ucna_fierz_mc; // Need construction.
//FierzHistogram mc(0,1500,bins);
*/
/// This needs to be static
UCNAFierzFitter ucna(bins, min_E, max_E);

void combined_chi2(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *p, Int_t /*iflag */  )
{
	double chi2=0, chi;
    int n = ucna.data.asymmetry.bins;
    double A=p[0], b=p[1], N=p[2]; // TODO make nPar correct here
	//double par[3] = {p[0],p[1],p[2]}; // A, b, N
	for (int i = 0; i < n; i++)
	{
		double E = ucna.data.asymmetry.histogram->GetBinCenter(i);
		double Y = ucna.data.asymmetry.histogram->GetBinContent(i);
        //double f = asymmetry_fit_func(&E,p);
        double f = A/(1 + b*m_e/(E+m_e));
		//double eY = ucna.data.asymmetry.histogram->GetBinError(i);
		double eY = 0.1;
        if (eY > 0) {
            chi = (Y-f)/eY;
            chi2 += chi*chi; 
        }
	}

    n = ucna.data.super_sum.bins;
	//double par[2] = {p[1], expected[0][1]};
	for (int i=0; i<n; i++) { 
		//E = ucna.data.super_sum.energy[i];
		/*double Y =      ucna.data .super_sum.values[i];
        double f = p[1]*ucna.sm   .super_sum.values[i] 
                 + p[2]*ucna.fierz.super_sum.values[i];
        double eY =     ucna.data .super_sum.errors[i];*/
		//double E      = ucna.data .super_sum.histogram->GetBinCenter(i);
		//chi = (fierzratio_values[i] - fierzratio_fit_func(&E,par)) / fierzratio_errors[i];
		double Y  = ucna.data .super_sum.histogram->GetBinContent(i);
        double eY = ucna.data .super_sum.histogram->GetBinError(i);
        double f  = N*ucna.sm .super_sum.histogram->GetBinContent(i) 
                  + N*b*ucna.fierz.super_sum.histogram->GetBinContent(i);
        if (eY > 0) {
		    chi = (Y-f)/eY;
		    chi2 += chi*chi; 
        }
	}
	fval = chi2; 
}



#if 1
//TF1* combined_fit(TH1D* asymmetry, TH1D* super_sum, double cov[nPar][nPar]) 
TF1* combined_fit(TH1D* asymmetry, TH1D* super_sum, TMatrixD &cov, TF1* func)
{ 
	/// create fit function
    if (not func) {
        func = new TF1("func", asymmetry_fit_func, ucna.min, ucna.max, nPar);
        func->SetParameters(iniParams);
        for (int i=0; i<nPar; i++)
            func->SetParName(i,iniParamNames[i]);
    }

	/// fill data structure for fit (coordinates + values + errors) 
	cout<<"Do global fit"<<endl;

	/// set up the minuit fitter
	TVirtualFitter::SetDefaultFitter("Minuit");
	TVirtualFitter *minuit = TVirtualFitter::Fitter(0,nPar);
	for (int i=0; i<nPar; ++i)
		minuit->SetParameter(i, func->GetParName(i), func->GetParameter(i), 1, 0, 0);
	minuit->SetFCN(combined_chi2);
	minuit->SetErrorDef(1);	        /// 1 for chi^2

	/// set print level
	double arglist[100];
	arglist[0] = 0;
	minuit->ExecuteCommand("SET PRINT", arglist, 1);

	/// minimize
	arglist[0] = 50;    /// number of function calls
	arglist[1] = 0.1;   /// tolerance
	minuit->ExecuteCommand("MIGRAD", arglist, nPar);

	double minParams[nPar];
	double parErrors[nPar];
	for (int i=0; i<nPar; ++i) {  
		minParams[i] = minuit->GetParameter(i);
		parErrors[i] = minuit->GetParError(i);
	}

	/// extract results from minuit
	double chi2, edm, errdef; 
	int nvpar, nparx;
	minuit->GetStats(chi2,edm,errdef,nvpar,nparx);

	func->SetParameters(minParams);
	func->SetParErrors(parErrors);
	func->SetChisquare(chi2);

	int ndf = asymmetry->GetNbinsX() + super_sum->GetNbinsX() - nvpar;
	func->SetNDF(ndf);
    
	TMatrixD matrix( nPar, nPar, minuit->GetCovarianceMatrix() );
	for (int i=0; i<nPar; i++)
		for (int j=0; j<nPar; j++)
			cov[i][j] = minuit->GetCovarianceMatrixElement(i,j);
	
	cout<<"\n    chi^2 = "<<chi2<<", ndf = "<<ndf<<", chi^2/ndf = "<<chi2/ndf<<".\n\n";

	return func; 
}
#endif

int fill_data(TString filename, TString title, 
              TString name, TH1D* histogram)
{
	/// load the files that contain data histograms
	TFile* tfile = new TFile(data_dir + filename);
	if (tfile->IsZombie()) {
		cout<<"Error loading "<<title<<":\n";
		cout<<"File not found: "<<filename<<".\n";
		exit(1);
	}

    if (histogram) {
		cout<<"Warning: histogram "<<title<<" already exists and is being deleted.\n";
        delete histogram;
    }

    histogram = (TH1D*)tfile->Get(name);
    if (not histogram) {
		cout<<"Error in file "<<filename<<":\n";
		cout<<"Error getting "<<title<<":\n";
		cout<<"Cannot find histogram named "<<name<<".\n";
        exit(1);
    }

	int entries=histogram->GetEntries();
	cout<<"Number of entries in "<<title<<" is "<<entries<<".\n";
    return entries;
}

int fill_simulation(TString filename, TString title, TString name, 
                    TH1D* histogram[2][2], TH1D* super_sum)
{
    for (int side=0; side<2; side++)
        for (int spin=0; spin<2; spin++)
            if (not histogram[side][spin]) {
                cout<<"Error: histogram for "<< name<<" is not constructed.\n";
                cout<<"Side: "<< side<<" Spin: "<<spin <<".\n";
                exit(1);
            }

    if (not super_sum) {
        cout<<"Error: super sum histogram for "<< name<<" is not constructed.\n";
        exit(1);
    }

	TFile* tfile = new TFile(mc_dir + filename);
	if (tfile->IsZombie()) {
		cout<<"Error loading "<<title<<":\n";
		cout<<"File not found: "<<filename<<".\n";
		exit(1);
	}

    TChain *chain = (TChain*)tfile->Get(name);
    if (not chain) {
		cout<<"Error in file "<<filename<<":\n";
		cout<<"Error getting "<<title<<":\n";
		cout<<"Cannot find chain or tree named "<<name<<".\n";
        exit(1);
    }

	int nEvents = chain->GetEntries();
	chain->SetBranchStatus("*",false);
	chain->SetBranchStatus("PID",true);
	chain->SetBranchStatus("side",true);
	chain->SetBranchStatus("type",true);
	chain->SetBranchStatus("Erecon",true);
	chain->SetBranchStatus("primMomentum",true);
    chain->SetBranchStatus("ScintPosAdjusted",true);
    chain->SetBranchStatus("ScintPosAdjusted",true);

	/* TFile* mc_tfile = new TFile("Fierz/mc.root", "recreate");
	if (mc_tfile->IsZombie())
	{
		cout << "Can't make MC file.\n";
		exit(1);
	}
	TNtuple* tntuple = new TNtuple("mc_ntuple", "MC NTuple", "s:load:energy"); */

	unsigned int nSimmed = 0;	/// events have been simulated after cuts
    int PID, side, type;
    double energy;
    double mwpcPosW[3], mwpcPosE[3], primMomentum[3];

    chain->SetBranchAddress("PID",&PID);
    chain->SetBranchAddress("side",&side);
    chain->SetBranchAddress("type",&type);
    chain->SetBranchAddress("Erecon",&energy);
	chain->SetBranchAddress("primMomentum",primMomentum);
    chain->GetBranch("ScintPosAdjusted")->GetLeaf("ScintPosAdjE")->SetAddress(mwpcPosE);
    chain->GetBranch("ScintPosAdjusted")->GetLeaf("ScintPosAdjW")->SetAddress(mwpcPosW);

    for (int evt=0; evt<nEvents; evt++) {
        chain->GetEvent(evt);

        /// cut out bad events
        if (PID!=1) 
            continue;

        /// radial fiducial cut
        double radiusW = mwpcPosW[0]*mwpcPosW[0] + mwpcPosW[1]*mwpcPosW[1]; 
        double radiusE = mwpcPosE[0]*mwpcPosE[0] + mwpcPosE[1]*mwpcPosE[1]; 
        if (radiusW > fidcut2 or radiusE > fidcut2) 
            continue;

        /// Type 0, Type I, Type II/III events 
        if (type<4) { 
            /// fill with loading efficiency 
            double p = random(0,1);
			double afp = (p < afp_off_prob)? -1 : +1;
            bool spin = (afp < 0)? EAST : WEST;
            //cout <<"energy: "<<energy<<" side: "<<"side: "<<side
            //          <<" spin: "<<spin<<" afp: "<<afp<<" p: "<<p<<".\n";
            histogram[side][spin]->Fill(energy, 1);
			//tntuple->Fill(side, spin, energy);
			nSimmed++;
        }
        /*  hEreconALL->Fill(Erecon);
        if (type==0) 
            hErecon0->Fill(Erecon);
        else if (type==1) 
            hErecon1->Fill(Erecon);
        else if (type==2 or type==3) 
            hErecon23->Fill(Erecon); */

		/// break when enough data has been generated.
		if(nSimmed >= nToSim)
			break;
    }    
     
    #if 0
	while (false) {
		/// perform energy calibrations/simulations to fill class 
        /// variables with correct values for this simulated event
		//G2P.recalibrateEnergy();
		
		/// check the event characteristics on each side
		for(Side s=EAST; s<=WEST; ++s) {
			/// get event classification type. 
            /// TYPE_IV_EVENT means the event didn't trigger this side.
			/// TODO EventType tp = G2P.fType;

			/// skip non-triggering events, or those outside 50mm position 
            /// cut (you could add your own custom cuts here, if you cared)
            /// TODO place all cuts here!
			//if (tp >= TYPE_I_EVENT or !G2P.passesPositionCut(s) or G2P.fSide != s)
			//	continue;
			
			/// print out event info, (simulated) reconstructed true energy
            /// and position, comparable to values in data
			#ifdef DEBUG
			printf("Event on side %c: type=%i, Erecon=%g @ position (%g,%g), %d\n",
				   sideNames(s), tp, G2P.getErecon(), G2P.wires[s][X_DIRECTION].center, 
				   G2P.wires[s][Y_DIRECTION].center, (unsigned)G2P.getAFP());

			/// print out event primary info, only available in simulation
			printf("\tprimary KE=%g, cos(theta)=%g\n", G2P.ePrim, G2P.costheta);
			#endif 

			/// fill with loading efficiency 
			/// NOT the way to do this anymore 
            //bool load = (nSimmed % 100 < loading_prob);

			/// calculate the energy with a distortion factor
			/// double energy = scale_x * G2P.getErecon();
            double energy = 0;  /// TODO XXXXXXXX
            histogram[s][load]->Fill(energy, 1);
			tntuple->Fill(s, load, energy);
			nSimmed++;
		}
		
		/// break when enough data has been generated.
		if(nSimmed >= nToSim)
			break;
	}
    #endif
    
	cout<<"Total number of Monte Carlo entries without cuts: "<<nEvents<<endl;
	cout<<"Total number of Monte Carlo entries with cuts: "<<nSimmed<<endl;

	//tntuple->SetDirectory(mc_tfile);
	//tntuple->Write();

	/// compute and normalize super sum
    compute_super_sum(histogram, super_sum);
    //normalize(super_sum, min_E, max_E);

    /*
    for (int side=0; side<2; side++)
        for (int spin=0; spin<2; spin++)
            normalize(histogram[side][spin], min_E, max_E);
            */

    return nSimmed;
}


/// compute little b factor using the Fierz ratio method
TH1D* compute_fierz_ratio(TH1D* data_histogram, TH1D* sm_histogram) {
    TH1D *fierz_ratio_histogram = new TH1D(*data_histogram);
	fierz_ratio_histogram->SetName("fierz_ratio_histogram");
    //fierz_ratio_histogram->Divide(ucna.data.super_sum.histogram, ucna.sm.super_sum.histogram);
    int bins = data_histogram->GetNbinsX();
    for (int bin = 1; bin < bins+1; bin++) {
		double X = data_histogram->GetBinContent(bin);
		double Y = sm_histogram->GetBinContent(bin);
		double Z = Y > 0 ? X/Y : 0;

		fierz_ratio_histogram->SetBinContent(bin, Z);

		double x = data_histogram->GetBinError(bin);
		double y = sm_histogram->GetBinError(bin);
		fierz_ratio_histogram->SetBinError(bin, Z*TMath::Sqrt(x*x/X/X + y*y/Y/Y));
	}
    fierz_ratio_histogram->GetYaxis()->SetRangeUser(0.9,1.1); // Set the range
    fierz_ratio_histogram->SetTitle("Ratio of UCNA data to Monte Carlo");
	cout<<data_histogram->GetNbinsX()<<endl;
	cout<<sm_histogram->GetNbinsX()<<endl;

    /// fit the Fierz ratio 
	char fit_str[1024];
    sprintf(fit_str, "1+[0]*(%f/(%f+x)-%f)", m_e, m_e, expected_fierz);
    TF1 *fierz_fit = new TF1("fierz_fit", fit_str, min_E, max_E);
    fierz_fit->SetParameter(0,0);
	fierz_ratio_histogram->Fit(fierz_fit, "Sr");

	/// A fit histogram for output to gnuplot
    TH1D *fierz_fit_histogram = new TH1D(*ucna.data.super_sum.histogram);
	for (int i = 0; i < fierz_fit_histogram->GetNbinsX(); i++)
		fierz_fit_histogram->SetBinContent(i, fierz_fit->Eval(fierz_fit_histogram->GetBinCenter(i)));

	/// compute chi squared
    double chisq = fierz_fit->GetChisquare();
    double NDF = fierz_fit->GetNDF();
	char b_str[1024];
	sprintf(b_str, "b = %1.3f #pm %1.3f", fierz_fit->GetParameter(0), fierz_fit->GetParError(0));
	char chisq_str[1024];
    printf("Chi^2 / (NDF-1) = %f / %f = %f\n", chisq, NDF-1, chisq/(NDF-1));
	sprintf(chisq_str, "#frac{#chi^{2}}{n-1} = %f", chisq/(NDF-1));

	/// draw the ratio plot
	fierz_ratio_histogram->SetStats(0);
    fierz_ratio_histogram->Draw();

	/// draw a legend on the plot
    TLegend* ratio_legend = new TLegend(0.3,0.85,0.6,0.65);
    ratio_legend->AddEntry(fierz_ratio_histogram, "Data ratio to Monte Carlo (Type 0)", "l");
    ratio_legend->AddEntry(fierz_fit, "Fierz term fit", "l");
    ratio_legend->AddEntry((TObject*)0, "1+b(#frac{m_{e}}{E} - #LT #frac{m_{e}}{E} #GT)", "");
    ratio_legend->AddEntry((TObject*)0, b_str, "");
    ratio_legend->AddEntry((TObject*)0, chisq_str, "");
    ratio_legend->SetTextSize(0.03);
    ratio_legend->SetBorderSize(0);
    ratio_legend->Draw();

	/// output for root
    TString fierz_ratio_pdf_filename = "mc/fierz_ratio.pdf";
    TCanvas *canvas = new TCanvas("fierz_canvas", "Fierz component of energy spectrum");
    if (not canvas) {
        cout<<"Can't open new canvas.\n";
        exit(0);
    }
    canvas->SaveAs(fierz_ratio_pdf_filename);

	/// output for gnuplot
	output_histogram("mc/super-sum-data.dat", data_histogram, 1, 1000);
	output_histogram("mc/super-sum-mc.dat", sm_histogram, 1, 1000);
	output_histogram("mc/fierz-ratio.dat", fierz_ratio_histogram, 1, 1);
	output_histogram("mc/fierz-fit.dat", fierz_fit_histogram, 1, 1);

	TFile* ratio_tfile = new TFile("Fierz/ratio.root", "recreate");
	if (ratio_tfile->IsZombie())
    {
		cout<<"Can't recreate MC file.\n";
		exit(1);
	}
	fierz_ratio_histogram->SetDirectory(ratio_tfile);
	fierz_ratio_histogram->Write();
	ratio_tfile->Close();

    return fierz_ratio_histogram;
}
    

int main(int argc, char *argv[])
{
	TApplication app("Extract Combined A + b Fitter", &argc, argv);
	//TH1::AddDirectory(kFALSE);
	/// Geant4 MC data scanner object
	/// G4toPMT G2P;
	/// use data from these MC files (most recent unpolarized beta decay, 
    /// includes Fermi function spectrum correction)
	/// note wildcard * in filename; 
    /// MC output is split up over many files, but G2P will TChain them together
	//G2P.addFile("/data2/mmendenhall/G4Out/2010/20120823_neutronBetaUnpol/analyzed_*.root");
	
	/// PMT Calibrator loads run-specific energy calibrations info for selected run
	/// and uses default Calibrations DB connection to most up-to-date though possibly unstable "mpm_debug"
	
	/// If you really want this to be random, 
    /// you will need to seed rand() with something 
    /// other than default.
	srand( time(NULL) );

    /// load the files that contain data histograms
    fill_data("Range_0-1000/CorrectAsym/CorrectedAsym.root",
              "2010 final official asymmetry",
              "hAsym_Corrected_C",
              ucna.data.asymmetry.histogram);

    fill_data("OctetAsym_Offic.root",
              "2010 final official supersum",
              "Total_Events_SuperSum",
              ucna.data.super_sum.histogram);

    fill_simulation("SimAnalyzed_Beta.root",
                    "Monte Carlo Standard Model beta spectrum",
                    "SimAnalyzed",
                    ucna.sm.raw,
					ucna.sm.super_sum.histogram);

    fill_simulation("SimAnalyzed_Beta_fierz.root",
                    "Monte Carlo Fierz beta spectrum",
                    "SimAnalyzed",
                    ucna.fierz.raw,
					ucna.fierz.super_sum.histogram);

	//histogram->SetDirectory(mc_tfile);
	//histogram->Write();

	//ucna.sm.super_sum.histogram->SetDirectory(mc_tfile);
	//ucna.sm.super_sum.histogram->Write();

	//ucna.fierz.super_sum.histogram->SetDirectory(mc_tfile);
	//ucna.fierz.super_sum.histogram->Write();

	//mc_tfile->Close();

    /* for (int side=0; side<2; side++)
        for (int spin=0; spin<2; spin++)
		{
            normalize(mc.fierz_histogram[side][spin], min_E, max_E);
            normalize(mc.sm_histogram[side][spin], min_E, max_E);
        }*/

    TCanvas *canvas = new TCanvas("fierz_canvas", "Fierz component of energy spectrum");
    if (not canvas) {
        cout<<"Can't open new canvas.\n";
        exit(0);
    }

	ucna.fierz.super_sum.histogram->SetStats(0);
    ucna.fierz.super_sum.histogram->SetLineColor(3);
    ucna.fierz.super_sum.histogram->Draw("");
    ucna.sm.super_sum.histogram->SetLineColor(1);
    ucna.sm.super_sum.histogram->Draw("Same");

    // fit a smooth model to the mc
    /*
    TF1 *mc_fit = new TF1("fierz_mc_fit", mc_model, 0, 1000, deg+1);
    mc_fit->SetParameter(0,-0.5);
    mc_fit->SetParameter(1,1.0);
    mc_fit->SetParameter(2,1.0);
    mc_fit->SetParameter(3,1.0);
    mc_fit->SetParameter(4,1.0);
    //mc_fit->SetParameter(5,1.0);
    //mc_fit->SetParameter(6,1.0);
    mc.sm_histogram->Fit("fierz_mc_fit");

    TString pdf_filename = "/data/kevinh/mc/fierz_models_to_fit.pdf";
    canvas->SaveAs(pdf_filename);
    */

    /**
     *  If you want a fast way to get the combined data spectrums for comparison, 
     *  I already have it extracted as a ROOT histogram for my own data/MC comparisons.
     *  You can read the ROOT TFile at:
     *      /media/hickerson/boson/Data/OctetAsym_Offic_2010_FINAL/OctetAsym_Offic.root
     *  and get the TH1D histograms named:
     *  Combined_Events_<S><afp><fg/bg><type> where 
     *      <S>='E' or 'W' is the side, 
     *      <afp>='0' or '1' for AFP Off/On, 
     *      <fg/bg>='0' for background or '1' for foreground data
     *      <type>='0','1','2' for type 0, I, or II/III backscatter events.
     */
#if 0
    TFile *ucna_data_tfile = new TFile(
        //"/home/mmendenhall/mpmAnalyzer/PostPlots/OctetAsym_div0/Combined/Combined.root");
	    //"/home/mmendenhall/mpmAnalyzer/PostPlots/OctetAsym_Offic_10keV_bins/Combined.root");
        //"/home/mmendenhall/mpmAnalyzer/PostPlots/OctetAsym_10keV_Bins/Combined");
		//"/home/mmendenhall/mpmAnalyzer/PostPlots/OctetAsym_10keV_Bins/OctetAsym_10keV_Bins.root");
		//"/home/mmendenhall/UCNA/PostPlots/OctetAsym_Offic/OctetAsym_Offic.root");
		//"/home/mmendenhall/Plots/OctetAsym_Offic/OctetAsym_Offic.root");
        "/media/hickerson/boson/Data/OctetAsym_Offic_2010_FINAL/OctetAsym_Offic.root");
	if (ucna_data_tfile->IsZombie())
	{
		//printf("File "+beta_filename+"not found.\n");
		cout << "Data file not found." << endl;
		exit(1);
	}

	/*
	TFile* ucna_correction_tfile = new TFile("Fierz/tree.root");
	if (ucna_correction_tfile->IsZombie())
	{
		//printf("File "+beta_filename+"not found.\n");
		cout << "Correction file not found." << endl;
		exit(1);
	}
	*/

	#define EVENT_TYPE -1 
    TH1D *ucna.data_raw[2][2] = {
	#if EVENT_TYPE == 0
        {   (TH1D*)ucna_data_tfile->Get("hEnergy_Type_0_E_Off"),
            (TH1D*)ucna_data_tfile->Get("hEnergy_Type_0_E_On")
        },{ (TH1D*)ucna_data_tfile->Get("hEnergy_Type_0_W_Off"),
            (TH1D*)ucna_data_tfile->Get("hEnergy_Type_0_W_On") }};
	#endif
	#if EVENT_TYPE == -1
        {   //(TH1D*)ucna_data_tfile->Get("Combined_Events_E010"),
            //(TH1D*)ucna_data_tfile->Get("Combined_Events_E110")
            (TH1D*)ucna_data_tfile->Get("hTotalEvents_E_Off;1"),
            (TH1D*)ucna_data_tfile->Get("hTotalEvents_E_On;1"),
        },{ //(TH1D*)ucna_data_tfile->Get("Combined_Events_W010"),
            //(TH1D*)ucna_data_tfile->Get("Combined_Events_W110")
            (TH1D*)ucna_data_tfile->Get("hTotalEvents_W_Off;1"),
            (TH1D*)ucna_data_tfile->Get("hTotalEvents_W_On;1") }};
	#endif
#endif
    /*
    ucna.data.raw[0][0]=(TH1D*)ucna_data_tfile->Get("hTotalEvents_E_Off;1");
    ucna.data.raw[0][1]=(TH1D*)ucna_data_tfile->Get("hTotalEvents_E_On;1");
    ucna.data.raw[1][0]=(TH1D*)ucna_data_tfile->Get("hTotalEvents_W_Off;1");
    ucna.data.raw[1][1]=(TH1D*)ucna_data_tfile->Get("hTotalEvents_W_On;1");
    */

    /* TODO figure out where these went.
    fill_data("OctetAsym_Offic.root",
              "2010 final official east afp off spectrum",
              "hTotalEvents_E_off;1",
              ucna.data.raw[0][0]);
    fill_data("OctetAsym_Offic.root",
              "2010 final official east afp on spectrum",
              "hTotalEvents_E_on;1",
              ucna.data.raw[0][1]);
    fill_data("OctetAsym_Offic.root",
              "2010 final official west afp off spectrum",
              "hTotalEvents_W_off;1",
              ucna.data.raw[1][0]);
    fill_data( "OctetAsym_Offic.root",
              "2010 final official west afp on spectrum",
              "hTotalEvents_W_on;1",
              ucna.data.raw[1][1]);
    printf("Number of bins in data %d\n", ucna.data.raw[0][0]->GetNbinsX());
    for (int side=EAST; side<=WEST; side++)
        for (int afp=EAST; afp<=WEST; afp++) {
            TString title = "2010 final official "+side?"west":"east";
            fill_data("OctetAsym_Offic.root",
                      "2010 final official west afp on spectrum",
                      "hTotalEvents_W_on;1",
                      ucna.data.raw[side][afp]);
        }
    */

    /* Already background subtracted...
        TH1D *background_histogram = (TH1D*)ucna_data_tfile->Get("Combined_Events_E000");
        ucna_data.raw->Add(background_histogram,-1);
        // normalize after background subtraction
        background_histogram->Draw("");
    */

/*
	for (int side = 0; side < 2; side++)
		for (int spin = 0; spin < 2; spin++)
		{
			cout << "Number of entries in (" 
					  << side << ", " << spin << ") is "
					  << (int)ucna.data.raw[side][spin]->GetEntries() << endl;
			if (ucna.data.raw[side][spin] == NULL)
			{
				puts("histogram is null. Aborting...");
				exit(1);
			}
		}
        */

	/*
    TH1D *ucna.correction_histogram = (TH1D*)ucna.correction_histogram->Get("Delta_3_C");
	if (not ucna.correction_histogram)
	{
		puts("Correction histogram is null. Aborting...");
		exit(1);
	}
	*/
    //TH1D *ucna.correction_histogram = new TH1D(*ucna.data.raw[0][0]);
	/*
	while (tfile has more entries)
	{
        double bin = (ucna.correction_file->GetBinContent(bin);
        double correction = ucna.correction_histogram->GetBinContent(bin);
        printf("Setting bin content for correction bin %d, to %f\n", bin, correction);
        ucna.data.super_sum.histogram->SetBinContent(bin, correction);
        ucna.data.super_sum.histogram->SetBinError(bin, correction_error);
    }
	*/

    /*
    TF1 *fit = new TF1("fierz_fit", theoretical_fierz_spectrum, 0, 1000, 3);
    fit->SetParameter(0,0.0);
    fit->SetParameter(1,0.0);
    fit->SetParameter(2,1.0);
    ucna.data.raw[0][0]->Fit("fierz_fit");
    double chisq = fit->GetChisquare();
    double N = fit->GetNDF();
    printf("Chi^2 / ( N - 1) = %f / %f = %f\n",chisq, N-1, chisq/(N-1));
    */

    TString fit_pdf_filename = "mc/fierz_fit_data.pdf";
    canvas->SaveAs(fit_pdf_filename);

    // compute and plot the super ratio
    /*
    TH1D *ucna.data.super_ratio.histogram = compute_super_ratio(ucna.data.raw);
    ucna.data.super_ratio.histogram->Draw();
    TString super_ratio_pdf_filename = "mc/super_ratio_data.pdf";
    canvas->SaveAs(super_ratio_pdf_filename);
    */

    // compute and plot the super ratio asymmetry 
    //TH1D *asymmetry_histogram = compute_corrected_asymmetry(ucna.data.raw, ucna.correction_histogram);

	// fit the Fierz term from the asymmetry
	/*
	char A_fit_str[1024];
    sprintf(A_fit_str, "[0]/(1+[1]*(%f/(%f+x)))", m_e, m_e);
    TF1 *A_fierz_fit = new TF1("A_fierz_fit", A_fit_str, min_E, max_E);
    A_fierz_fit->SetParameter(-.12,0);
    A_fierz_fit->SetParameter(1,0);
	asymmetry_histogram->Fit(A_fierz_fit, "Sr");
	asymmetry_histogram->SetMaximum(0);
	asymmetry_histogram->SetMinimum(-0.2);
	*/

	// compute chi squared
	/*
    double chisq = A_fierz_fit->GetChisquare();
    double NDF = A_fierz_fit->GetNDF();
	char A_str[1024];
	sprintf(A_str, "A = %1.3f #pm %1.3f", A_fierz_fit->GetParameter(0), A_fierz_fit->GetParError(0));
	char A_b_str[1024];
	sprintf(A_b_str, "b = %1.3f #pm %1.3f", A_fierz_fit->GetParameter(1), A_fierz_fit->GetParError(1));
	char A_b_chisq_str[1024];
    printf("Chi^2 / ( NDF - 1) = %f / %f = %f\n", chisq, NDF-1, chisq/(NDF-1));
	sprintf(A_b_chisq_str, "#frac{#chi^{2}}{n-1} = %f", chisq/(NDF-1));
	*/

	// draw the ratio plot
	//asymmetry_histogram->SetStats(0);
    //asymmetry_histogram->Draw();

	// draw a legend on the plot
	/*
    TLegend* asym_legend = new TLegend(0.3,0.85,0.6,0.65);
    asym_legend->AddEntry(asymmetry_histogram, "Asymmetry data", "l");
    asym_legend->AddEntry(A_fierz_fit, "Fierz term fit", "l");
    asym_legend->AddEntry((TObject*)0, A_str, "");
    asym_legend->AddEntry((TObject*)0, A_b_str, "");
    asym_legend->AddEntry((TObject*)0, A_b_chisq_str, "");
    asym_legend->SetTextSize(0.03);
    asym_legend->SetBorderSize(0);
    asym_legend->Draw();
	*/

    TString asymmetry_pdf_filename = "mc/asymmetry_data.pdf";
    canvas->SaveAs(asymmetry_pdf_filename);

    /// Compute the super sums
    // TODO make seperate from the loaded one: 
    // ucna.data.super_sum.histogram = compute_super_sum(ucna.data.raw);
    //normalize(ucna.data.super_sum.histogram, min_E, max_E);
    ucna.data.super_sum.histogram->SetLineColor(2);
	ucna.data.super_sum.histogram->SetStats(0);
    ucna.data.super_sum.histogram->Draw("");

    /// Draw Monte Carlo
    ucna.sm.super_sum.histogram->SetLineColor(1);
    ucna.sm.super_sum.histogram->SetMarkerStyle(4);
    ucna.sm.super_sum.histogram->Draw("same p0");

    /// make a pretty legend
    TLegend * legend = new TLegend(0.6,0.8,0.7,0.6);
    legend->AddEntry(ucna.data.super_sum.histogram, "Type 0 super sum", "l");
    legend->AddEntry(ucna.sm.super_sum.histogram, "Monte Carlo super sum", "p");
    //legend->AddEntry(bonehead_sum_histogram, "Bonehead sum", "l");
    legend->SetTextSize(0.03);
    legend->SetBorderSize(0);
    legend->Draw();

    /// save the data and Mote Carlo plots
    TString super_sum_pdf_filename = "mc/super_sum_data.pdf";
    canvas->SaveAs(super_sum_pdf_filename);


    /// CODE BREAK


	//double cov[nPar][nPar]; 
    TMatrixD cov(nPar,nPar);
	double entries = ucna.data.super_sum.histogram->GetEffectiveEntries();
	double N = GetEntries(ucna.data.super_sum.histogram, min_E, max_E);

	/// set all expectation values for this range
    double nSpec = 4;
    TMatrixD expected(nSpec,nSpec);
	for (int m=0; m<nSpec; m++)
		for (int n=0; n<nSpec; n++)
			expected[m][n] = evaluate_expected_fierz(m,n,min_E,max_E);
	
	/// find the predicted inverse covariance matrix for this range
	double A = -0.12;
	TMatrixD p_cov_inv(nPar,nPar);
	for (int i=0; i<nPar; i++)
		for (int j=0; j<nPar; j++)
	        p_cov_inv[i][j] = 0;
    if (nPar > 0)
	    p_cov_inv[0][0] =  N/4*expected[2][0];
    if (nPar > 1) {
        p_cov_inv[1][0] = 
        p_cov_inv[0][1] = -N*A/4*expected[2][1];
        p_cov_inv[1][1] =  N*(expected[0][2] - expected[0][1]*expected[0][1]);
    }
    if (nPar > 2)
	    p_cov_inv[2][2] =  N;

	/// find the covariance matrix
	double det = 0;
	//double p_var_A = 1/p_cov_inv[0][0];
	//double p_var_b = 1/p_cov_inv[1][1];
	TMatrixD p_cov = p_cov_inv.Invert(&det);

	/// actually do the fitting
	TF1* func = combined_fit(ucna.data.asymmetry.histogram, 
                             ucna.data.super_sum.histogram, cov, 0);

	/// output the data info
    cout<<setprecision(5);
	cout<<" ENERGY RANGE:\n";
	cout<<"    Energy range is "<<min_E<<" - "<<max_E<<" keV.\n";
	cout<<"    Number of counts in full data is "<<(int)entries<<".\n";
	cout<<"    Number of counts in energy range is "<<(int)N<<".\n";
	cout<<"    Efficiency energy cut is "<< N/entries*100<<"%.\n";

	/// output the fit covariance details
	cout<<"\n FIT COVARIANCE MATRIX\n";
	for (int i=0; i<nPar; i++) {
        cout<<"    ";
		for (int j=0; j<nPar; j++)
			cout<<setw(14)<<cov[i][j];
	    cout<<endl;
	}

    /// get the fit standard errors
	//double sig_A = sqrt(cov[0][0]);
	//double sig_b = sqrt(cov[1][1]);

	/// output the predicted covariance details	
	cout<<"\n PREDICTED COVARIANCE MATRIX\n";
	for (int i=0; i<nPar; i++) {
        cout<<"    ";
		for (int j=0; j<nPar; j++)
			cout<<setw(14)<<p_cov[i][j];
		cout<<"\n";
	}

    /// get the predicted standard errors
	//double p_sig_A = sqrt(p_cov[0][0]);
	//double p_sig_b = sqrt(p_cov[1][1]);

    /*
	cout<<endl;
	cout<<" FOR UNCOMBINED FITS:\n";
	cout<<"    Expected statistical error for A in this range is without b is " 
					<< sqrt(p_var_A)<<endl;
	cout<<"    Expected statistical error for b in this range is without A is "
					<< sqrt(p_var_b)<<endl;
	cout<<endl;
	cout<<" FOR COMBINED FITS:\n";
	cout<<"    Expected statistical error for A in this range is "<<p_sig_A<<endl;
	cout<<"    Actual statistical error for A in this range is "<<sig_A<<endl;
	cout<<"    Ratio for A error is "<<sig_A / p_sig_A<<endl;
	cout<<"    Expected statistical error for b in this range is "<<p_sig_b<<endl;
	cout<<"    Actual statistical error for b in this range is "<<sig_b<<endl;
	cout<<"    Ratio for b error is "<<sig_b / p_sig_b<<endl;
	cout<<"    Expected cor(A,b) = "<<p_cov[1][0] / (p_sig_A * p_sig_b)<<endl;
	cout<<"    Actual cor(A,b) = "<<cov[1][0] / sqrt(cov[0][0] * cov[1][1])<<endl;
    */

    /// Compute independent standard errors
	cout<<"\n FOR UNCOMBINED FITS:\n";
	for (int i=0; i<nPar; i++) {
        TString name = func->GetParName(i);
        double sigma = 1/sqrt(p_cov_inv[i][i]);
	    cout<<"    Expected independent statistical error for "<<name<<" is "<<sigma<<".\n";
    }

    /// Compare predicted and actual standard errors
	cout<<"\n FOR COMBINED FITS:\n";
	for (int i=0; i<nPar; i++) {
        TString name = func->GetParName(i);
        //double param = func->GetParameter(i);
	    double sigma = sqrt(cov[i][i]);
	    double expected_sigma = sqrt(p_cov[i][i]);
        double factor = sigma/expected_sigma;
        cout<<"    Expected statistical error for "<<name<<" in this range is "<<expected_sigma<<".\n";
        cout<<"    Actual statistical error for "<<name<<" in this range is "<<sigma<<".\n";
        cout<<"    Ratio for "<<name<<" error is "<<factor<<".\n";
    }

    /// Compare predicted and actual correlations
	cout<<"\n CORRELATIONS FACTORS FOR COMBINED FITS:\n";
	for (int i=0; i<nPar; i++) {
        TString name_i = func->GetParName(i);
	    for (int j = i+1; j<nPar; j++) {
            TString name_j = func->GetParName(j);
	        double p_cor_ij = p_cov[j][i]/sqrt(p_cov[i][i]*p_cov[j][j]);
            double cor_ij = cov[j][i]/sqrt(cov[i][i]*cov[j][j]);
            cout<<"    Expected cor("<<name_i<<","<<name_j<<") = "<<p_cor_ij<<".\n";
            cout<<"    Actual cor("<<name_i<<","<<name_j<<") = "<<cor_ij<<".\n";
        }
    }

	/*
	// A fit histogram for output to gnuplot
    TH1D *fierz_fit_histogram = new TH1D(*asymmetry_histogram);
	for (int i = 0; i < fierz_fit_histogram->GetNbinsX(); i++)
		fierz_fit_histogram->SetBinContent(i, fierz_fit->Eval(fierz_fit_histogram->GetBinCenter(i)));
	

	/// output for root
    TString pdf_filename = "/data/kevinh/mc/asymmetry_fierz_term_fit.pdf";
    canvas->SaveAs(pdf_filename);

	/// output for gnuplot
	//output_histogram("/data/kevinh/mc/super-sum-data.dat", ucna.data.super_sum.histogram, 1, 1000);
	//output_histogram("/data/kevinh/mc/super-sum-mc.dat", ucna.sm.super_sum.histogram, 1, 1000);
	//output_histogram("/data/kevinh/mc/fierz-ratio.dat", fierz_ratio_histogram, 1, 1);
	//output_histogram("/data/kevinh/mc/fierz-fit.dat", fierz_fit_histogram, 1, 1);

	*/

#if 0
	// Create a new canvas.
	TCanvas * c1 = new TCanvas("c1","Two Histogram Fit example",100,10,900,800);
	c1->Divide(2,2);
	gStyle->SetOptFit();
	gStyle->SetStatY(0.6);

	c1->cd(1);
	ucna.data.asymmetry.histogram->Draw();
	func->SetRange(min_E, max_E);
	func->DrawCopy("cont1 same");
	/*
	c1->cd(2);
	asymmetry->Draw("lego");
	func->DrawCopy("surf1 same");
	*/
	c1->cd(3);
	func->SetRange(min_E, max_E);
	///TODO fierzratio_histogram->Draw();
	func->DrawCopy("cont1 same");
	/*
	c1->cd(4);
	fierzratio->Draw("lego");
	gPad->SetLogz();
	func->Draw("surf1 same");
	*/
#endif
	app.Run();

	return 0;
}
