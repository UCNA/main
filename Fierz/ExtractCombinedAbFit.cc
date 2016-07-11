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
#include <limits>


/// C includes
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/// name spaces
using std::setw;
using std::cout;
using namespace TMath;
double infity = std::numeric_limits<double>::infinity();

/// cuts and settings
double KEmin = 0;                 /// min kinetic energy for plots
double KEmax = 800;               /// max kinetic range for plots
int KEbins=(KEmax-KEmin)/10;      /// number of bins to use fit spectral plots
double fit_min = 120;             /// min kinetic energy for plots
double fit_max = 630;             /// max kinetic range for plots
int fit_bins=(fit_max-fit_min)/10;/// number of bins to use fit spectral plots
double fedutial_cut = 50;         /// radial cut in millimeters TODO!! HARD CODED IN MODEL

/// set up free fit parameters with best guess

#if 0
static TString FIT_TYPE = "AbN";
static const int nPar = 3;
double afp_ratio = 0.40;
TString paramNames[3] = {"A", "b", "N"};
double paramInits[3] = {-0.12, 0.0, 1};
#else
static TString FIT_TYPE = "bN";
static const int nPar = 2;
double afp_ratio = 0.40;
TString paramNames[2] = {"b", "N"};
double paramInits[2] = {0.0, 1};
int A_index=-1, b_index=0, N_index=1;
#endif

/// path to experiment data files
TString data_dir = "/media/hickerson/boson/Data/OctetAsym_Offic_2010_FINAL/"; 

/// path to Monte Carlo files
TString mc_dir = "/home/xuansun/Documents/SimProcessedFiles/100mill_beta/";
TString mc_syst_dir = "/home/xuansun/Documents/SimProcessedFiles/100mill_both_twiddled/";

/// path to save output plots
TString plots_dir = "/home/hickerson/Dropbox/Root/";

/// path to save output root structures 
TString root_output_dir = "/home/hickerson/Documents/";




/// GLOBAL MODELS
//ucna.fidcut2 = fedutial_cut*fedutial_cut;

/// This needs to be static and global for MINUIT to work
UCNAFierzFitter* global_ff = 0;

void combined_chi2(Int_t & n, Double_t * /*grad*/ , Double_t &chi2, Double_t *p, Int_t /*iflag */  )
{
    assert(global_ff);
    if (FIT_TYPE=="AbN") {
        double A=p[0], b=p[1], N=p[2];
	    chi2 = global_ff->combined_chi2(A,b,N);
    }
    else if (FIT_TYPE=="bN") {
        double b=p[0], N=p[1];
        chi2 = global_ff->supersum_chi2(b,N);
    }
}


int cl = 14;
void output_matrix(TString title, TMatrixD matrix)
{
	cout<<"\n "<<title<<"\n";
    cout<<"     ";
	for (int i=0; i<nPar; i++)
		cout<<setw(cl)<<paramNames[i];
	cout<<"\n";
	for (int i=0; i<nPar; i++) {
        cout<<"    "<<paramNames[i];
		for (int j=0; j<nPar; j++)
			cout<<setw(cl)<<matrix[i][j];
	    cout<<"\n";
	}
}


/// FITTING 
void fit(UCNAFierzFitter &ff) {
    /// Set up the fit functions parameters.
    TF1 asymmetry_func("asymmetry_fit_func", &asymmetry_fit_func, fit_min, fit_max, nPar);
    asymmetry_func.SetParameters(paramInits);
    for (int i=0; i<nPar; i++)
        asymmetry_func.SetParName(i, paramNames[i]);

    /*TF1 supersum_func("supersum_fit_func", &supersum_fit_func, fit_min, fit_max, nPar);
    supersum_func.SetParameters(paramInits);
    for (int i=0; i<nPar; i++)
        supersum_func.SetParName(i, paramNames[i]);*/

	/// Actually do the combined fitting.
    global_ff = &ff; /// TODO cannot be called in a member
    TMatrixD cov(nPar,nPar);
	ff.combined_fit(cov, &asymmetry_func, /*supersum_func,*/ &combined_chi2); /// TODO no member
	ff.compute_fit(&asymmetry_func);
    global_ff = 0;

    /// Look up sizes
    double all_entries = ff.data.super_sum.GetEntries();
	double eff_entries = ff.data.super_sum.GetEffectiveEntries(KEmin, KEmax);
	double fit_entries = ff.data.super_sum.GetEffectiveEntries(fit_min, fit_max);

    /// Set up reasonable guesses 
    double A = -0.12;
    double b = 0;
    double N = fit_entries;
    double ex_cos = 0.5; //evaluate_expected_cos_theta(fit_min,fit_max);

    /// PRINT OUT REPORT OF FITS, CORRELATIONS AND ERRORS

	/// Output the data and cut info.
    cout<<setprecision(5);
	cout<<" ENERGY RANGE:\n";
	cout<<"    Full energy range is "<<KEmin<<" - "<<KEmax<<" keV.\n";
	cout<<"    Fit energy range cut is "<<fit_min<<" - "<<fit_max<<" keV.\n";
	cout<<"    Number of counts in data is "<<int(all_entries)<<".\n";
	cout<<"    Effective number of counts in full energy range is "<<int(eff_entries)<<".\n";
	cout<<"    Effective number of counts in fit energy range cut is "<<int(fit_entries)<<".\n";
	cout<<"    Efficiency of energy cut is "<< int(fit_entries/eff_entries*1000)/10<<"%.\n";

	/// Set all expectation values for this range.
    double nSpec = 4;
    TMatrixD ex(nSpec,nSpec);
	for (int m=0; m<nSpec; m++)
		for (int n=0; n<nSpec; n++)
			ex[m][n] = evaluate_expected_fierz(m,n,fit_min,fit_max);
	
	/// Calculate the predicted inverse covariance matrix for this range.
	TMatrixD est_cov_inv(nPar,nPar);
	for (int i=0; i<nPar; i++)
		for (int j=0; j<nPar; j++)
	        est_cov_inv[i][j] = 0;

    if (A_index >= 0)
        est_cov_inv[A_index][A_index] += N*(A*A*ex_cos*ex[2][2]);
    if (b_index >= 0) {
	    est_cov_inv[b_index][b_index] = N*ex_cos*ex[2][0];
        if (A_index >= 0) {
            est_cov_inv[A_index][b_index] = 
            est_cov_inv[b_index][A_index] = -N*A*ex_cos*ex[2][1];
            est_cov_inv[A_index][A_index] += N*(ex[0][2] - ex[0][1]*ex[0][1]);
        }
    }
    if (N_index >= 0) {
	    est_cov_inv[N_index][N_index] = N;
        if (A_index >= 0) {
            est_cov_inv[A_index][N_index] =
            est_cov_inv[N_index][A_index] = N*ex[0][1];
        }
    }

/*  if (b_index > -1)
	    est_cov_inv[0][0] = N*ex_cos*ex[2][0];
    if (nPar > 1) {
        est_cov_inv[1][0] = 
        est_cov_inv[0][1] = -N*A*ex_cos*ex[2][1];
        //est_cov_inv[1][1] =  N*(A*A*ex[2][2]/4 + ex[0][2] - ex[0][1]*ex[0][1]);
        est_cov_inv[1][1] = N*(A*A*ex_cos*ex[2][2] + ex[0][2] - ex[0][1]*ex[0][1]);
    }
    if (nPar > 2) {
	    est_cov_inv[1][2] =
	    est_cov_inv[2][1] = N*ex[0][1];
	    est_cov_inv[2][2] = N;
    } */

	/// Output the fit covariance matrix from the fit.
	output_matrix("FIT COVARIANCE MATRIX", cov);

	/// Compute the estimated covariance matrix by inverse.
	double det = 0;
	TMatrixD est_cov = est_cov_inv;
	output_matrix("ESTIMATED COVARIANCE MATRIX", est_cov.Invert(&det));

	/// Output the estimated covariance details.
    TMatrixD ratio_cov(nPar,nPar);
	for (int i=0; i<nPar; i++)
		for (int j=0; j<nPar; j++)
			ratio_cov[i][j] = est_cov[i][j]/cov[i][j];
	output_matrix("ESTIMATED/ACTUAL COVARIANCE RATIOS", ratio_cov);

	/// Output the fit covariance inverse matrix time 2.
	TMatrixD cov_inv = cov;
	output_matrix("FIT HALF CHI^2 DERIVATIVE MATRIX", cov_inv.Invert(&det));

	/// Compute the estimated covariance matrix by inverse.
	output_matrix("ESTIMATED HALF CHI^2 DERIVATIVE MATRIX", est_cov_inv);

	/// Output the estimated covariance details.
    TMatrixD ratio_cov_inv(nPar,nPar);
	for (int i=0; i<nPar; i++)
		for (int j=0; j<nPar; j++)
			ratio_cov_inv[i][j] = est_cov_inv[i][j]/cov_inv[i][j];
	output_matrix("ESTIMATED/ACTUAL HALF CHI^2 DERIVATIVE RATIOS", ratio_cov_inv);

#if 0
    /// Compute independent standard errors.
	cout<<"\n FOR UNCOMBINED FITS:\n";
	cout<<"    Expected independent statistical sigma and error:\n";
    cout<<"    "<<setw(1)<<" " 
                <<setw(cl)<<"value"
                <<setw(cl)<<"sigma" 
                <<setw(cl+1)<<"error\n";
	for (int i=0; i<nPar; i++) {
        TString name = paramNames[i];
        double value = asymmetry_func.GetParameter(i);
        double sigma = 1/Sqrt(est_cov_inv[i][i]);
        double error = 100*sigma/value;
        cout<<"    "<<setw(1) <<name
                    <<setw(cl)<<value
                    <<setw(cl)<<sigma
                    <<setw(cl-1)<<error<<"%\n";
    }
#endif

    /// Compare predicted and actual standard errors.
	cout<<"\n FOR COMBINED FITS:\n";
    cout<<"    Actual and estimated combined statistical sigmas and errors:\n";
    cout<<"    "<<setw(1)<<" "
                <<setw(cl)<<"value" 
                <<setw(cl)<<"actual sigma"
                <<setw(cl)<<"est. sigma" 
                <<setw(cl)<<"actual error" 
                <<setw(cl)<<"est. error" 
                <<setw(cl+1)<<"est./actual\n";
	for (int i=0; i<nPar; i++) {
        TString param = paramNames[i];
        double value = asymmetry_func.GetParameter(i);
	    double sigma = Sqrt(cov[i][i]);
	    double error = 100*sigma/value;
	    double est_sigma = Sqrt(est_cov[i][i]);
	    double est_error = 100*est_sigma/value;
        double ratio = est_sigma/sigma;
        cout<<"    "<<setw(1)   <<param 
                    <<setw(cl)  <<value 
                    <<setw(cl)  <<sigma 
                    <<setw(cl)  <<est_sigma 
                    <<setw(cl-1)<<error <<"%" 
                    <<setw(cl-1)<<est_error<<"%"
                    <<setw(cl)  <<ratio<<"\n";
    }

    /// Compare predicted and actual correlations.
	cout<<"\n CORRELATIONS FACTORS FOR COMBINED FITS:\n";
    cout<<"    "<<setw(8)<<" "
                <<setw(cl)<<"actual"
                <<setw(cl)<<"estimate" 
                <<setw(cl+1)<<"est./actual\n";
	for (int i=0; i<nPar; i++) {
        TString name_i = paramNames[i];
	    for (int j = i+1; j<nPar; j++) {
            TString name_j = paramNames[j];
            TString cor_name = "cor("+name_i+","+name_j+")";
            double cor_ij = cov[j][i]/Sqrt(cov[i][i]*cov[j][j]);
	        double est_cor_ij = est_cov[j][i]/Sqrt(est_cov[i][i]*est_cov[j][j]);
            double ratio = est_cor_ij/cor_ij;
            cout<<"    "<<setw(8)<<cor_name
                        <<setw(cl)<<cor_ij
                        <<setw(cl)<<est_cor_ij
                        <<setw(cl)<<ratio<<"\n";
        }
    }
}




/// MAIN APPLICATION
int main(int argc, char *argv[])
{
	TApplication app("Extract Combined A+b Fitter", &argc, argv);
	srand( time(NULL) );    /// set this to make random or repeatable

    UCNAFierzFitter ucna("ucna", "UCNA", KEbins, KEmin, KEmax, fit_bins, fit_min, fit_max);
    UCNAFierzFitter fake("fake", "Fake UCNA", KEbins, KEmin, KEmax, fit_bins, fit_min, fit_max);

    /// LOAD 2010 UCNA DATA

    /// Load the files that already contain data super histogram.
    /*
    for (int side=EAST; side<=WEST; side++)
        for (int afp=EAST; afp<=WEST; afp++) {
            TString s = side? "W":"E", a = afp? "On":"Off";
            TString title = "2010 final official "+s+" afp "+a;
            TString name = "hTotalEvents_"+s+"_"+a+";1";
            int entries = ff.data.counts[side][afp]->fill(data_dir+"OctetAsym_Offic.root", name, title);
            if (entries) {
                cout<<"Status: Number of entries in "<<(side? "west":"east")
                    <<" side with afp "<<a<<" is "<<entries<<".\n";
            }
            else
                cout<<"Error: found no events in "<<title<<".\n";
        }
        */

    /// LOAD PRECOMPUTED HISTOGRAMS AND OVERWRITE 

/*
    /// Load the files that already contain data asymmetry histogram.
    ucna.data.asymmetry.fill(
        data_dir+"Range_0-1000/CorrectAsym/CorrectedAsym.root",
        "hAsym_Corrected_C",
        "2010 final official asymmetry");

    /// Load the files that already contain data super histogram.
    ucna.data.super_sum.fill(
        data_dir+"OctetAsym_Offic.root",
        "Total_Events_SuperSum",
        "2010 final official supersum");

    /// Load Monte Carlo simulated Standard Model events
    ucna.sm.fill(mc_dir+"SimAnalyzed_Beta_7.root",
               "SimAnalyzed",
               "Monte Carlo Standard Model beta spectrum");

    /// Load Monte Carlo simulated Fierz events
    ucna.fierz.fill(mc_dir+"SimAnalyzed_Beta_fierz_7.root",
                  "SimAnalyzed",
                  "Monte Carlo Fierz beta spectrum");

    //fit(ucna);
    //display(ucna);
    */


    /// LOAD FAKE DATA FROM MONTE CARLO
    /// Load Monte Carlo simulated Standard Model events
    fake.vector.fill(
        mc_syst_dir+"SimAnalyzed_2010_Beta_paramSet_100_0.root",
        "SimAnalyzed",
        "Vector Standard Model Monte Carlo beta spectrum", afp_ratio, 0, 0);

    fake.axial.fill(
        mc_syst_dir+"SimAnalyzed_2010_Beta_paramSet_100_0.root",
        "SimAnalyzed",
        "Axial-vector Standard Model Monte Carlo beta spectrum", afp_ratio, 1, 0);

    /// Load Monte Carlo simulated Fierz events
    fake.fierz.fill(
        mc_syst_dir+"SimAnalyzed_2010_Beta_fierz_paramSet_100_0.root",
        "SimAnalyzed",
        "Fierz Monte Carlo beta spectrum", afp_ratio, 0, 1); // TODO this is supressing the errors

    /// For now load real asymmetry data as fake histogram. TODO Fix.
    /// Load Monte Carlo simulated Standard Model events
    double A = -0.12;
    double b = 0;
    fake.data.fill(
        mc_dir+"SimAnalyzed_Beta_7.root",
        "SimAnalyzed",
        "Monte Carlo Standard Model beta spectrum", afp_ratio, A, b);

    /*
    fake.data.asymmetry.fill(
        data_dir+"Range_0-1000/CorrectAsym/CorrectedAsym.root",
        "hAsym_Corrected_C",
        "2010 final official asymmetry");
    */

    fit(fake);
    fake.display(plots_dir);

	app.Run();
	return 0;
}




