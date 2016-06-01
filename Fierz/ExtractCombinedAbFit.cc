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

/// name spaces
using std::setw;
using std::cout;
using namespace TMath;

/// cuts and settings
double KEmin = 0;                 /// min kinetic energy for plots
double KEmax = 800;               /// max kinetic range for plots
int KEbins=(KEmax-KEmin)/10;      /// number of bins to use fit spectral plots
double fit_min = 120;             /// min kinetic energy for plots
double fit_max = 630;             /// max kinetic range for plots
int fit_bins=(fit_max-fit_min)/10;/// number of bins to use fit spectral plots
double fedutial_cut = 50;         /// radial cut in millimeters TODO!! HARD CODED IN MODEL


/// set up free fit parameters with best guess
static const int nPar = 3;
TString paramNames[3] = {"A", "b", "N"};
double paramInits[3] = {-0.12, 0, 1};

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
UCNAFierzFitter ucna(KEbins, KEmin, KEmax, fit_bins, fit_min, fit_max);
UCNAFierzFitter fake(KEbins, KEmin, KEmax, fit_bins, fit_min, fit_max);
//ucna.fidcut2 = fedutial_cut*fedutial_cut;

/// This needs to be static and global for MINUIT to work
void combined_chi2(Int_t & n, Double_t * /*grad*/ , Double_t &chi2, Double_t *p, Int_t /*iflag */  )
{
    assert(n==3);
    double A=p[0], b=p[1], N=p[2]; // TODO make nPar correct here
	chi2 = ucna.combined_chi2(A,b,N);
}





/// MAIN APPLICATION
int main(int argc, char *argv[])
{
	TApplication app("Extract Combined A+b Fitter", &argc, argv);
	srand( time(NULL) );    /// set this to make random or repeatable

    /// LOAD 2010 UCNA DATA

    /// Load the files that already contain data super histogram.
    /*
    for (int side=EAST; side<=WEST; side++)
        for (int afp=EAST; afp<=WEST; afp++) {
            TString s = side? "W":"E", a = afp? "On":"Off";
            TString title = "2010 final official "+s+" afp "+a;
            TString name = "hTotalEvents_"+s+"_"+a+";1";
            int entries = ucna.data.counts[side][afp]->fill(data_dir+"OctetAsym_Offic.root", name, title);
            if (entries) {
                cout<<"Status: Number of entries in "<<(side? "west":"east")
                    <<" side with afp "<<a<<" is "<<entries<<".\n";
            }
            else
                cout<<"Error: found no events in "<<title<<".\n";
        }
        */

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


    /// LOAD FAKE DATA FROM MONTE CARLO

    /// Load Monte Carlo simulated Standard Model events
    fake.sm.fill(mc_syst_dir+"SimAnalyzed_2010_Beta_paramSet_100_0.root",
                 "SimAnalyzed",
                 "Monte Carlo Standard Model beta spectrum");

    /// Load Monte Carlo simulated Fierz events
    fake.fierz.fill(mc_syst_dir+"SimAnalyzed_2010_Beta_fierz_paramSet_100_0.root",
                    "SimAnalyzed",
                    "Monte Carlo Fierz beta spectrum");

    /// For now load real asymmetry data as fake histogram.
    fake.data.asymmetry.fill(
        data_dir+"Range_0-1000/CorrectAsym/CorrectedAsym.root",
        "hAsym_Corrected_C",
        "2010 final official asymmetry");

    /// Just overwrite
    //fake

    /// LOAD MONTE CARLO SIMULATION EVENTS

    /// Load Monte Carlo simulated Standard Model events
    ucna.sm.fill(mc_dir+"SimAnalyzed_Beta_7.root",
                 "SimAnalyzed",
                 "Monte Carlo Standard Model beta spectrum");

    /// Load Monte Carlo simulated Fierz events
    ucna.fierz.fill(mc_dir+"SimAnalyzed_Beta_fierz_7.root",
                    "SimAnalyzed",
                    "Monte Carlo Fierz beta spectrum");

    /// SAVE ALL PROCESSED HISTOGRAMS 

    // TODO ???


    /// FITTING 

    /// Set up the fit parameters.
    TF1 asymmetry_func("asymmetry_fit_func", &asymmetry_fit_func, fit_min, fit_max, nPar);
    //TF1 supersum_func("supersum_fit_func", &supersum_fit_func, KEmin_b, KEmax_b, nPar);
    asymmetry_func.SetParameters(paramInits);
    //super_sum_func.SetParameters(paramInits);
    for (int i=0; i<nPar; i++) {
        asymmetry_func.SetParName(i, paramNames[i]);
        //super_sum_func.SetParName(i, paramNames[i]);
    }

    /// Set up constants and vars
    TMatrixD cov(nPar,nPar);
    double all_entries = ucna.data.super_sum.GetEntries();
	double eff_entries = ucna.data.super_sum.GetEffectiveEntries(KEmin, KEmax);
	double fit_entries = ucna.data.super_sum.GetEffectiveEntries(fit_min, fit_max);

	/// Actually do the combined fitting.
	//ucna.combined_fit(cov, &asymmetry_func, &supersum_func, &combined_chi2);
	ucna.combined_fit(cov, &asymmetry_func, &combined_chi2);
	ucna.compute_fit(&asymmetry_func);

    /// Set up reasonable guesses 
    double A = -0.12;
    //double b = 0;
    double N = fit_entries;
    double exCos = 0.25; //evaluate_expected_cos_theta(fit_min,fit_max);

	/// Set all expectation values for this range.
    double nSpec = 4;
    TMatrixD ex(nSpec,nSpec);
	for (int m=0; m<nSpec; m++)
		for (int n=0; n<nSpec; n++)
			ex[m][n] = evaluate_expected_fierz(m,n,fit_min,fit_max);
	
	/// Calculate the predicted inverse covariance matrix for this range.
	TMatrixD p_cov_inv(nPar,nPar);
	for (int i=0; i<nPar; i++)
		for (int j=0; j<nPar; j++)
	        p_cov_inv[i][j] = 0;
    if (nPar > 0)
	    p_cov_inv[0][0] = N*exCos*ex[2][0];
    if (nPar > 1) {
        p_cov_inv[1][0] = 
        p_cov_inv[0][1] = -N*A*exCos*ex[2][1];
        //p_cov_inv[1][1] =  N*(A*A*ex[2][2]/4 + ex[0][2] - ex[0][1]*ex[0][1]);
        p_cov_inv[1][1] = N*(A*A*exCos*ex[2][2] + ex[0][2] - ex[0][1]*ex[0][1]);
    }
    if (nPar > 2)
	    p_cov_inv[2][2] = N;

	/// Compute the predicted covariance matrix by inverse.
	double det = 0;
	TMatrixD p_cov = p_cov_inv.Invert(&det);



    /// PRINT OUT REPORT OF FITS, CORRELATIONS AND ERRORS

	/// Output the data info.
    int cl = 14;
    cout<<setprecision(5);
	cout<<" ENERGY RANGE:\n";
	cout<<"    Full energy range is "<<KEmin<<" - "<<KEmax<<" keV.\n";
	cout<<"    Fit energy range cut is "<<fit_min<<" - "<<fit_max<<" keV.\n";
	cout<<"    Number of counts in data is "<<int(all_entries)<<".\n";
	cout<<"    Effective number of counts in full energy range is "<<int(eff_entries)<<".\n";
	cout<<"    Effective number of counts in fit energy range cut is "<<int(fit_entries)<<".\n";
	cout<<"    Efficiency of energy cut is "<< int(fit_entries/eff_entries*1000)/10<<"%.\n";
	cout<<"    Number of fit counts in full range is "<<int(N)<<".\n";

	/// Output the fit covariance details.
	cout<<"\n FIT COVARIANCE MATRIX\n";
    cout<<"     ";
	for (int i=0; i<nPar; i++)
		cout<<setw(cl)<<paramNames[i];
	cout<<"\n";
	for (int i=0; i<nPar; i++) {
        cout<<"    "<<paramNames[i];
		for (int j=0; j<nPar; j++)
			cout<<setw(cl)<<cov[i][j];
	    cout<<"\n";
	}

	/// Output the predicted covariance details.
	cout<<"\n PREDICTED COVARIANCE MATRIX\n";
    cout<<"     ";
	for (int i=0; i<nPar; i++)
		cout<<setw(cl)<<paramNames[i];
	cout<<"\n";
	for (int i=0; i<nPar; i++) {
        cout<<"    "<<paramNames[i];
		for (int j=0; j<nPar; j++)
			cout<<setw(cl)<<p_cov[i][j];
		cout<<"\n";
	}

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
        double sigma = 1/Sqrt(p_cov_inv[i][i]);
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
                <<setw(cl+1)<<"actual/est.\n";
	for (int i=0; i<nPar; i++) {
        TString param = paramNames[i];
        double value = asymmetry_func.GetParameter(i);
	    double sigma = Sqrt(cov[i][i]);
	    double error = 100*sigma/value;
	    double est_sigma = Sqrt(p_cov[i][i]);
	    double est_error = 100*est_sigma/value;
        double ratio = sigma/est_sigma;
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
                <<setw(cl+1)<<"actual/est.\n";
	for (int i=0; i<nPar; i++) {
        TString name_i = paramNames[i];
	    for (int j = i+1; j<nPar; j++) {
            TString name_j = paramNames[j];
            TString cor_name = "cor("+name_i+","+name_j+")";
            double cor_ij = cov[j][i]/Sqrt(cov[i][i]*cov[j][j]);
	        double p_cor_ij = p_cov[j][i]/Sqrt(p_cov[i][i]*p_cov[j][j]);
            double ratio = cor_ij/p_cor_ij;
            cout<<"    "<<setw(8)<<cor_name
                        <<setw(cl)<<cor_ij
                        <<setw(cl)<<p_cor_ij
                        <<setw(cl)<<ratio<<"\n";
        }
    }


    /// DISPLAYING AND OUTPUTTING
    TCanvas canvas("fierz_fitter_canvas",
                   "Combined Fierz component of energy spectrum");
    TLegend legend(0.55,0.65,0.85,0.85);

    ucna.data.asymmetry.draw(
                   "data_asymmetry", 
                   "UCNA 2010 #Lambda(E)", 
                   &canvas,&legend,"",1,0);
    ucna.fit.asymmetry.draw(
                   "fit_asymmetry", 
                   "Fit #Lambda(E)", 
                   &canvas,&legend,"Same",6,0);
    canvas.SaveAs(plots_dir+"data_asymmetry.pdf");

    ucna.sm.super_sum.draw(
                   "sm_supersum", 
                   "Standard Model Monte Carlo #Sigma", 
                   &canvas,&legend,"",4,0);
    ucna.fierz.super_sum.draw(
                   "fierz_supersum", 
                   "Fierz Monte Carlo #Sigma", 
                   &canvas,&legend,"Same",6,0);
    canvas.SaveAs(plots_dir+"monte_carlo.pdf");

    ucna.data.super_sum.draw(
                   "data_supersum", 
                   "UCNA 2010 #Sigma(E)", 
                   &canvas,&legend,"",1,0);
    ucna.fit.super_sum.draw(
                   "fit_supersum", 
                   "Fit #Sigma(E)", 
                   &canvas,&legend,"Same",6,0);
    canvas.SaveAs(plots_dir+"data_supersum.pdf");

    /*
	/// Output for gnuplot
	//save_data(plots_dir+"super-sum-data.dat", ucna.data.super_sum, 1, 1000);
	//save_data(plots_dir+"super-sum-mc.dat", ucna.sm.super_sum, 1, 1000);
	//save_data(plots_dir+"fierz-ratio.dat", fierz_ratio_histogram, 1, 1);
	//save_data(plots_dir+"fierz-fit.dat", fierz_fit_histogram, 1, 1);

	*/

	app.Run();

	return 0;
}




/// OLD BUT USEFUL... LIKE SANTA
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
    ucna.data.counts[0][0]=(TH1D*)ucna_data_tfile->Get("hTotalEvents_E_Off;1");
    ucna.data.counts[0][1]=(TH1D*)ucna_data_tfile->Get("hTotalEvents_E_On;1");
    ucna.data.counts[1][0]=(TH1D*)ucna_data_tfile->Get("hTotalEvents_W_Off;1");
    ucna.data.counts[1][1]=(TH1D*)ucna_data_tfile->Get("hTotalEvents_W_On;1");
    */

    /* TODO figure out where these went.
    fill_data("OctetAsym_Offic.root",
              "2010 final official east afp off spectrum",
              "hTotalEvents_E_off;1",
              ucna.data.counts[0][0]);
    fill_data("OctetAsym_Offic.root",
              "2010 final official east afp on spectrum",
              "hTotalEvents_E_on;1",
              ucna.data.counts[0][1]);
    fill_data("OctetAsym_Offic.root",
              "2010 final official west afp off spectrum",
              "hTotalEvents_W_off;1",
              ucna.data.counts[1][0]);
    fill_data( "OctetAsym_Offic.root",
              "2010 final official west afp on spectrum",
              "hTotalEvents_W_on;1",
              ucna.data.counts[1][1]);


    /// TODO figure out where these went.
    for (int side=EAST; side<=WEST; side++)
        for (int afp=EAST; afp<=WEST; afp++) {
            TString sw = side? "west":"east", s = side? "W":"E", a = afp? "on":"off";
            TString title = "2010 final official "+s+" afp "+a;
            TString cut = "hTotalEvents_"+s+"_"+a+";1";
            int entries = fill_data("OctetAsym_Offic.root", title, cut, ucna.data.counts[side][afp]);
            if (entries) 
                cout<<"Status: Number of entries in "<<sw<<" side with afp "<<a<<" is "<<entries<<".\n";
            else
                cout<<"Error: found no events in "<<title<<".\n";
        }

    Already background subtracted...
        TH1D *background_histogram = (TH1D*)ucna_data_tfile->Get("Combined_Events_E000");
        ucna_data.counts->Add(background_histogram,-1);
        // normalize after background subtraction
        background_histogram->Draw("");
    */

    /*
	for (int side = 0; side < 2; side++)
		for (int spin = 0; spin < 2; spin++)
		{
			cout << "Number of entries in (" 
					  << side << ", " << spin << ") is "
					  << (int)ucna.data.counts[side][spin]->GetEntries() << endl;
			if (ucna.data.counts[side][spin] == NULL)
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
    //TH1D *ucna.correction_histogram = new TH1D(*ucna.data.counts[0][0]);
	/*
	while (tfile has more entries)
	{
        double bin = (ucna.correction_file->GetBinContent(bin);
        double correction = ucna.correction_histogram->GetBinContent(bin);
        printf("Setting bin content for correction bin %d, to %f\n", bin, correction);
        ucna.data.super_sum.SetBinContent(bin, correction);
        ucna.data.super_sum.SetBinError(bin, correction_error);
    }
	*/

    /*
    TF1 *fit = new TF1("fierz_fit", theoretical_fierz_spectrum, 0, 1000, 3);
    fit->SetParameter(0,0.0);
    fit->SetParameter(1,0.0);
    fit->SetParameter(2,1.0);
    ucna.data.counts[0][0]->Fit("fierz_fit");
    double chisq = fit->GetChisquare();
    double N = fit->GetNDF();
    printf("Chi^2 / ( N - 1) = %f / %f = %f\n",chisq, N-1, chisq/(N-1));

    TString fit_pdf_filename = "mc/fierz_fit_data.pdf";
    canvas->SaveAs(fit_pdf_filename);

    // compute and plot the super ratio
    TH1D *ucna.data.super_ratio = compute_super_ratio(ucna.data.counts);
    ucna.data.super_ratio.Draw();
    TString super_ratio_pdf_filename = "mc/super_ratio_data.pdf";
    canvas->SaveAs(super_ratio_pdf_filename);
    */

    // compute and plot the super ratio asymmetry 
    //TH1D *asymmetry_histogram = compute_corrected_asymmetry(ucna.data.counts, ucna.correction_histogram);

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
    /*/*
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
    ucna.data.counts[0][0]=(TH1D*)ucna_data_tfile->Get("hTotalEvents_E_Off;1");
    ucna.data.counts[0][1]=(TH1D*)ucna_data_tfile->Get("hTotalEvents_E_On;1");
    ucna.data.counts[1][0]=(TH1D*)ucna_data_tfile->Get("hTotalEvents_W_Off;1");
    ucna.data.counts[1][1]=(TH1D*)ucna_data_tfile->Get("hTotalEvents_W_On;1");
    */

    /* TODO figure out where these went.
    fill_data("OctetAsym_Offic.root",
              "2010 final official east afp off spectrum",
              "hTotalEvents_E_off;1",
              ucna.data.counts[0][0]);
    fill_data("OctetAsym_Offic.root",
              "2010 final official east afp on spectrum",
              "hTotalEvents_E_on;1",
              ucna.data.counts[0][1]);
    fill_data("OctetAsym_Offic.root",
              "2010 final official west afp off spectrum",
              "hTotalEvents_W_off;1",
              ucna.data.counts[1][0]);
    fill_data( "OctetAsym_Offic.root",
              "2010 final official west afp on spectrum",
              "hTotalEvents_W_on;1",
              ucna.data.counts[1][1]);


    /// TODO figure out where these went.
    for (int side=EAST; side<=WEST; side++)
        for (int afp=EAST; afp<=WEST; afp++) {
            TString sw = side? "west":"east", s = side? "W":"E", a = afp? "on":"off";
            TString title = "2010 final official "+s+" afp "+a;
            TString cut = "hTotalEvents_"+s+"_"+a+";1";
            int entries = fill_data("OctetAsym_Offic.root", title, cut, ucna.data.counts[side][afp]);
            if (entries) 
                cout<<"Status: Number of entries in "<<sw<<" side with afp "<<a<<" is "<<entries<<".\n";
            else
                cout<<"Error: found no events in "<<title<<".\n";
        }

    Already background subtracted...
        TH1D *background_histogram = (TH1D*)ucna_data_tfile->Get("Combined_Events_E000");
        ucna_data.counts->Add(background_histogram,-1);
        // normalize after background subtraction
        background_histogram->Draw("");
    */

    /*
	for (int side = 0; side < 2; side++)
		for (int spin = 0; spin < 2; spin++)
		{
			cout << "Number of entries in (" 
					  << side << ", " << spin << ") is "
					  << (int)ucna.data.counts[side][spin]->GetEntries() << endl;
			if (ucna.data.counts[side][spin] == NULL)
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
    //TH1D *ucna.correction_histogram = new TH1D(*ucna.data.counts[0][0]);
	/*
	while (tfile has more entries)
	{
        double bin = (ucna.correction_file->GetBinContent(bin);
        double correction = ucna.correction_histogram->GetBinContent(bin);
        printf("Setting bin content for correction bin %d, to %f\n", bin, correction);
        ucna.data.super_sum.SetBinContent(bin, correction);
        ucna.data.super_sum.SetBinError(bin, correction_error);
    }
	*/

    /*
    TF1 *fit = new TF1("fierz_fit", theoretical_fierz_spectrum, 0, 1000, 3);
    fit->SetParameter(0,0.0);
    fit->SetParameter(1,0.0);
    fit->SetParameter(2,1.0);
    ucna.data.counts[0][0]->Fit("fierz_fit");
    double chisq = fit->GetChisquare();
    double N = fit->GetNDF();
    printf("Chi^2 / ( N - 1) = %f / %f = %f\n",chisq, N-1, chisq/(N-1));

    TString fit_pdf_filename = "mc/fierz_fit_data.pdf";
    canvas->SaveAs(fit_pdf_filename);

    // compute and plot the super ratio
    TH1D *ucna.data.super_ratio = compute_super_ratio(ucna.data.counts);
    ucna.data.super_ratio.Draw();
    TString super_ratio_pdf_filename = "mc/super_ratio_data.pdf";
    canvas->SaveAs(super_ratio_pdf_filename);
    */

    // compute and plot the super ratio asymmetry 
    //TH1D *asymmetry_histogram = compute_corrected_asymmetry(ucna.data.counts, ucna.correction_histogram);

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
