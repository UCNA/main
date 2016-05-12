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
//static unsigned nToSim = 5e7;				/// how many triggering events to simulate
//static double afp_off_prob = 1/1.68; 	    /// afp off probability per neutron (0.68/1.68 for on)
int KEbins = 150;                           /// number of bins to use fit spectral plots

double KEmin = 120;                         /// min kinetic energy for plots
double KEmax = 650;                         /// max kinetic range for plots
double KEmin_A = 120;                       /// min kinetic energy for asymmetry fit
double KEmax_A = 650;                       /// max kinetic range for asymmetry fit
double KEmin_b = 120;                       /// min kinetic energy for Fierz fit
double KEmax_b = 650;                       /// max kinetic range for Fierz fit
double fedutial_cut = 50;                   /// radial cut in millimeters 
double fidcut2 = 50*50;                     /// mm^2 radial cut

/// set up free fit parameters with best guess
static const int nPar = 3;
TString paramNames[3] = {"A", "b", "N"};
double paramInits[3] = {-0.12, 0, 0.02};

/// path to experiment data files
TString data_dir = "/media/hickerson/boson/Data/OctetAsym_Offic_2010_FINAL/"; 

/// path to Monte Carlo files
TString mc_dir = "/home/xuansun/Documents/SimProcessedFiles/1mill_beta/";

/// path to save output plots
TString plots_dir = "/home/hickerson/Dropbox/Root/";

/// path to save output root structures 
TString root_output_dir = "/home/hickerson/Documents/";



void output_data_file(TString name, TH1D* h, double ax, double ay)
{
	TString filename = plots_dir + name + ".dat";
    ofstream ofs;
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





/// This needs to be static and global for MINUIT to work
UCNAFierzFitter ucna(KEbins, KEmin, KEmax);
void combined_chi2(Int_t & n, Double_t * /*grad*/ , Double_t &chi2, Double_t *p, Int_t /*iflag */  )
{
    assert(n==3);
    double A=p[0], b=p[1], N=p[2]; // TODO make nPar correct here
	chi2 = ucna.combined_chi2(A,b,N);
}

#if 0 
/** 
 *load the files that contain data histograms
 * XXX MOVED TO UCNAhistogram::fill
 */
int fill_data(TString filename, TString title, 
              TString name, TH1D* histogram)
{
	TFile* tfile = new TFile(data_dir + filename);
	if (tfile->IsZombie()) {
		cout<<"Error loading "<<title<<":\n";
		cout<<"File not found: "<<filename<<".\n";
		return 0;
	}

    if (histogram) {
		cout<<"Warning: Histogram "<<title<<" already exists and is being ignored.\n";
        //cout<<histogram<<endl;
        //delete histogram;
    }

    histogram = (TH1D*)tfile->Get(name);
    if (not histogram) {
		cout<<"Error: In file "<<filename<<":\n";
		cout<<"       Error getting "<<title<<".\n";
		cout<<"       Cannot find histogram named "<<name<<".\n";
        return 0;
    }

	int entries=histogram->GetEntries();
	cout<<"Number of entries in "<<title<<" is "<<entries<<".\n";
    return entries;
}


/**
 * Save root data
 * XXX MOVED TO UCNAhistogram::save() and UCNAmodel::save()
 */
void save_root_data(TString filename, TString title, TString name, 
                   TH1D* rates[2][2], TH1D* super_sum, TH1D* asymmetry)
{
    TString filepath = root_output_dir + filename;
	TFile* tfile = new TFile(filepath, "recreate");
	if (tfile->IsZombie()) {
		cout<<"Error: Problem saving "<<title<<":\n";
		cout<<"       Cannot create file "<<filename<<".\n";
		cout<<"       in path "<<mc_dir<<".\n";
        cout<<"Aborting...\n";
		exit(1);
	}

    if (test_construction(rates, super_sum)) {
        for (int side=0; side<2; side++)
            for (int spin=0; spin<2; spin++) {
                rates[side][spin]->SetDirectory(tfile);
                rates[side][spin]->Write();
            }
    } else {
        cout<<"Error: Rates or super sum for "<<name<<":\n";
        cout<<"       Histogram not constructed.\n";
        cout<<"Aborting...\n";
        exit(1);
    }

    if (super_sum) {
        super_sum->SetDirectory(tfile);
        super_sum->Write();
    } else {
        cout<<"Error: Super sum for "<<name<<":\n";
        cout<<"       Histogram not constructed.\n";
        cout<<"Aborting...\n";
        exit(1);
    }

    if (asymmetry) {
        asymmetry->SetDirectory(tfile);
        asymmetry->Write();
    } else {
        cout<<"Error: Asymmetry for "<<name<<":\n";
        cout<<"       Histogram not constructed.\n";
        cout<<"Aborting...\n";
        exit(1);
    }

    if (ntuple) {
        tntuple->SetDirectory(mc_tfile);
        tntuple->Write();
    } else {
        cout<<"Warning: Ntuple not set. Can't save data.\n";
    }
	tfile->Close();
}
#endif
/**
 * Fill in asymmetry and super_ratio, and super sums from simulation data.
 * Use wild card * in filename where data is split up over many files
 * and they get Tchained together.
 * XXX MOVED TO UCNAhistogram::save() and UCNAmodel::save()
 */
#if 0
int fill_simulation(TString filename, TString title, TString name, 
                    TH1D* rates[2][2], TH1D* super_sum, TH1D* asymmetry)
{
    if (not test_construction(rates, super_sum)) {
        cout<<"Error: Rates or super sum for "<<name<<":\n";
        cout<<"       Histogram not constructed.\n";
        cout<<"Aborting...\n";
        exit(1);
    }

    if (not asymmetry) {
        cout<<"Error: Asymmetry for "<<name<<":\n";
        cout<<"       Histogram not constructed.\n";
        cout<<"Aborting...\n";
        exit(1);
    }

    TString filepath = mc_dir + filename;
	TFile* tfile = new TFile(filepath);
	if (tfile->IsZombie()) {
		cout<<"Error: Problem filling "<<title<<":\n";
		cout<<"       File "<<filename<<" not found\n";
		cout<<"       in path "<<mc_dir<<".\n";
        cout<<"Aborting...\n";
		exit(1);
	}

    TChain *chain = (TChain*)tfile->Get(name);
    if (not chain) {
		cout<<"Error: In file "<<filename<<":\n";
		cout<<"       Cannot get "<<title<<":\n";
		cout<<"       Cannot find chain or tree named "<<name<<".\n";
        cout<<"Aborting...\n";
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

	//TNtuple* tntuple = new TNtuple("mc_ntuple", "MC NTuple", "s:load:energy"); */
	TNtuple* tntuple = 0;

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
            //cout<<"energy: "<<energy<<" side: "<<"side: "<<side
            //         <<" spin: "<<spin<<" afp: "<<afp<<" p: "<<p<<".\n";
            rates[side][spin]->Fill(energy, 1);
            if (tntuple)
			    tntuple->Fill(side, spin, energy);
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
     
	cout<<"Total number of Monte Carlo entries:\n";
	cout<<"      Entries without cuts:  "<<nEvents<<endl;
	cout<<"      Entries with cuts:     "<<nSimmed<<endl;
	cout<<"      Efficiencies of cuts:  "<<(100.0*nSimmed)/double(nEvents)<<"%\n";

	/// compute and normalize super sum
    compute_super_sum(rates, super_sum);
    compute_asymmetry(rates, asymmetry);

    //super_sum.normalize(min_E, max_E);
    //normalize(asymmetry, min_E, max_E);
    //for (int side=0; side<2; side++)
    //    for (int spin=0; spin<2; spin++)
    //        normalize(rates[side][spin], min_E, max_E);

    return nSimmed;
}
#endif

#if 0
// TODO move to UCNAmodel
/// compute little b factor using the Fierz ratio method
TH1D* compute_fierz_ratio(TH1D* data_histogram, TH1D* sm_histogram) {
    TH1D *fierz_ratio_histogram = new TH1D(*data_histogram);
	fierz_ratio_histogram->SetName("fierz_ratio_histogram");
    //fierz_ratio_histogram->Divide(ucna.data.super_sum, ucna.sm.super_sum);
    int bins = data_histogram->GetNbinsX();
    for (int bin = 1; bin < bins+1; bin++) {
		double X = data_histogram->GetBinContent(bin);
		double Y = sm_histogram->GetBinContent(bin);
		double Z = Y>0? X/Y : 0;

		fierz_ratio_histogram->SetBinContent(bin, Z);

		double x = data_histogram->GetBinError(bin);
		double y = sm_histogram->GetBinError(bin);
		fierz_ratio_histogram->SetBinError(bin, Z*Sqrt(x*x/X/X + y*y/Y/Y));
	}
    fierz_ratio_histogram->GetYaxis()->SetRangeUser(0.9,1.1); // Set the range
    fierz_ratio_histogram->SetTitle("Ratio of UCNA data to Monte Carlo");
	cout<<data_histogram->GetNbinsX()<<endl;
	cout<<sm_histogram->GetNbinsX()<<endl;

    /// fit the Fierz ratio 
	char fit_str[1024];
	double expected_fierz = evaluate_expected_fierz(0,0,KEmin,KEmax,18112); /// TODO I'm not sure the 0,0 is right!
    sprintf(fit_str, "1+[0]*(%f/(%f+x)-%f)", m_e, m_e, expected_fierz);
    TF1 *fierz_fit = new TF1("fierz_fit", fit_str, KEmin_b, KEmax_b);
    fierz_fit->SetParameter(0,0);
	fierz_ratio_histogram->Fit(fierz_fit, "Sr");

	/// A fit histogram for output to gnuplot
    TH1D *fierz_fit_histogram = new TH1D(ucna.data.super_sum);
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
    TString fierz_ratio_pdf_filename = plots_dir + "fierz_ratio.pdf";
    TCanvas *canvas = new TCanvas("fierz_canvas", "Fierz component of energy spectrum");
    if (not canvas) {
        cout<<"Can't open new canvas.\n";
        exit(0);
    }
    canvas->SaveAs(fierz_ratio_pdf_filename);

	/// output for gnuplot
	output_data_file("super-sum-data", data_histogram, 1, 1000);
	output_data_file("super-sum-mc", sm_histogram, 1, 1000);
	output_data_file("fierz-ratio", fierz_ratio_histogram, 1, 1);
	output_data_file("fierz-fit", fierz_fit_histogram, 1, 1);

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
    
#endif

/// DISPLAYING AND OUTPUTTING
void draw_histogram(TH1D* histogram, TString name, TString title,
                    TCanvas* canvas = 0, /*TLegend* legend,*/ TString draw = "", int color = 0, int marker = 0)
{
    if (not canvas)
        canvas = new TCanvas(name + "_canvas", title + " Canvas");
    if (not canvas) {
        cout<<"Error: Can't construct a new canvas for "<<name<<".\n";
        exit(1);
    }

    /// Draw a histogram.
	histogram->SetStats(0);
    histogram->SetLineColor(color);
    histogram->SetMarkerStyle(marker);
    histogram->Draw(draw);

    /// Make a pretty legend.
    /*
    if (not legend) {
        TLegend * legend = new TLegend(0.6,0.8,0.7,0.6);
        //legend->AddEntry(&ucna.data.super_sum, "Type 0 super sum", "l");
        //legend->AddEntry(&ucna.sm.super_sum, "Monte Carlo super sum", "p");
        legend->SetTextSize(0.03);
        legend->SetBorderSize(0);
    }
    legend->AddEntry(histogram, title, "p");
    legend->Draw();
    */

    /// Save the data and Mote Carlo plots.
    TString filename = plots_dir + name + ".pdf";
    canvas->SaveAs(filename);
}



/// MAIN APPLICATION
int main(int argc, char *argv[])
{
	TApplication app("Extract Combined A + b Fitter", &argc, argv);
	srand( time(NULL) );    /// set this to make random or repeatable

    /// LOAD 2010 UCNA DATA

    #if 0
    /// Load the files that already contain data asymmetry histogram.
    fill_data("Range_0-1000/CorrectAsym/CorrectedAsym.root",
              "2010 final official asymmetry",
              "hAsym_Corrected_C",
              &ucna.data.asymmetry);

    /// Load the files that already contain data super histogram.
    fill_data("OctetAsym_Offic.root",
              "2010 final official supersum",
              "Total_Events_SuperSum",
              &ucna.data.super_sum);

    /// Load the files that already contain data super histogram.
    for (int side=EAST; side<=WEST; side++)
        for (int afp=EAST; afp<=WEST; afp++) {
            TString s = side? "W":"E", a = afp? "on":"off";
            TString title = "2010 final official "+s+" afp "+a;
            TString cut = "hTotalEvents_"+s+"_"+a+";1";
            int entries = fill_data("OctetAsym_Offic.root", 
                                    title, cut, ucna.data.counts[side][afp]);
            if (entries) {
                cout<<"Status: Number of entries in "<<(side? "west":"east")
                    <<" side with afp "<<a<<" is "<<entries<<".\n";
            }
            else
                cout<<"Error: found no events in "<<title<<".\n";
                /// TODO figure out where these went.
        }
    #endif
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

    /// Load the files that already contain data super histogram.
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
                /// TODO figure out where these went.
        }


    /// LOAD MONTE CARLO SIMULATION EVENTS

    #if 0
    /// Load Monte Carlo simulated Standard Model events
    fill_simulation("SimAnalyzed_Beta_0.root",
                    "Monte Carlo Standard Model beta spectrum",
                    "SimAnalyzed",
                    ucna.sm.counts,
					&ucna.sm.super_sum, 
                    &ucna.sm.asymmetry);

    /// Load Monte Carlo simulated Fierz events
    fill_simulation("SimAnalyzed_Beta_fierz_0.root",
                    "Monte Carlo Fierz beta spectrum",
                    "SimAnalyzed",
                    ucna.fierz.counts,
					&ucna.fierz.super_sum, 
                    &ucna.fierz.asymmetry);
    #endif

    /// Load Monte Carlo simulated Standard Model events
    ucna.sm.fill(mc_dir+"SimAnalyzed_Beta_0.root",
                 "SimAnalyzed",
                 "Monte Carlo Standard Model beta spectrum");

    /// Load Monte Carlo simulated Fierz events
    ucna.fierz.fill(mc_dir+"SimAnalyzed_Beta_fierz_0.root",
                    "SimAnalyzed",
                    "Monte Carlo Fierz beta spectrum");

    /// SAVE ALL HISTOGRAMS 

    // TODO ???


    /// FITTING 

    /// Set up constants and vars
    TMatrixD cov(nPar,nPar);
	double entries = ucna.data.super_sum.GetEffectiveEntries();
	//double A = -0.12;
	//double N = GetEntries(&ucna.data.super_sum, KEmin, KEmax);
    //double N = ucna.data.super_sum.GetEntries();

    /// Set up the fit parameters.
    TF1 asymmetry_func("asymmetry_fit_func", &asymmetry_fit_func, KEmin_A, KEmax_A, nPar);
    //TF1 supersum_func("asymmetry_fit_func", &supersum_fit_func, KEmin_b, KEmax_b, nPar);
    
    /// Setup the parameters.
    asymmetry_func.SetParameters(paramInits);
    //super_sum_func.SetParameters(paramInits);
    for (int i=0; i<nPar; i++) {
        asymmetry_func.SetParName(i, paramNames[i]);
        //super_sum_func.SetParName(i, paramNames[i]);
    }

	/// Actually do the combined fitting.
	//ucna.combined_fit(cov, &asymmetry_func, &supersum_func, &combined_chi2);
	ucna.combined_fit(cov, &asymmetry_func, &combined_chi2);
    double A = asymmetry_func.GetParameter(0);
    double b = asymmetry_func.GetParameter(1);
    double N = asymmetry_func.GetParameter(2);
	ucna.compute_fit(b,N);

	/// Set all expectation values for this range.
    double nSpec = 4;
    TMatrixD expected(nSpec,nSpec);
	for (int m=0; m<nSpec; m++)
		for (int n=0; n<nSpec; n++)
			expected[m][n] = evaluate_expected_fierz(m,n,KEmin_b,KEmax_b,5112);
	
	/// Calculate the predicted inverse covariance matrix for this range.
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

	/// Compute the predicted covariance matrix by inverse.
	double det = 0;
	TMatrixD p_cov = p_cov_inv.Invert(&det);



    /// PRINT OUT REPORT OF FITS, CORRELATIONS AND ERRORS

	/// Output the data info.
    int cl = 14;
    cout<<setprecision(5);
	cout<<" ENERGY RANGE:\n";
	cout<<"    Full Energy range is "<<KEmin<<" - "<<KEmax<<" keV.\n";
	cout<<"    Energy range for A fit is "<<KEmin_A<<" - "<<KEmax_A<<" keV.\n";
	cout<<"    Energy range for b fit is "<<KEmin_b<<" - "<<KEmax_b<<" keV.\n";
	cout<<"    Number of counts in full data is "<<(int)entries<<".\n";
	cout<<"    Number of counts in energy range is "<<(int)N<<".\n";
	cout<<"    Efficiency energy cut is "<< N/entries*100<<"%.\n";

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
    TCanvas *canvas = new TCanvas("fierz_fitter_canvas",
                                  "Combined Fierz component of energy spectrum");
    if (not canvas) {
        cout<<"Error: Can't open new canvas.\n";
        exit(1);
    }

/*
	ucna.fierz.super_sum.SetStats(0);
    ucna.fierz.super_sum.SetLineColor(3);
    ucna.fierz.super_sum.Draw("");

    TString fit_pdf_filename = plots_dir + "combined_fierz_fit_data.pdf";
    ucna.sm.super_sum.SetLineColor(1);
    ucna.sm.super_sum.Draw("Same");
    canvas->SaveAs(fit_pdf_filename);
    */

    /*
    draw_histogram(&ucna.sm.asymmetry, 
                   "sm_asymmetry", 
                   "Combined Asymmetry", 
                   canvas, "", 3, 0);
    */
    draw_histogram(&ucna.data.asymmetry, 
                   "data_asymmetry", 
                   "UCNA 2010 Asymmetry", 
                   canvas,"",1,0);

    TString asymmetry_pdf_filename = plots_dir + "data_asymmetry.pdf";
    canvas->SaveAs(asymmetry_pdf_filename);

    draw_histogram(&ucna.sm.super_sum, 
                   "sm_supersum", 
                   "Standard Model Monte Carlo super sum", 
                   canvas,"",4,0);
    draw_histogram(&ucna.fierz.super_sum, 
                   "fierz_supersum", 
                   "Fierz Monte Carlo super sum", 
                   canvas,"Same",1,0);

    TString monte_carlo_pdf_filename = plots_dir + "monte_carlo.pdf";
    canvas->SaveAs(monte_carlo_pdf_filename);

    draw_histogram(&ucna.data.super_sum, 
                   "data_supersum", 
                   "Data super sum", 
                   canvas,"",2,0);

    draw_histogram(&ucna.fit.super_sum, 
                   "fit_supersum", 
                   "Fit super sum", 
                   canvas,"Same",6,0);

    TString data_supersum_pdf_filename = plots_dir + "data_supersum.pdf";
    canvas->SaveAs(data_supersum_pdf_filename);

    /*
    /// Draw the data super sums
    ucna.data.super_sum.Scale(200);
	ucna.data.super_sum.SetStats(0);
    ucna.data.super_sum.SetLineColor(2);
    ucna.data.super_sum.Draw("");

    /// Draw Monte Carlo super sum
	ucna.data.super_sum.SetStats(0);
    ucna.sm.super_sum.SetLineColor(1);
    ucna.sm.super_sum.SetMarkerStyle(4);
    ucna.sm.super_sum.Draw("same p0");

    /// Make a pretty legend.
    TLegend * legend = new TLegend(0.6,0.8,0.7,0.6);
    legend->AddEntry(ucna.data.super_sum, "Type 0 super sum", "l");
    legend->AddEntry(ucna.sm.super_sum, "Monte Carlo super sum", "p");
    legend->SetTextSize(0.03);
    legend->SetBorderSize(0);
    legend->Draw();

    /// Save the data and Mote Carlo plots.
    TString super_sum_pdf_filename = plots_dir + "super_sum_data.pdf";
    canvas->SaveAs(super_sum_pdf_filename);
*/

	/*
	/// A fit histogram for output to gnuplot
    TH1D *fierz_fit_histogram = new TH1D(*asymmetry_histogram);
	for (int i = 0; i < fierz_fit_histogram->GetNbinsX(); i++)
		fierz_fit_histogram->SetBinContent(i, fierz_fit->Eval(fierz_fit_histogram->GetBinCenter(i)));
	

	/// output for root
    TString pdf_filename = "/data/kevinh/mc/asymmetry_fierz_term_fit.pdf";
    canvas->SaveAs(pdf_filename);

	/// output for gnuplot
	//output_data_file("/data/kevinh/mc/super-sum-data", ucna.data.super_sum, 1, 1000);
	//output_data_file("/data/kevinh/mc/super-sum-mc", ucna.sm.super_sum, 1, 1000);
	//output_data_file("/data/kevinh/mc/fierz-ratio", fierz_ratio_histogram, 1, 1);
	//output_data_file("/data/kevinh/mc/fierz-fit", fierz_fit_histogram, 1, 1);

	*/

#if 0
	// Create a new canvas.
	TCanvas * c1 = new TCanvas("c1","Two Histogram Fit example",100,10,900,800);
	c1->Divide(2,2);
	gStyle->SetOptFit();
	gStyle->SetStatY(0.6);

	c1->cd(1);
	ucna.data.asymmetry.Draw();
	func.SetRange(min_E, max_E);
	func.DrawCopy("cont1 same");
	/*
	c1->cd(2);
	asymmetry->Draw("lego");
	func.DrawCopy("surf1 same");
	*/
	c1->cd(3);
	func.SetRange(min_E, max_E);
	///TODO fierzratio_histogram->Draw();
	func.DrawCopy("cont1 same");
	/*
	c1->cd(4);
	fierzratio->Draw("lego");
	gPad->SetLogz();
	func.Draw("surf1 same");
	*/
#endif
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
