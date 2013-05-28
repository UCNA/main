// UCNA includes
#include "G4toPMT.hh"
#include "PenelopeToPMT.hh"
#include "CalDBSQL.hh"

// ROOT includes
#include <TH1F.h>
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
const double    L           = -1.27590;             /// gA/gV               (450)
const double    M           = 1 + 3*L*L;            /// matrix element
const double    I_0         = 1.63632;              /// 0th moment          (25)
const double    I_1         = 1.07017;              /// 1st moment          (15)
const double    x_1         = I_1/I_0;              /// first m/E moment    (9)


static double electron_mass = 510.9989; 	// needed for the physics of Fierz interference
double min_E = 230;
double max_E = 660;
static double expected_fierz = 0.6540;	// full range
static int integral_size = 1000;
//static double expected_fierz = 0.6111;		// for range 150 - 600
//static double expected_gluck = 11.8498;     // for range 150 - 600
//static unsigned nToSim = 5E7;				// how many triggering events to simulate
//static double loading_prob = 40; 		// ucn loading probability (percent)
//static int bins = 150;						// replace with value from data or smoothing fit
//double scale_x = 1.015;
//double scale_x = 1.0;

#if 1
using namespace std;
int output_histogram(string filename, TH1F* h, double ax, double ay)
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



void compute_fit(TH1F* histogram, TF1* fierz_fit) 
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
	double E = x[0] + electron_mass;

	return A / (1 + b * electron_mass / E);
}



double fierzratio_fit_func(double *x, double *par)
{
	double b = par[1];
	double E = x[0] + electron_mass;

	return (1 + b * electron_mass / E) / (1 + b * expected_fierz);
}



// data need to be globals to be visible by fcn 

vector<double> asymmetry_energy;        
vector<double> asymmetry_values;        
vector<double> asymmetry_errors;        
vector<double> fierzratio_energy;        
vector<double> fierzratio_values;        
vector<double> fierzratio_errors;        



void combined_chi2(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *p, Int_t /*iflag */  )
{
	double chi2 = 0; 
	double chi,	E; 

	int n = asymmetry_energy.size();
	for (int i = 0; i < n; ++i )
	{
		E = asymmetry_energy[i];
		chi = (asymmetry_values[i] - asymmetry_fit_func(&E,p)) / asymmetry_errors[i];
		chi2 += chi*chi; 
	}

	n = fierzratio_energy.size();
	for (int i = 0; i < n; ++i ) { 
		E = fierzratio_energy[i];
		chi = (fierzratio_values[i] - fierzratio_fit_func(&E,p)) / fierzratio_errors[i];
		chi2 += chi*chi; 
	}
	fval = chi2; 
}



TF1* combined_fit(TH1F* asymmetry, TH1F* fierzratio) 
{ 
	double iniParams[2] = { -0.15, 0 };
	const char * iniParamNames[2] = { "A", "b" };
	// create fit function
	TF1 * func = new TF1("func", asymmetry_fit_func, min_E, max_E, 2);
	func->SetParameters(iniParams);
	for (int i = 0; i < 2; i++)
		func->SetParName(i, iniParamNames[i]);

	//if (true) { 
		// fill data structure for fit (coordinates + values + errors) 
		std::cout << "Do global fit" << std::endl;
		// fit now all the function together

		// fill data structure for fit (coordinates + values + errors) 
		TAxis *xaxis1  = asymmetry->GetXaxis();
		TAxis *xaxis2  = fierzratio->GetXaxis();

		int nbinX1 = asymmetry->GetNbinsX(); 
		int nbinX2 = fierzratio->GetNbinsX(); 

		/// reset data structure
		asymmetry_energy = vector<double>();
		asymmetry_values = vector<double>();
		asymmetry_errors = vector<double>();
		fierzratio_energy = vector<double>();
		fierzratio_values = vector<double>();
		fierzratio_errors = vector<double>();

		for (int ix = 1; ix <= nbinX1; ++ix)
		{
			double E = xaxis1->GetBinCenter(ix);
			if (min_E < E and E < max_E)
			{
				asymmetry_energy.push_back( E );
				asymmetry_values.push_back( asymmetry->GetBinContent(ix) );
				asymmetry_errors.push_back( asymmetry->GetBinError(ix) );
				//cout << xaxis1->GetBinCenter(ix) << endl;
			}
		}

		for (int ix = 1; ix <= nbinX2; ++ix)
		{
			double E = xaxis2->GetBinCenter(ix);
			if (min_E < E and E < max_E)
			{
				fierzratio_energy.push_back( E );
				fierzratio_values.push_back( fierzratio->GetBinContent(ix) );
				fierzratio_errors.push_back( fierzratio->GetBinError(ix) );
			}
		}


		TVirtualFitter::SetDefaultFitter("Minuit");
		TVirtualFitter * minuit = TVirtualFitter::Fitter(0,2);
		for (int i = 0; i < 2; ++i) {  
			minuit->SetParameter(i, func->GetParName(i), func->GetParameter(i), 1, 0, 0);
		}
		minuit->SetFCN(combined_chi2);
		minuit->SetErrorDef(1);	// 1 for chi^2

		double arglist[100];
		arglist[0] = 0;
		// set print level
		minuit->ExecuteCommand("SET PRINT",arglist,1);

		// minimize
		arglist[0] = 50; // number of function calls
		arglist[1] = 0.1; // tolerance
		minuit->ExecuteCommand("MIGRAD",arglist,2);

		//get result
		double minParams[2];
		double parErrors[2];
		for (int i = 0; i < 2; ++i) {  
			minParams[i] = minuit->GetParameter(i);
			parErrors[i] = minuit->GetParError(i);
		}
		double chi2, edm, errdef; 
		int nvpar, nparx;
		minuit->GetStats(chi2,edm,errdef,nvpar,nparx);

		func->SetParameters(minParams);
		func->SetParErrors(parErrors);
		func->SetChisquare(chi2);
		int ndf = asymmetry_energy.size() + fierzratio_energy.size()- nvpar;
		func->SetNDF(ndf);

		cout << "chi^2 = " << chi2 << ", ndf = " << ndf << ", chi^2/ndf = " << chi2/ndf << endl;


		// add to list of functions
		asymmetry->GetListOfFunctions()->Add(func);
		fierzratio->GetListOfFunctions()->Add(func);
		
	/*
	}
	else {     
		// fit independently
		asymmetry->Fit(func);
		fierzratio->Fit(func);
	}
	*/


	return func; 
}




int main(int argc, char *argv[]) {
	TApplication app("Combined Fit", &argc, argv);
	expected_fierz = evaluate_expected_fierz(min_E, max_E);
    //TCanvas *canvas = new TCanvas("fierz_canvas", "Fierz component of energy spectrum");

	//asym_histogram->SetStats(0);
    //asym_histogram->SetLineColor(3);
    //asym_histogram->Draw("");

    TFile *asymmetry_data_tfile = new TFile(
		"/home/mmendenhall/Plots/OctetAsym_Offic/Range_0-1000/CorrectAsym/CorrectedAsym.root");
	if (asymmetry_data_tfile->IsZombie())
	{
		std::cout << "File not found." << std::endl;
		exit(1);
	}

    TFile *ucna_data_tfile = new TFile(
		"/home/mmendenhall/Plots/OctetAsym_Offic/OctetAsym_Offic.root");
	if (ucna_data_tfile->IsZombie())
	{
		std::cout << "File not found." << std::endl;
		exit(1);
	}

    TFile *fierzratio_data_tfile = new TFile(
		"Fierz/ratio.root");
	if (fierzratio_data_tfile->IsZombie())
	{
		std::cout << "File not found." << std::endl;
		exit(1);
	}

    TH1F *asymmetry_histogram = 
            (TH1F*)asymmetry_data_tfile->Get("hAsym_Corrected_C");
    TH1F *supersum_histogram = 
            (TH1F*)ucna_data_tfile->Get("Total_Events_SuperSum");
    TH1F *fierzratio_histogram = 
            (TH1F*)fierzratio_data_tfile->Get("fierz_ratio_histogram");

	/*
	// fit the Fierz ratio 
	char fit_str[1024];
    sprintf(fit_str, "[0]/(1+[1]*(%f/(%f+x)))", electron_mass, electron_mass);
    TF1 *fierz_fit = new TF1("fierz_fit", fit_str, min_E, max_E);
    fierz_fit->SetParameter(0,-0.12);
    fierz_fit->SetParameter(1,0);
	asymmetry_histogram->Fit(fierz_fit, "Sr");

	compute_fit(asymmetry_histogram, fierz_fit);
	*/

	TF1* func = combined_fit(asymmetry_histogram, fierzratio_histogram);

	double entries = supersum_histogram->GetEffectiveEntries();
	double part_int = supersum_histogram->Integral(
					  supersum_histogram->FindBin(min_E),
					  supersum_histogram->FindBin(max_E));
	double full_int = supersum_histogram->Integral();
	double N = entries * part_int / full_int;
	cout << "The energy range is " << min_E << " - " << max_E << " keV" << endl;
	cout << "Number of counts in full data is " << (int)entries << endl;
	cout << "Number of counts in energy range is " <<  (int)N << endl;
	cout << "The expected statistical error for A in this range is " << 2.7 / sqrt(N) << endl;
	cout << "The actual statistical error for A in this range is " << func->GetParError(0) << endl;
	cout << "The ratio for A error is " << func->GetParError(0) * sqrt(N) / 2.7 << endl;
	cout << "The expected statistical error for b in this range is " << 14.8 / sqrt(N) << endl;
	cout << "The actual statistical error for b in this range is " << func->GetParError(1) << endl;
	cout << "The ratio for b error is " << func->GetParError(1) * sqrt(N) / 14.8 << endl;

	/*
	// A fit histogram for output to gnuplot
    TH1F *fierz_fit_histogram = new TH1F(*asymmetry_histogram);
	for (int i = 0; i < fierz_fit_histogram->GetNbinsX(); i++)
		fierz_fit_histogram->SetBinContent(i, fierz_fit->Eval(fierz_fit_histogram->GetBinCenter(i)));

	

	// output for root
    TString pdf_filename = "/data/kevinh/mc/asymmetry_fierz_term_fit.pdf";
    canvas->SaveAs(pdf_filename);

	// output for gnuplot
	//output_histogram("/data/kevinh/mc/super-sum-data.dat", super_sum_histogram, 1, 1000);
	//output_histogram("/data/kevinh/mc/super-sum-mc.dat", mc.sm_super_sum_histogram, 1, 1000);
	//output_histogram("/data/kevinh/mc/fierz-ratio.dat", fierz_ratio_histogram, 1, 1);
	//output_histogram("/data/kevinh/mc/fierz-fit.dat", fierz_fit_histogram, 1, 1);

	*/

	// Create a new canvas.
	TCanvas * c1 = new TCanvas("c1","Two HIstogram Fit example",100,10,900,800);
	c1->Divide(2,2);
	gStyle->SetOptFit();
	gStyle->SetStatY(0.6);

	c1->cd(1);
	asymmetry_histogram->Draw();
	func->SetRange(min_E, max_E);
	func->DrawCopy("cont1 same");
	/*
	c1->cd(2);
	asymmetry->Draw("lego");
	func->DrawCopy("surf1 same");
	*/
	c1->cd(3);
	func->SetRange(min_E, max_E);
	fierzratio_histogram->Draw();
	func->DrawCopy("cont1 same");
	/*
	c1->cd(4);
	fierzratio->Draw("lego");
	gPad->SetLogz();
	func->Draw("surf1 same");
	*/

	app.Run();

	return 0;
}
