// UCNA includes
#include "G4toPMT.hh"
#include "PenelopeToPMT.hh"
#include "CalDBSQL.hh"
#include "FierzFitter.hh"

// ROOT includes
#include <TH1F.h>
#include <TLegend.h>
#include <TF1.h>
#include <TVirtualFitter.h>
#include <TList.h>
#include <TStyle.h>
#include <TApplication.h>
#include <TMatrixD.h>

// C++ includes
#include <iostream>
#include <fstream>
#include <string>

// C includes
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

double min_E = 220;
double max_E = 670;
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


// data need to be globals to be visible by fcn 
vector<double> asymmetry_energy;        
vector<double> asymmetry_values;        
vector<double> asymmetry_errors;        
vector<double> fierzratio_energy;        
vector<double> fierzratio_values;        
vector<double> fierzratio_errors;        

double expected_fierz_0_1;
double expected_fierz_0_2;
double expected_fierz_2_0;
double expected_fierz_2_1;


void combined_chi2(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *p, Int_t /*iflag */  )
{
	double chi2 = 0; 
	double chi,	E; 

	int n = asymmetry_energy.size();
	for (int i = 0; i < n; ++i )
	{
		double par[2] = {p[0],p[1]};
		E = asymmetry_energy[i];
		chi = (asymmetry_values[i] - asymmetry_fit_func(&E,par)) / asymmetry_errors[i];
		chi2 += chi*chi; 
	}

	n = fierzratio_energy.size();
	for (int i = 0; i < n; ++i ) { 
		double par[2] = {p[1], expected_fierz_0_1};
		E = fierzratio_energy[i];
		chi = (fierzratio_values[i] - fierzratio_fit_func(&E,par)) / fierzratio_errors[i];
		chi2 += chi*chi; 
	}
	fval = chi2; 
}



TF1* combined_fit(TH1F* asymmetry, TH1F* fierzratio, double cov[2][2]) 
{ 
	int nPar = 2;
	double iniParams[2] = { -0.15, 0 };
	const char * iniParamNames[2] = { "A", "b" };
	// create fit function
	TF1 * func = new TF1("func", asymmetry_fit_func, min_E, max_E, nPar);
	func->SetParameters(iniParams);
	for (int i = 0; i < nPar; i++)
		func->SetParName(i, iniParamNames[i]);

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
	TVirtualFitter * minuit = TVirtualFitter::Fitter(0,nPar);
	for (int i = 0; i < nPar; ++i) {  
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
	minuit->ExecuteCommand("MIGRAD",arglist,nPar);

	//get result
	double minParams[nPar];
	double parErrors[nPar];
	for (int i = 0; i < nPar; ++i) {  
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

	TMatrixD matrix( nPar, nPar, minuit->GetCovarianceMatrix() );
	for (int i = 0; i < nPar; i++)
		for (int j = 0; j < nPar; j++)
			cov[i][j] = minuit->GetCovarianceMatrixElement(i,j);

	cout << "chi^2 = " << chi2 << ", ndf = " << ndf << ", chi^2/ndf = " << chi2/ndf << endl;

	return func; 
}




int main(int argc, char *argv[]) {
	TApplication app("Combined Fit", &argc, argv);
	expected_fierz_0_1 = evaluate_expected_fierz(0, 1, min_E, max_E);
	expected_fierz_0_2 = evaluate_expected_fierz(0, 2, min_E, max_E);
	expected_fierz_2_0 = evaluate_expected_fierz(2, 0, min_E, max_E);
	expected_fierz_2_1 = evaluate_expected_fierz(2, 1, min_E, max_E);

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

	double cov[2][2]; 
	double entries = supersum_histogram->GetEffectiveEntries();
	double N = GetEntries(supersum_histogram, min_E, max_E);

	TF1* func = combined_fit(asymmetry_histogram, fierzratio_histogram, cov);

	cout << " COVARIANCE MATRIX cov(A,b) =\n";
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			cout << "\t\t" << cov[i][j];
		}
		cout << "\n";
	}

	cout << " The energy range is " << min_E << " - " << max_E << " keV" << endl;
	cout << " Number of counts in full data is " << (int)entries << endl;
	cout << " Number of counts in energy range is " <<  (int)N << endl;
	cout << " The expected statistical error for A in this range is " << 2.7 / sqrt(N) << endl;
	cout << " The actual statistical error for A in this range is " << func->GetParError(0) << endl;
	cout << " The ratio for A error is " << func->GetParError(0) * sqrt(N) / 2.7 << endl;
	cout << " The expected statistical error for b in this range is " << 14.8 / sqrt(N) << endl;
	cout << " The actual statistical error for b in this range is " << func->GetParError(1) << endl;
	cout << " The ratio for b error is " << func->GetParError(1) * sqrt(N) / 14.8 << endl;
	cout << " cor(A,b) = " << cov[1][0] / sqrt(cov[0][0] * cov[1][1]) << endl;

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
