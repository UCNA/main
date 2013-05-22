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

// c++ includes
#include <iostream>
#include <fstream>
#include <string>

// c includes
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

static double electron_mass = 510.9989; 	// needed for the physics of Fierz interference
double min_E = 230;
double max_E = 660;
static double expected_fierz = 0.6540;	// full range
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
	double E = x[0];

	return A * (1 + b * electron_mass / E);
}



double fierz_ratio_fit_func(double *x, double *par)
{
	double b = par[1];
	double E = x[0];

	return (1 + b * electron_mass / E) / (1 + b * expected_fierz);
}



// data need to be globals to be visible by fcn 

vector<double> energy;        
vector<double> values;        
vector<double> errors;        

void combined_chi2(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *p, Int_t /*iflag */  )
{
	int n = energy.size();
	double chi2 = 0; 
	double chi,	E; 
	for (int i = 0; i <n; ++i ) { 
		E = energy[i];
		chi = (values[i] - asymmetry_fit_func(&E,p)) / errors[i];
		chi2 += chi*chi; 
		chi = (values[i] - fierz_ratio_fit_func(&E,p)) / errors[i];
		chi2 += chi*chi; 
	}
	fval = chi2; 
}



int combined_fit(TH1F* asymmetry, TH1F* fierz_ratio) 
{ 
	double iniParams[2] = { -0.12, 0 };
	// create fit function
	TF1 * func = new TF1("func", combined_chi2, min_E, max_E, 2);
	func->SetParameters(iniParams);

	//if (true) { 
		// fill data structure for fit (coordinates + values + errors) 
		std::cout << "Do global fit" << std::endl;
		// fit now all the function together

		// fill data structure for fit (coordinates + values + errors) 
		TAxis *xaxis1  = asymmetry->GetXaxis();
		TAxis *xaxis2  = fierz_ratio->GetXaxis();

		int nbinX1 = asymmetry->GetNbinsX(); 
		int nbinX2 = fierz_ratio->GetNbinsX(); 

		/// reset data structure
		energy = vector<double>();
		values = vector<double>();
		errors = vector<double>();

		for (int ix = 1; ix <= nbinX1; ++ix)
			if (asymmetry->GetBinContent(ix) > 0)
			{
				energy.push_back( xaxis1->GetBinCenter(ix) );
				values.push_back( asymmetry->GetBinContent(ix) );
				errors.push_back( asymmetry->GetBinError(ix) );
			}

		for (int ix = 1; ix <= nbinX2; ++ix)
			if (fierz_ratio->GetBinContent(ix) > 0)
			{
				energy.push_back( xaxis2->GetBinCenter(ix) );
				values.push_back( fierz_ratio->GetBinContent(ix) );
				errors.push_back( fierz_ratio->GetBinError(ix) );
			}


		TVirtualFitter::SetDefaultFitter("Minuit");
		TVirtualFitter * minuit = TVirtualFitter::Fitter(0,2);
		for (int i = 0; i < 10; ++i) {  
			minuit->SetParameter(i, func->GetParName(i), func->GetParameter(i), 0.01, 0,0);
		}
		minuit->SetFCN(combined_chi2);

		double arglist[100];
		arglist[0] = 0;
		// set print level
		minuit->ExecuteCommand("SET PRINT",arglist,2);

		// minimize
		arglist[0] = 5000; // number of function calls
		arglist[1] = 0.01; // tolerance
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
		int ndf = energy.size()-nvpar;
		func->SetNDF(ndf);

		cout << "Chi2 Fit = " << chi2 << " ndf = " << ndf << "  " << func->GetNDF() << endl;

		// add to list of functions
		asymmetry->GetListOfFunctions()->Add(func);
		fierz_ratio->GetListOfFunctions()->Add(func);
	/*
	}
	else {     
		// fit independently
		asymmetry->Fit(func);
		fierz_ratio->Fit(func);
	}
	*/



	// Create a new canvas.
	TCanvas * c1 = new TCanvas("c1","Two HIstogram Fit example",100,10,900,800);
	c1->Divide(2,2);
	gStyle->SetOptFit();
	gStyle->SetStatY(0.6);

	c1->cd(1);
	asymmetry->Draw();
	func->SetRange(min_E, max_E);
	func->DrawCopy("cont1 same");
	c1->cd(2);
	asymmetry->Draw("lego");
	func->DrawCopy("surf1 same");
	c1->cd(3);
	func->SetRange(min_E, max_E);
	fierz_ratio->Draw();
	func->DrawCopy("cont1 same");
	c1->cd(4);
	fierz_ratio->Draw("lego");
	gPad->SetLogz();
	func->Draw("surf1 same");

	return 0; 
}




int main(int argc, char *argv[]) {
    TCanvas *canvas = new TCanvas("fierz_canvas", "Fierz component of energy spectrum");

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

    TH1F *asymmetry_histogram = 
            (TH1F*)asymmetry_data_tfile->Get("hAsym_Corrected_C");
    //printf("Number of bins in data %d\n", ucna_data_histogram->GetNbinsX());

    TH1F *supersum_histogram = 
            (TH1F*)ucna_data_tfile->Get("Total_Events_SuperSum");

	// fit the Fierz ratio 
	char fit_str[1024];
    sprintf(fit_str, "[0]/(1+[1]*(%f/(%f+x)))", electron_mass, electron_mass);
    TF1 *fierz_fit = new TF1("fierz_fit", fit_str, min_E, max_E);
    fierz_fit->SetParameter(0,-0.12);
    fierz_fit->SetParameter(1,0);
	asymmetry_histogram->Fit(fierz_fit, "Sr");

	compute_fit(asymmetry_histogram, fierz_fit);

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

	return 0;
}




#if 0
// --------------------------------------------------
//
// Root example
// + Example to fit two histograms at the same time 
// Author: Rene Brun

#include "TH2D.h"
#include "TF2.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TRandom3.h"
#include "TVirtualFitter.h"
#include "TList.h"

#include <vector>
#include <map>
#include <iostream>

double gauss2D(double *x, double *par) {
	double z1 = double((x[0]-par[1])/par[2]);
	double z2 = double((x[1]-par[3])/par[4]);
	return par[0]*exp(-0.5*(z1*z1+z2*z2));
}   

double my2Dfunc(double *x, double *par) {
	double *p1 = &par[0];
	double *p2 = &par[5];
	return gauss2D(x,p1) + gauss2D(x,p2);
}



// data need to be globals to be visible by fcn 

std::vector<std::pair<double, double> > coords;        
std::vector<double > values;        
std::vector<double > errors;        

void myFcn(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *p, Int_t /*iflag */  )
{
	int n = coords.size();
	double chi2 = 0; 
	double tmp,x[2]; 
	for (int i = 0; i <n; ++i ) { 
		x[0] = coords[i].first;
		x[1] = coords[i].second;
		tmp = ( values[i] - my2Dfunc(x,p))/errors[i];
		chi2 += tmp*tmp; 
	}
	fval = chi2; 
}
TRandom3 rndm; 
void FillHisto(TH2D * h, int n, double * p) { 


	const double mx1 = p[1]; 
	const double my1 = p[3]; 
	const double sx1 = p[2]; 
	const double sy1 = p[4]; 
	const double mx2 = p[6]; 
	const double my2 = p[8]; 
	const double sx2 = p[7]; 
	const double sy2 = p[9]; 
	//const double w1 = p[0]*sx1*sy1/(p[5]*sx2*sy2); 
	const double w1 = 0.5; 

	double x, y; 
	for (int i = 0; i < n; ++i) {
		// generate randoms with larger gaussians
		rndm.Rannor(x,y);

		double r = rndm.Rndm(1);
		if (r < w1) { 
			x = x*sx1 + mx1; 
			y = y*sy1 + my1; 
		}
		else { 
			x = x*sx2 + mx2; 
			y = y*sy2 + my2; 
		}      
		h->Fill(x,y);

	}
}




int TwoHistoFit2D(bool global = true) { 

	// create two histograms 

	int nbx1 = 50;
	int nby1 = 50;
	int nbx2 = 50;
	int nby2 = 50;
	double xlow1 = 0.; 
	double ylow1 = 0.; 
	double xup1 = 10.; 
	double yup1 = 10.; 
	double xlow2 = 5.; 
	double ylow2 = 5.; 
	double xup2 = 20.; 
	double yup2 = 20.; 

	TH2D * h1 = new TH2D("h1","core",nbx1,xlow1,xup1,nby1,ylow1,yup1);
	TH2D * h2 = new TH2D("h2","tails",nbx2,xlow2,xup2,nby2,ylow2,yup2);

	double iniParams[10] = { 100, 6., 2., 7., 3, 100, 12., 3., 11., 2. };
	// create fit function
	TF2 * func = new TF2("func",my2Dfunc,xlow2,xup2,ylow2,yup2, 10);
	func->SetParameters(iniParams);

	// fill Histos
	int n1 = 1000000;
	int n2 = 1000000; 
	//  h1->FillRandom("func", n1);
	//h2->FillRandom("func",n2);
	FillHisto(h1,n1,iniParams);
	FillHisto(h2,n2,iniParams);

	// scale histograms to same heights (for fitting)
	double dx1 = (xup1-xlow1)/double(nbx1); 
	double dy1 = (yup1-ylow1)/double(nby1);
	double dx2 = (xup2-xlow2)/double(nbx2);
	double dy2 = (yup2-ylow2)/double(nby2);
	//   h1->Sumw2();
	//   h1->Scale( 1.0 / ( n1 * dx1 * dy1 ) );
	// scale histo 2 to scale of 1 
	h2->Sumw2();
	h2->Scale(  ( double(n1) * dx1 * dy1 )  / ( double(n2) * dx2 * dy2 ) );


	if (global) { 
		// fill data structure for fit (coordinates + values + errors) 
		std::cout << "Do global fit" << std::endl;
		// fit now all the function together

		// fill data structure for fit (coordinates + values + errors) 
		TAxis *xaxis1  = h1->GetXaxis();
		TAxis *yaxis1  = h1->GetYaxis();
		TAxis *xaxis2  = h2->GetXaxis();
		TAxis *yaxis2  = h2->GetYaxis();

		int nbinX1 = h1->GetNbinsX(); 
		int nbinY1 = h1->GetNbinsY(); 
		int nbinX2 = h2->GetNbinsX(); 
		int nbinY2 = h2->GetNbinsY(); 

		/// reset data structure
		coords = std::vector<std::pair<double,double> >();
		values = std::vector<double>();
		errors = std::vector<double>();


		for (int ix = 1; ix <= nbinX1; ++ix) { 
			for (int iy = 1; iy <= nbinY1; ++iy) { 
				if ( h1->GetBinContent(ix,iy) > 0 ) { 
					coords.push_back( std::make_pair(xaxis1->GetBinCenter(ix), yaxis1->GetBinCenter(iy) ) );
					values.push_back( h1->GetBinContent(ix,iy) );
					errors.push_back( h1->GetBinError(ix,iy) );
				}
			}
		}
		for (int ix = 1; ix <= nbinX2; ++ix) { 
			for (int iy = 1; iy <= nbinY2; ++iy) { 
				if ( h2->GetBinContent(ix,iy) > 0 ) { 
					coords.push_back( std::make_pair(xaxis2->GetBinCenter(ix), yaxis2->GetBinCenter(iy) ) );
					values.push_back( h2->GetBinContent(ix,iy) );
					errors.push_back( h2->GetBinError(ix,iy) );
				}
			}
		}

		TVirtualFitter::SetDefaultFitter("Minuit");
		TVirtualFitter * minuit = TVirtualFitter::Fitter(0,10);
		for (int i = 0; i < 10; ++i) {  
			minuit->SetParameter(i, func->GetParName(i), func->GetParameter(i), 0.01, 0,0);
		}
		minuit->SetFCN(myFcn);

		double arglist[100];
		arglist[0] = 0;
		// set print level
		minuit->ExecuteCommand("SET PRINT",arglist,2);

		// minimize
		arglist[0] = 5000; // number of function calls
		arglist[1] = 0.01; // tolerance
		minuit->ExecuteCommand("MIGRAD",arglist,2);

		//get result
		double minParams[10];
		double parErrors[10];
		for (int i = 0; i < 10; ++i) {  
			minParams[i] = minuit->GetParameter(i);
			parErrors[i] = minuit->GetParError(i);
		}
		double chi2, edm, errdef; 
		int nvpar, nparx;
		minuit->GetStats(chi2,edm,errdef,nvpar,nparx);

		func->SetParameters(minParams);
		func->SetParErrors(parErrors);
		func->SetChisquare(chi2);
		int ndf = coords.size()-nvpar;
		func->SetNDF(ndf);

		std::cout << "Chi2 Fit = " << chi2 << " ndf = " << ndf << "  " << func->GetNDF() << std::endl;

		// add to list of functions
		h1->GetListOfFunctions()->Add(func);
		h2->GetListOfFunctions()->Add(func);
	}
	else {     
		// fit independently
		h1->Fit(func);
		h2->Fit(func);
	}



	// Create a new canvas.
	TCanvas * c1 = new TCanvas("c1","Two HIstogram Fit example",100,10,900,800);
	c1->Divide(2,2);
	gStyle->SetOptFit();
	gStyle->SetStatY(0.6);

	c1->cd(1);
	h1->Draw();
	func->SetRange(xlow1,ylow1,xup1,yup1);
	func->DrawCopy("cont1 same");
	c1->cd(2);
	h1->Draw("lego");
	func->DrawCopy("surf1 same");
	c1->cd(3);
	func->SetRange(xlow2,ylow2,xup2,yup2);
	h2->Draw();
	func->DrawCopy("cont1 same");
	c1->cd(4);
	h2->Draw("lego");
	gPad->SetLogz();
	func->Draw("surf1 same");

	return 0; 
}
#endif
