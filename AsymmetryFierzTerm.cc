/// \file DataScannerExample.cc example code for using MC simulation data
#include "G4toPMT.hh"
#include "PenelopeToPMT.hh"
#include "CalDBSQL.hh"
#include <TH1F.h>
#include <TLegend.h>
//#include <TFitResult.h> // v5.27
#include <TF1.h>

#include <iostream>
#include <fstream>
#include <string>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

static double electron_mass = 510.9989; 	// needed for the physics of Fierz interference
double min_E = 220;
double max_E = 670;
//static double expected_fierz = 0.6540;	// full range
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

int main(int argc, char *argv[]) {
    TCanvas *canvas = new TCanvas("fierz_canvas", "Fierz component of energy spectrum");

	//asym_histogram->SetStats(0);
    //asym_histogram->SetLineColor(3);
    //asym_histogram->Draw("");

    TFile *ucna_data_tfile = new TFile(
		"/home/mmendenhall/Plots/OctetAsym_Offic/Range_0-1000/CorrectAsym/CorrectedAsym.root");
	if (ucna_data_tfile->IsZombie())
	{
		std::cout << "File not found." << std::endl;
		exit(1);
	}

    TH1F *asymmetry_histogram = 
            (TH1F*)ucna_data_tfile->Get("hAsym_Corrected_C");
    //printf("Number of bins in data %d\n", ucna_data_histogram->GetNbinsX());


	// fit the Fierz ratio 
	char fit_str[1024];
    sprintf(fit_str, "[0]/(1+[1]*(%f/(%f+x)))", electron_mass, electron_mass);
    TF1 *fierz_fit = new TF1("fierz_fit", fit_str, min_E, max_E);
    fierz_fit->SetParameter(0,-0.12);
    fierz_fit->SetParameter(1,0);
	asymmetry_histogram->Fit(fierz_fit, "Sr");

	// A fit histogram for output to gnuplot
    TH1F *fierz_fit_histogram = new TH1F(*asymmetry_histogram);
	for (int i = 0; i < fierz_fit_histogram->GetNbinsX(); i++)
		fierz_fit_histogram->SetBinContent(i, fierz_fit->Eval(fierz_fit_histogram->GetBinCenter(i)));

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
	asymmetry_histogram->SetStats(0);
    asymmetry_histogram->Draw();

	// draw a legend on the plot
    TLegend* ratio_legend = new TLegend(0.3,0.85,0.6,0.65);
    ratio_legend->AddEntry(asymmetry_histogram, "Asymmetry data", "l");
    ratio_legend->AddEntry(fierz_fit, "Fierz term fit", "l");
    //ratio_legend->AddEntry((TObject*)0, "1+b(#frac{m_{e}}{E} - #LT #frac{m_{e}}{E} #GT)", "");
    ratio_legend->AddEntry((TObject*)0, A_str, "");
    ratio_legend->AddEntry((TObject*)0, b_str, "");
    ratio_legend->AddEntry((TObject*)0, chisq_str, "");
    ratio_legend->SetTextSize(0.03);
    ratio_legend->SetBorderSize(0);
    ratio_legend->Draw();

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
