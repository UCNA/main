// c includes 
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
//#include <gsl/gsl_poly.h>

// ROOT includes
#include <TF1.h>
#include <TChain.h>
#include <TString.h>
#include <TCut.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH2F.h>
#include <TStyle.h>
#include <TApplication.h>

#include <iostream>
#include <fstream>
#include <TPaveStats.h>
#include <TPave.h>
#include <TAttText.h>
#include <TList.h>
using namespace std;

// Fitter includes
/*#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Fit/Chi2FCN.h"
#include "TList.h"
#include "Math/WrappedMultiTF1.h"
#include "HFitInterface.h"
*/

#include <TMinuit.h>

/**
 *
 * Authors: Kevin Peter Hickerson
 *          Michael Mendenhall
 *          Simon Slutsky
 * Date: Aug 2010
 * Modified: July 16, 2011
 * Modified: Feb 25, 2013
 * Modified: Sep 30, 2013 (SS)
 * Branched from pd_led_pmt.cc Apr 14, 2015 (SS)
 *    Will implement a combined fit to 405 and 465 graphs
 *
 * Build instructions:
USE THIS:
 *  g++ `root-config --cflags` pd_led_pmt_combinedfit.cc `root-config --libs` -lMinuit -o pd_led_pmt_combinedfit_analysis
FOLLOWING DOESN'T WORK:
 *  g++ `root-config --cflags` `root-config --libs` pd_led_pmt.cc -o pd_led_pmt_analysis
 *
 * If GSL is needed use:
 *  g++ `root-config --cflags` `root-config --libs` -lgsl -lgslcblas -lm pd_led_pmt.cc -o pd_led_pmt_analysis
 *
 */

#define NUM_CHANNELS 8
#define LED_TYPE DOWN
#define USE_ROOT_APPLICATION false
#define OUTPUT_IMAGE true
#define OUTPUT_IMAGE_DIR "/data1/saslutsky/LEDPulser/images_05_04_2015_16way_linearPMT_fit_longer_range_raw_output_18745_18768/"  // DON'T OMIT THE TRAILING SLASH
#define VERBOSE true
#define LINEARIZE false
#define ORDER 2 // Power law fit
//#define ORDER 3 // Kevin's default
#define FIT2D false
#define DO_LED_FIT 1
#define POLYFIT 0
#define OFFSETPOLYFIT 1
//#define CONSTRAINFIT 1
#define BETAADCCOUNTS 782
#define PLOTBOTHLEDS 1
#define RANGE_MIN 5.0
#define MOVEDUPGRAPH false
#define KEVSCALED false
#define RANGE_MAX_OVERRIDE false
#define RANGE_MAX_VALUE 100.0 // only if RANGE_MAX_OVERRIDE = true
#define FIXBETAENDPOINT true
#define RELATIVEBETAPLOTS false  // supersedes FIXBETAENDPOINT (and everything else)
#define COMBINEDLED 1 // replace later with a loop
#define USEBETAOFFSETS false
#define RANGE_EXTENSION 3.0

const int pulser_steps = 64;

//Declare data-arrays globally for fits
vector<float> gPMT[2][NUM_CHANNELS];
vector<float> gPD[2][NUM_CHANNELS];
vector<float> gPMTerr[2][NUM_CHANNELS];
vector<float> gPDerr[2][NUM_CHANNELS];
 // need a fixed parameter for origin of PD or fits won't make sense; value is arbitrary "average by eye" of PMT beta endpoints in PD value
//float gPDoff[2] = {30.0, 100.0};
float gPDoff[2] = {0.0, 100.0};
//float gPDoff = 0.0;

// need global array of beta-endpoints so we can include these in fits
vector<float> gPDBetaEP[2]; // one vector for each wavelength 
vector<float> gPMTBetaEP;  // tube property, same for both wavelengths
// for some reason need a global to indicate which LED we're using in the fit
//float gLED = 0;

// temp arrays to hold data until cuts can be made
vector<float> _gPMT[2][NUM_CHANNELS];
vector<float> _gPD[2][NUM_CHANNELS];
vector<float> _gPMTerr[2][NUM_CHANNELS];
vector<float> _gPDerr[2][NUM_CHANNELS];

Double_t func(float gPDval, Double_t *par, Int_t i, Int_t led);
//Double_t PDfunc(float gPDval, Double_t *par, Int_t led, Int_t i);
Double_t subfcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t led, Int_t iflag);
Double_t combiErr(Double_t * par, Int_t i, Int_t led, Int_t k);
Double_t PDInterperr(Double_t * par, Int_t i, Int_t led, Int_t k);

// calculate chi^2 
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{ 
  Double_t value = 0;
  for (int ded = 0; ded < 2; ded++){   // sum chi^2 for each led 
    value += subfcn(npar, gin, f, par, ded, iflag);
  }

  f = value;
}  

Double_t subfcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t led, Int_t iflag)
{
  Double_t PMT_term = 0;
  Double_t PD_term = 0;
  Double_t delta = 0;
  Double_t chisq = 0;
  
  //  int led = COMBINEDLED; // replace later with loop over LEDs
  //int led = gLED; // no longer needed now that we're looping
  for (int i = 0; i < NUM_CHANNELS; i++){
    // for (int i = 0; i < 1; i++){
    Double_t _chisq_temp = 0;
    for (int k=0; k<gPD[led][i].size(); k++){
      // only uses PMT errors - should redo to include PD errors (see TGraph::Fit() Doc)
      if (gPMTerr[led][i][k] < 0.000000001) continue;
      if (gPDerr[led][i][k] < 0.000000001) continue;
      //      delta = (gPMT[led][i][k] - func(gPD[led][i][k], par, i, led))/gPMTerr[led][i][k];
      //_chisq_temp = delta*delta;
      // include PD errors:
      delta = (gPMT[led][i][k] - func(gPD[led][i][k], par, i, led));
      _chisq_temp += delta*delta/combiErr(par, i, led, k); // err already squared;
      //      if (combiErr(par, led, i, k) < 0.0000000001){
      //	cout << "Error <= 0:  " << combiErr(par, led, i, k) << endl;
      //	continue;}
    }
    chisq += _chisq_temp;
  }
  //  cout << "Chi^2 = " << chisq << endl;
  return chisq;
}

Double_t combiErr(Double_t * par, Int_t i, Int_t led, Int_t k){
  Double_t combinederr, PDerr, PMTerr;
  PMTerr = gPMTerr[led][i][k];
  PDerr = PDInterperr(par, i, led, k);
  combinederr = PMTerr*PMTerr + PDerr*PDerr;
  return combinederr;
}

// from TGraph doc page: total error term when x and y both have errors is = 
// #frac{(y-f1(x))^{2}}{ey^{2}+(#frac{1}{2}(exl+exh)f1'(x))^{2}}

// works for a quadratic PD with linear PMT
Double_t PDInterperr(Double_t * par, Int_t i, Int_t led, Int_t k){ 
  Double_t pde = gPDerr[led][i][k];
  //  cout << "led: " << led << endl;
  if (pde == 0) cout << "0000000 " << led << " " << i << " " << k << endl;
  Double_t scale = 0;
  if (led == 0) scale = par[18];
  if (led == 1) scale = 1.0;
  Double_t deriv = (1./par[1+2*i]) * scale*(par[16] + 2*par[17]*gPD[led][i][k]);
  return pde*deriv;
}  

// par[0] - par[15] --> PMTs offset, lin;  par[16]-par[18] --> PD lin, quad, scale
Double_t func(float gPDval, Double_t *par, Int_t i, Int_t led)
{
  Double_t scale = 0;
  if (led == 0) scale = par[18];
  if (led == 1) scale = 1.0;
  Double_t value = (1./par[1+2*i]) * ( -par[0+2*i] + scale*(par[16]*gPDval + par[17]*gPDval*gPDval) );
  return value;
}

/* // par[0] - par[23] --> PMTs, par[24]-par[26] --> PD, gPDBetaEP[led][i] = beta endpoint in PD units for tube 
Double_t func(float gPDval, Double_t *par, Int_t i, Int_t led)
{
  Double_t PDval = PDfunc(gPDval, par, led, i);
  Double_t value = par[0+3*i] + par[1+3*i]*(PDval - gPMTBetaEP[i]) + par[2+3*i]*(PDval - gPMTBetaEP[i])*(PDval - gPMTBetaEP[i]); // Use PMTEP since PDval is approx linear with PMT response (??)
  return value;
  }*/

  //  Double_t value = par[0+3*i] + par[1+3*i]*(par[24] + par[25]*(gPD-BetaEP[led][i]) + par[26]*(gPD-BetaEP[led][i])*(gPD-BetaEP[led][i])) + par[2+3*i]*(par[24] + par[25]*(gPD-BetaEP[led][i]) + par[26]*(gPD-BetaEP[led][i])*(gPD-BetaEP[led][i]))*(par[24] + par[25]*(gPD-BetaEP[led][i]) + par[26]*(gPD-BetaEP[led][i])*(gPD-BetaEP[led][i]));
  //  Double_t value = par[0+3*i] + par[1+3*i]*gPD + par[2+3*i]*gPD*gPD;

 /*Double_t PDfunc(float gPDval, Double_t *par, Int_t led, Int_t i)
{
  // a quadratic function of PD
  //  Double_t PDvalue = gPD;
  //  Double_t PDvalue = par[24] + par[25]*(gPDval - gPDoff[led]) + par[26]*(gPDval - gPDoff[led])*(gPDval - gPDoff[led]);
  Double_t PDvalue = par[24] + par[25]*(gPDval - gPDBetaEP[led][i]) + par[26]*(gPDval - gPDBetaEP[led][i])*(gPDval - gPDBetaEP[led][i]);
  if (led == 0) PDvalue *= par[27]; // scaling factor between two LED wavelengths
  //  if (led == 0) PDvalue *= 3; // scaling factor between two LED wavelengtsh
  return PDvalue;
  } */


TF1* FitGaussian(const char *name, TTree *tree, TCut* cut)
{
  char gaussian_name[1024];
  sprintf(gaussian_name, "pedestal_histogram_%s", name);
  char gaussian_draw[1024];
  sprintf(gaussian_draw, "%s >> %s", name, gaussian_name);
  TH1F* gaussian_histogram = new TH1F(gaussian_name, "Gaussian Events", 2000, 0, 2000);
  tree->Draw(gaussian_draw, *cut, "Ngoff"); 
  int max_bin = gaussian_histogram->GetMaximumBin();
  float max_bin_x = gaussian_histogram->GetBinCenter(max_bin);
  TF1 *fit = new TF1("gauss_fit", "gaus", max_bin_x-12, max_bin_x+12);
  if (not gaussian_histogram->Fit(fit, "RN"))
    {
      printf("Gaussian fit success: mu = %g, sigma = %g\n", fit->GetParameter(1), fit->GetParameter(2));
      return fit;
    } 
  else 
    {
      printf("Couldn't fit Gaussian to %s\n", name);
      return 0;
    }
}


/**
 * main
 */
int main (int argc, char **argv)
{
  // first run number
  if (argc < 2)
    {
      printf("Usage: %s <daq run number> [<second run number>]\n", argv[0]);
      exit(1);
    }
  int run = atoi(argv[1]);
  int run_end = run;
  if (not run)
    {
      printf("Need a valid run number for first argument.\n");
      printf("Usage: %s <daq run number> [<second run number>]\n", argv[0]);
      exit(1);
    }

#if VERBOSE
  cout << "Starting run number " << run << endl;
#endif
  
  // next run number for a range
  if (argc > 2)
    {
      run_end = atoi(argv[2]);
      if (not run_end)
        {
	  printf("Need a valid run number for second argument.\n");
	  printf("Usage: %s <daq run number> [<second run number>]\n", argv[0]);
	  exit(1);
        }
    }
  
  
  // constants (maybe should be set by options too?)
  //enum {DOWN, UP, GAIN, TIME, TIMEUP, TIMEDOWN, TIMEGAIN};
  enum {DOWN, UP, GAIN, TIME};
  const int max_cycles = 100; // more than an hour
  int n = max_cycles * pulser_steps * 100;  // this can be updated later to be more accurate
  //	const float max_pdc_channel = 400;
  const float max_pdc_channel = 800;
  const float max_adc_channel = 3400;
  const float mid_adc_channel = 3000;
  const float min_adc_channel = 100;
  const float max_npe = 400;
  const float epsilon = 8e4; // us
  const float window = 3e4; // us
  const float avg_period = 60.03e6; // us
  const float max_period = avg_period + window; // us
  const float min_period = avg_period - window; // us
  
  
  // some strings we will use
  TString wavelength[2] = { "405", "465" };
  TString detector[8] = { "E1", "E2", "E3", "E4", "W1", "W2", "W3", "W4" };
  TString Qadc[8] = { "Qadc0", "Qadc1", "Qadc2", "Qadc3", "Qadc4", "Qadc5", "Qadc6", "Qadc7"};


  vector <float> best_beta_endpt;
  vector <float> best_beta_endpt_PD[2]; // one for each led
  //  float extended_range_max[2][NUM_CHANNELS]; // deprecated
  float range_max[2][NUM_CHANNELS], range_min[2][NUM_CHANNELS];

  // run this as a ROOT application
#if USE_ROOT_APPLICATION
  TApplication app("LED Scans", &argc, argv);
#endif
  

  // Plot options
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1);
  

  // Open files
  // example paths for various computers:
  // ferrari.caltech.edu: /home/kevinh/Data/UCN/UCNA2011/rootfiles
  // pcucn19.caltech.edu: /home/data_analyzed/2011/rootfiles
  // ucn.caltech.edu: 	/data/ucnadata/2011/rootfiles
  TChain h1("h1");
  for (int n = 1; n < argc; n++) {
    char * rootfiles;
    char filename[1024];
    rootfiles = getenv("UCNADATADIR");
    if (not rootfiles)
      {
	printf("Environment variable UCNADATADIR not set. Can't find root files.\n");
	printf("try setting it with export. Here are some examples for different machines\n");
	printf("On the machine pcucn18.caltech.edu, try\n");
	printf(" export UCNADATADIR=/home/data_analyzed/2011/rootfiles\n");
	printf("On the machine ucn.caltech.edu, tryi\n");
	printf(" export UCNADATADIR=/data/ucnadata/2011/rootfiles\n");
	exit(1);
      }
    sprintf(filename, "%s/full%s.root", rootfiles, argv[1]);
    if (not h1.Add(filename))
      {
	printf("Could not find root file full%s.root\n", argv[1]);
	printf("The path set by UCNADATADIR could be set to the wrong dir\n");
	printf("or the run number does not exist in that dir.");
	exit(1);
      }
  }
  
  
  // Initialize the cycles peaks
  float max_val[max_cycles];
  float max_time[max_cycles];
  for (int cycle = 0; cycle < max_cycles; cycle++)
    {
      max_val[cycle] = 0;
      max_time[cycle] = 0;
    }
  
  
  // Sync peaks and find periods
  {
    // set the branches we need
    h1.SetBranchStatus("*",0);
    h1.SetBranchStatus("Sis00", 1);
    h1.SetBranchStatus("S83028", 1);
    h1.SetBranchStatus("Pdc36", 1);
    
    n = h1.GetEntries();
    h1.SetEstimate(n);
    int vn = h1.Draw("Pdc36:S83028",
		     "(int(Sis00) & 128) > 0", "goff");
    
    double *pdc36v = h1.GetV1();
    double *s83028v = h1.GetV2();
    std::cout << "Number of LED events: " << vn << std::endl;
    
    double start_time = 0;//s83028v[0];
    int start_index = 0;
    const int repeat = 8; 
    for (int i = vn - 1; i >= start_index + repeat; i--) 
      {
	double time = s83028v[i] - start_time;
	int cycle = time/avg_period;
	bool max_found = true;
	
	// hack to fix first cycle syncing problem
	if (cycle == 0) 
	  break;
	
	// need a few large values in a row
	for (int j = 0; j < repeat; j++)
	  {
	    if (pdc36v[i-j] <= max_val[cycle])
	      max_found = false;
	  }
	
	// then set the max found for this cycle
	if (max_found)
	  {
	    max_val[cycle] = pdc36v[i]; 
	    max_time[cycle] = s83028v[i]; 
	  }
      }
    
#if VERBOSE
    std::cout << "Start time: " << start_time << "\n";
    std::cout << "1st max PD (Pdc36): " << max_val[0] << "\n";
    std::cout << "1st start time (S83028): " << max_time[0] / 1e6 << "\n";
    for (int cycle = 1; cycle < max_cycles; cycle++)
      {
	if (max_val[cycle] == 0)
	  break;
	
	std::cout << cycle + 1 << "th max PD (Pdc36): " << max_val[cycle] << "\n";
	std::cout << cycle + 1 << "th start time (S83028): " << max_time[cycle]/ 1e6  << "\n";
	std::cout << cycle << "th Period: " << (max_time[cycle] - max_time[cycle-1]) / 1e6 << "\n";
      }
#endif
  }
  
  // Define a histo for the true PD signal history (before cycle parsing)
  TString true_time_title = "True PD time sequence for run ";
  true_time_title += run;
  TH2F* true_time_his2D = new TH2F("true_time_his2D", true_time_title,
				   360000., 0, 3600.,
				   1<<8, -100, 1000);
  h1.Draw("Pdc36:S83028/1000000.>>true_time_his2D", "(int(Sis00) & 128) > 0", "goff");
  
  // Define cuts
  TCut *led_cut = new TCut("(int(Sis00) & 128) > 0");  // 129 if east-PMTs on, 161 if GMS-ref also on
  TCut *pedestal_cut = new TCut("!(int(Sis00) & 1)");  
  
  
  // The histograms and canvases we will use 
  TH2F* time_his2D[4];
  TH1F* pmt_gain_his1D[8];
  TH2F* pmt_gain_his2D[8];
  TH1F * time_his1D;
  TH2F* pd_pmt_his2D[2][8];
  TGraphErrors* graph[2][8];
  TGraphErrors* graph_465_moved_up[8];
  TGraphErrors* PMT_keV_graph[2][8];
  TGraphErrors* constrained_graph[2][8];
  TCanvas* canvas[8];
  TGraph* g[8];
  TGraphErrors* g_PE_PMT[8];
  TGraphErrors* resg[8];
  TCanvas * PE_PMT_canvas = new TCanvas();
  PE_PMT_canvas->Divide(2,4);

  // global Style options
  gStyle->SetPalette(1);
  gStyle->SetOptStat("");
  
  
  // find PD pedestal
  h1.SetBranchStatus("Pdc36", 1);
  TF1 *pd_pedestal_fit = FitGaussian("Pdc36", &h1, pedestal_cut);
  float pd_pedestal = 0;
  if (pd_pedestal_fit)
    pd_pedestal = pd_pedestal_fit->GetParameter(1);
  else
    printf("couldn't fit pedestal to Pdc36\n");
  printf("PD pedestal = %f\n", pd_pedestal);
  
  
  // Making our histograms for the time sequence
  TString time_title = "PD time sequence for run ";
  time_title += run;
  time_his2D[TIME] = new TH2F("H2F_time_name", time_title, 
			      3600*10, 0, 3600,
			      1<<8, -pd_pedestal, 500-pd_pedestal);
  
  time_his2D[UP] = new TH2F("H2F_time_up_name", time_title, 
			    3600*10, 0, 3600,
			    1<<8, -pd_pedestal, 500-pd_pedestal);
  
  time_his2D[DOWN] = new TH2F("H2F_time_down_name", time_title, 
			      3600*10, 0, 3600,
			      1<<8, -pd_pedestal, 500-pd_pedestal);
  
  time_his2D[GAIN] = new TH2F("H2F_time_gain_name", time_title, 
			      3600*10, 0, 3600,
			      1<<8, -pd_pedestal, 500-pd_pedestal);
  time_his1D = new TH1F("H1F_time_gain_name_1D", time_title,
			100, 200, 300);
  
  // prepare a file to store linearity fit parameters
  //  double * fitpars = 0; 
  //  double * fitparerrors = 0;
  double ** fitpars = new double*[2];
  double ** fitparerrors = new double*[2];
  double fitchisq = 0;
  int Npar;

  double * pefitpars = 0; 
  double * pefitparerrors = 0;
  double pefitchisq = 0;
  int peNpar;

  TString fitfilename = "FitResults.txt";
  //	#if CONSTRAINFIT
  TString constrainfitfilename = "FitResults_constrained.txt";
  //	#endif
  fitfilename = OUTPUT_IMAGE_DIR + fitfilename;
  constrainfitfilename = OUTPUT_IMAGE_DIR + constrainfitfilename;

  ofstream fitfile;
  fitfile.open(fitfilename, std::ofstream::out | std::ofstream::app);
  //  ofstream constrainfitfile;
  //constrainfitfile.open(constrainfitfilename, std::ofstream::out | std::ofstream::app);
  // Can't get the string to write for only first file
  /*TString fit_header_string = "Run\t";
    fit_header_string += "Channel\t";
    fit_header_string += "Wavelength\t";
    fit_header_string += "p0\tp0Err\t";
    fit_header_string += "p1\tp1Err\t";
    fit_header_string += "p2\tp2Err\t";
    fit_header_string += "Chi2\t";
    fit_header_string += "\n";
    ftfile << fit_header_string;
  */
  
  TString pepmtfitfilename = "FitResults_PE_PMT.txt";
  pepmtfitfilename = OUTPUT_IMAGE_DIR + pepmtfitfilename;
  ofstream pepmtfitfile; 
  pepmtfitfile.open(pepmtfitfilename, std::ofstream::out | std::ofstream::app);

  TString combfitfilename = "FitResults_Combined.txt";
  combfitfilename = OUTPUT_IMAGE_DIR + combfitfilename;
  ofstream combfitfile; 
  combfitfile.open(combfitfilename, std::ofstream::out | std::ofstream::app);
  
  TString convfactorfilename = "PD_keV_conversion_factors.txt";
  convfactorfilename = OUTPUT_IMAGE_DIR + convfactorfilename;
  ofstream convfactorfile;
  convfactorfile.open(convfactorfilename, std::ofstream::out | std::ofstream::app);

  TString rawfilename = "RawGraphData.txt";
  rawfilename = OUTPUT_IMAGE_DIR + rawfilename;
  ofstream rawfile;
  rawfile.open(rawfilename, std::ofstream::out | std::ofstream::app);
      
  float gain_sum = 0;
  float gain2_sum = 0;
  int gain_cnt = 0;
  
  // Max Beta ADC Channels for given run segments 
  // will use later for determing run ranges
  //  float beta_Cd_ratio = 63./782.;
  //  float beta_Bi_ratio = 1047./782.;
  float beta_Cd_ratio = 63./BETAADCCOUNTS;
  float beta_Bi_ratio = 1047./BETAADCCOUNTS;

  //< 20500
  float BetaADC_below_20500[NUM_CHANNELS] = {500, 500, 850, 500, 700, 550, 700, 1050};
  // 20500-21250
  float BetaADC_20500_21250[NUM_CHANNELS] = {550, 400, 850, 500, 450, 600, 500, 1050};
  // >21250
  float BetaADC_21250_above[NUM_CHANNELS] = {650, 650, 1200, 550, 750, 700, 900, 1000};
  
  // go through each PMT channel
  for (unsigned i = 0; i < NUM_CHANNELS; i++) 
    {
      // set up only the branches we use
      h1.SetBranchStatus("*", 0);
      h1.SetBranchStatus(Qadc[i], 1);
      h1.SetBranchStatus("Sis00", 1);
      h1.SetBranchStatus("S83028", 1);
      h1.SetBranchStatus("Pdc36", 1);
      
      
      // Find PMT Qadc Pedestal
      TF1 *pmt_pedestal_fit = FitGaussian(Qadc[i], &h1, pedestal_cut);
      float pmt_pedestal = 0;
      if (pmt_pedestal_fit)
	pmt_pedestal = pmt_pedestal_fit->GetParameter(1);
      
      
      // Define PMT:PD scan histograms
      TString H2F_down_name = "H2F_down_name_";
      H2F_down_name += detector[i];
      pd_pmt_his2D[DOWN][i] = new TH2F(H2F_down_name, "", 
				       200, 0, 200,
				       1<<8, -pmt_pedestal, 4096-pmt_pedestal);

      TString H2F_up_name = "H2F_up_name_";
      H2F_up_name += detector[i];
      pd_pmt_his2D[UP][i] = new TH2F(H2F_up_name, "", 
				     //				400, 0, 400,	
				     800, 0, 800,
				     1<<8, -pmt_pedestal, 4096-pmt_pedestal);
      
      TString time_gain_his1D_name = "time_gain_1D_";
      TString time_gain_his2D_name = "time_gain_2D_";
      time_gain_his1D_name += detector[i];
      time_gain_his2D_name += detector[i];
      
      pmt_gain_his1D[i] = new TH1F (time_gain_his1D_name, "", 
				    300, 0, 3000);
      pmt_gain_his2D[i] = new TH2F (time_gain_his2D_name, "", 
				    3600*10, 0, 3600,
				    300, 0, 3000);
      
      pd_pmt_his2D[DOWN][i]->GetXaxis()->SetTitle("Photodiode");
      pd_pmt_his2D[DOWN][i]->GetYaxis()->SetTitle("ADC");
      pd_pmt_his2D[UP][i]->GetXaxis()->SetTitle("Photodiode");
      pd_pmt_his2D[UP][i]->GetYaxis()->SetTitle("ADC");
      pmt_gain_his1D[i]->GetXaxis()->SetTitle("PMT ADC counts");
      
      TString draw_time_cmd = "S83028:Pdc36:";
      draw_time_cmd += Qadc[i];
      
      h1.SetEstimate(n);
      int vn = h1.Draw(draw_time_cmd, "(int(Sis00) & 128) > 0", "goff");
      double *s83028v = h1.GetV1();
      double *pdc36v = h1.GetV2();
      double *qadcv = h1.GetV3();
      std::cout << "Number of LED events: " << vn << std::endl;
      
      //int start_time = max_time[i];
      int start_time = 0;
      int last_time = start_time;
      int cycle = 0;
      int subcycle = 0;
      
      float x_sum[2][pulser_steps];
      float y_sum[2][pulser_steps];
      float x2_sum[2][pulser_steps];
      float y2_sum[2][pulser_steps];
      int num_pulses[2][pulser_steps];
      
      for (int led = 0; led < 2; led++)
	{
	  for (int pulse = 0; pulse < pulser_steps; pulse++)
	    {
	      x_sum[led][pulse] = 0;
	      y_sum[led][pulse] = 0;
	      x2_sum[led][pulse] = 0;
	      y2_sum[led][pulse] = 0;
	      num_pulses[led][pulse] = 0;
	    }
	}
      
      for (int j = 0; j < vn; j++) 
	{
	  if (cycle+1 == max_cycles or max_val[cycle+1] == 0)
	    break;
	  
	  float period = max_time[cycle+1] - max_time[cycle];
	  float time = s83028v[j] - start_time;
	  float cycle_time = s83028v[j] - last_time;
	  if (cycle_time > period or period > max_period or period < min_period)
	    {
	      cycle++;
	      subcycle = 0;
	      cycle_time -= last_time;
	      last_time = max_time[cycle];
	    }
	  
	  if (min_period < period and period < max_period)
	    {
	      float subperiod = period / pulser_steps;
	      float subcycle_time = cycle_time - subperiod * subcycle;
	      if (subcycle_time > subperiod)
		{
		  subcycle++;
		  subcycle_time -= subperiod;
		}
	      
	      float x = pdc36v[j] - pd_pedestal;
	      float y = qadcv[j] - pmt_pedestal;
	      
	      /*#if LINEARIZE
		#if ORDER
		// back out linearity from pd
		//double a[3] = {0, 7.601/7.601, -0.005296/7.601};
		double a[4] = {0, 8.19/8.19, -0.008783/8.19, 4.14e-6/8.19};
		
		float _x = a[ORDER];
		for (int k = ORDER; k > 0; k--)
		_x = a[k-1] + _x * x; 
		x = _x;
		
		#endif
		#endif
	      */
	      if (epsilon < subcycle_time and subcycle_time < subperiod/2 - epsilon) // ramps
		{
		  if (not i){               // only fill once
		    time_his2D[GAIN]->Fill(time/1e6, x);
		    time_his1D->Fill(x);
		    gain_sum += x;
		    gain2_sum += x*x;
		    gain_cnt++;
		  }
		  // PMT response during gain period
		  pmt_gain_his1D[i]->Fill(y);
		  pmt_gain_his2D[i]->Fill(time/1e6, y);
		}
	      if (subperiod / 2 + epsilon < subcycle_time and subcycle_time < subperiod - epsilon) // gain
		{
		  int led = (cycle_time > period/2);
		  if (cycle_time > period / 2) // ramp up (465nm?)
		    {
		      if (not i)
			time_his2D[UP]->Fill(time/1e6,x);
		      pd_pmt_his2D[UP][i]->Fill(x, y);
		      x_sum[UP][subcycle] += x;
		      y_sum[UP][subcycle] += y;
		      x2_sum[UP][subcycle] += x*x;
		      y2_sum[UP][subcycle] += y*y;
		      num_pulses[UP][subcycle] ++;
		    }
		  if (cycle_time < period / 2) // ramp down (405nm?)
		    {
		      if (not i)        // Only Fill once
			time_his2D[DOWN]->Fill(time/1e6,x);
		      pd_pmt_his2D[DOWN][i]->Fill(x, y);
		      x_sum[DOWN][subcycle] += x;
		      y_sum[DOWN][subcycle] += y;
		      x2_sum[DOWN][subcycle] += x*x;
		      y2_sum[DOWN][subcycle] += y*y;
		      num_pulses[DOWN][subcycle] ++;
		    }
		}
	      time_his2D[TIME]->Fill( time/1e6, x);
	    }
	}

      // set the graph points to the averages of the led cloud
      // TODO might be better to actually fit the peak better
#if FIT2D
      printf("Coming Soon");
      return 0;
#else
      
      for (int led = 0; led < 2; led++)
	{
	  graph[led][i] = new TGraphErrors(pulser_steps);
	  for (int pulse = 0; pulse < pulser_steps; pulse++)
	    {
	      float d = num_pulses[led][pulse];
	      if (d > 0)
		{
		  float x_avg = x_sum[led][pulse]/d;
		  float y_avg = y_sum[led][pulse]/d;
		  float x2_avg = x2_sum[led][pulse]/d; 
		  float y2_avg = y2_sum[led][pulse]/d;
		  float sx = sqrt((x2_avg - x_avg*x_avg)/(d-1));
		  float sy = sqrt((y2_avg - y_avg*y_avg)/(d-1));
		  //		  cout << "sx " << x2_avg << ", " << x_avg*x_avg << ", " << sx << endl;
		  //cout << sy << endl;
		  graph[led][i]->SetPoint(pulse, x_avg, y_avg);
		  graph[led][i]->SetPointError(pulse, sx, sy);
		  // pipe values out to global so fit function can access
		  //if (!led){ 		  // just do 405 nm for now, fix later
		  _gPD[led][i].push_back(x_avg);
		  _gPMT[led][i].push_back(y_avg);
		  _gPDerr[led][i].push_back(sx);
		  _gPMTerr[led][i].push_back(sy);
		  //}
		}
	    }
	}
#endif
                  
      //Use approx Maximum beta-endpoints over each run interval to 
      // determine the best run range. Beta endpoint = 782 keV, Cd-109 = 63 keV, Bi-207 = 1047 keV	  
      //      float best_beta_endpt; // make a vector to store for later loops
      float ADC_max, ADC_min;
      float nPE_max, nPE_min;
      double slope_DOWN; double slope_UP;
      float pd_to_keV_factor[2];

      /*      if (run < 20500) best_beta_endpt = BetaADC_below_20500[i];
      else if (run > 20500 && run < 21250) best_beta_endpt = BetaADC_20500_21250[i];
      else if (run > 21250) best_beta_endpt = BetaADC_21250_above[i];*/
      if (run < 20500) best_beta_endpt.push_back(BetaADC_below_20500[i]);
      else if (run > 20500 && run < 21250) best_beta_endpt.push_back(BetaADC_20500_21250[i]);
      else if (run > 21250) best_beta_endpt.push_back(BetaADC_21250_above[i]);
      else cout << "RUN NOT IN RANGE" << endl;
      ADC_min = best_beta_endpt[i]*beta_Cd_ratio;
      ADC_max = best_beta_endpt[i]*beta_Bi_ratio;
      cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
      cout << ADC_min << "=ADC_min " << ADC_max << "=ADC_max" << endl;
      cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

#if USEBETAOFFSETS
      gPMTBetaEP.push_back(best_beta_endpt[i]);
#else
      gPMTBetaEP.push_back(0);
#endif

      // find the best range for the fit
      //      float range_max[2], range_min[2];// extended_range_max[2]; needs wider scope
      
      /*		float prefit_range_max[2], prefit_range_min[2];
			prefit_range_min[0] = 25; prefit_range_min[1] = 50;
			prefit_range_max[0] = 150; prefit_range_max[1] = 300;
      */
      
      
      //	printf("found range (%f,%f)\n", range_min, range_max); 
      
      
      // Find # photoelectrons (nPE) corresponding to range min and max
      for (int led = 0; led < 2; led++) 
	{
	  // Find range in PD counts corresponding to ADC range
	  printf("Finding range... ");
	  float _max = 0, _min = 0;
	  range_max[led][i] = 0; range_min[led][i] = 0;
	  for (int k = 0; k < pulser_steps; k++)
	    {
	      double x,y;
	      if (led){
		graph[led][i]->GetPoint(k, x, y);
	      }
	      if (!led){ // 405 pulser counts backwards from large signal to small
		graph[led][i]->GetPoint(pulser_steps - k - 1, x, y);
	      }
	      //	      cout << x << "  " << y << endl;
	      if (y < ADC_min){_min = x;}
	      if (y < ADC_max){_max = x;}
	      if (y > ADC_max){break;}

	      range_max[led][i] = _max;
	      range_min[led][i] = _min;
	    }
	  if (range_min[led][i] < RANGE_MIN ) range_min[led][i] =RANGE_MIN;  // avoid bad clump of below-threshold LED shots
	  // Fudge factor to slightly extend range for better fits
	  //	  range_min[led] = range_min[led]*0.9; 
	  range_max[led][i] = range_max[led][i]*1.05;
	  //	  if (led == 0) extended_range_max[led][i] = range_max[led]*2;
	  //	  else extended_range_max[led][i] = extended_range_max[0][i]; //see if lowering 465 helps
	  //	  extended_range_max[led][i] = RANGE_EXTENSION*range_max[led];
	  range_max[led][i] = RANGE_EXTENSION*range_max[led][i];
	  
#if RANGE_MAX_OVERRIDE
	  range_max[led][i] = RANGE_MAX_VALUE;
#endif
#if RELATIVEBETAPLOTS
	  range_min[1][i] = range_min[0][i]; // fix range_min to be same for both LED when comparing
#endif
	  printf("found range (%f,%f)\n", range_min[led][i], range_max[led][i]); 
	      
	  // Use best value from ADC max to create a scaling factor 
	  // and make a plot of PMT counts vs "keV", where "keV" is 
	  // estimated based on Beta endpoint in PMT counts
	  // Write this factor to file for later use

	  pd_to_keV_factor[led] = beta_Bi_ratio*BETAADCCOUNTS/range_max[led][i]; 
	  
	  TString out_factor_string = "";
	  out_factor_string += run;                out_factor_string += "\t"; 
	  out_factor_string += i;                  out_factor_string += "\t"; 
	  out_factor_string += wavelength[led];    out_factor_string += "\t"; 
	  out_factor_string += pd_to_keV_factor[led];   out_factor_string += "\t"; 
	  out_factor_string += "\n";
	  
	  convfactorfile << out_factor_string;

	  // = (1047./782.)*782 = 1047 = upper Bi peak in keV; applying factor to range_max gives 1047
	  PMT_keV_graph[led][i] = new TGraphErrors(pulser_steps);
	  double _pdc_in_keV = 0;
	  double _pdc, _pmt, _pdc_e, _pmt_e; 
	  for (int pulse = 0; pulse < pulser_steps; pulse++)
	    {
	      graph[led][i]->GetPoint(pulse, _pdc, _pmt);
	      _pdc_e = graph[led][i]->GetErrorX(pulse);
	      _pmt_e = graph[led][i]->GetErrorY(pulse);
	      
	      _pdc_in_keV = _pdc*pd_to_keV_factor[led];
	      PMT_keV_graph[led][i]->SetPoint(pulse, _pdc_in_keV, _pmt);
	      PMT_keV_graph[led][i]->SetPointError(pulse, _pdc_e*pd_to_keV_factor[led], _pmt_e);
	    }

	  // convert range from PD to keV // Replaced with Energy Scaled graph 1/13/2015 
	  cout << range_min[led][i] << " HIH " << range_max[led][i] << endl;
#if KEVSCALED
	  range_min[led][i] = range_min[led][i]*pd_to_keV_factor[led][i];
	  range_max[led][i] = range_max[led][i]*pd_to_keV_factor[led][i];
#endif
	  cout << range_min[led][i] << " HIH" << range_max[led][i] << endl;

	  // cull global vectors to match ranges
	  for (int s=0; s < _gPD[led][i].size(); s++){
	    // write all the data to file
	    TString out_raw_string = "";
	    out_raw_string += run;                   out_raw_string += "\t"; 
	    out_raw_string += i;                     out_raw_string += "\t"; 
	    out_raw_string += led;                   out_raw_string += "\t"; 
	    out_raw_string += _gPD[led][i][s];       out_raw_string += "\t"; 
	    out_raw_string += _gPDerr[led][i][s];    out_raw_string += "\t"; 
	    out_raw_string += _gPMT[led][i][s];      out_raw_string += "\t"; 
	    out_raw_string += _gPMTerr[led][i][s];   out_raw_string += "\t"; 
	    out_raw_string += "\n";
	    rawfile << out_raw_string;
	    if (_gPD[led][i][s] > range_min[led][i] && _gPD[led][i][s] < range_max[led][i]){
	    //	    if (_gPD[led][i][s] > range_min[led] && _gPD[led][i][s] < extended_range_max[led][i]){
	      gPD[led][i].push_back(_gPD[led][i][s]);
	      gPMT[led][i].push_back(_gPMT[led][i][s]);
	      gPDerr[led][i].push_back(_gPDerr[led][i][s]);
	      gPMTerr[led][i].push_back(_gPMTerr[led][i][s]);
	      /*if (led){
	      cout << " WORD " <<  gPD[led][i].back() << " " <<
		gPDerr[led][i].back() << " " <<
		gPMT[led][i].back() << " " <<
		gPMTerr[led][i].back() << endl;
	      //  }*/
	    }
	  }  
	  
	  //	  if (led) return 0;
	  
	  // Do fits

	  TString fit_string;
 	  
#if OFFSETPOLYFIT
	  //	  for (int constrainfit = 0; constrainfit < 2; constrainfit++){
	    // PREFIT to get approximate gain
	    TString prefit_string = "[0] + [1]*x + [2]*x*x";
	    //			  TF1 *prefit = new TF1("prefit", prefit_string, prefit_range_min[led], prefit_range_max[led]);
	    //	    TF1 *prefit = new TF1("prefit", prefit_string, range_min[led], range_max[led]);
	    TF1 *prefit = new TF1("prefit", prefit_string, range_min[led][i], range_max[led][i]);
	    //	    graph[led][i]->Fit(prefit, "R");  // save time, don't do fit
	    float prefitoffset = prefit->GetParameter(0);
	    float prefitgain = prefit->GetParameter(1);
	    if (!prefitgain) prefitgain = 1;
	    float PD_Bept = (BETAADCCOUNTS - prefitoffset)/prefitgain; // ignore quadratic term to get apprx PD value for Beta endpt
		//			  prefit->Delete();
	    
	    //   fit_string = "[0] + [1]*x + [2]*(x    )**2";
#if ORDER
	    /*    printf("Fitting LED linearity...\n");	
	    fit_string = "[0] + [1]*x";
	    for (int k = 1; k < ORDER; k++)
	      {
		fit_string += " + [";
		fit_string += k+1;
		fit_string += "]*x**";
		fit_string += k+1;
	      }

	    */
	    
	    fit_string = "[0] + [1]*x + [2]*(x)**2";

#endif

#if KEVSCALED 
	    fit_string = "[0] + [1]*(x-[3]) + [2]*(x-[3])**2";
#endif
#if FIXBETAENDPOINT 
	    fit_string = "[0] + [1]*(x-[3]) + [2]*(x-[3])**2";
	    //	    fit_string = "[0] + [1]*([4] + [5]*(x-[3]) + [6]*(x-[3])**2) + [2]*([4] + [5]*(x-[3]) + [6]*(x-[3])**2)**2";
#endif
	    cout << "Fit String: " <<  fit_string << endl;
	    // Carry out the fit
	    
	    cout << "Carrying out fit " << i << " for LED " << led << endl;
	    //TF1 * fit = new TF1("fit", fit_string, range_min[led], range_max[led]);
	    TF1 * fit = new TF1("fit", fit_string, range_min[led][i], range_max[led][i]);
	    // initialize fit function
	    
	    //Replaced with Energy Scaled graph 1/13/2015 
	    //	    if (led) fit->SetParameters(10, 0.5, -0.05); 
	    //	    else fit->SetParameters(10, 10.0, -0.05);
	    if (led) fit->SetParameters(750, 0.75, -0.0005, 10, 0.1, 0.0); 
	    else fit->SetParameters(0.0, 0.75, -0.0005, 10, 0.1, 0.0);
   
	    //	    float best_beta_endpt_PD = range_max[led]/beta_Bi_ratio; // make common to all loops
	    best_beta_endpt_PD[led].push_back(range_max[led][i]/beta_Bi_ratio);
#if RANGE_MAX_OVERRIDE // need actual bi peak position if range_max gets altered
	    // not using beta offsets for fits at the moment 4/30/15 SS
	    // best_beta_endpt_PD[led][i] = bi_peak_pd/beta_Bi_ratio;
#endif

	    //	    if (constrainfit){  // don't bother with this fit.
	    //	    cout << "Fixing x-offset to be apprx Beta endpoint " << PD_Bept << endl;
#if FIXBETAENDPOINT 
	    cout << "Fixing x-offset to be apprx Beta endpoint " << best_beta_endpt_PD[led][i] << endl;
	    if (led) fit->SetParameters(750., 5., -0.0005); 
	    //else fit->SetParameters(3000.0, 10., -0.1);
	    //	    else  fit->SetParameters(600, 10, -0.0005); 
	    else  fit->SetParameters(800, 20, -0.0005); 
	    fit->FixParameter(3, best_beta_endpt_PD[led][i]);  // fix x-offset to be the calculated PD value for beta endpoint 
#endif
	    // }

#if USEBETAOFFSETS
	    gPDBetaEP[led].push_back(best_beta_endpt_PD[led][i]); // push beta endpoints to global array for fit use	    
#else
	    gPDBetaEP[led].push_back(0);
#endif

#if DO_LED_FIT
#if KEVSCALED 
	    PMT_keV_graph[led][i]->Fit(fit, "R");
#else	    
	    graph[led][i]->Fit(fit, "R");  
#endif
	    //Testing
	    /*if (led>0)  {
	      PMT_keV_graph[0][i]->Draw("AP");
	      PMT_keV_graph[1][i]->Draw("Psame");
	      app.Run();
	      return 0;
	      }*/

	    std::cout << "A fit for run " << run << " finished." << std::endl;
	    
	    //Write fit parameters to file - include channel and led info
	    //	    fitpars = fit->GetParameters();
	    //fitparerrors = fit->GetParErrors();
	    fitpars[led] = fit->GetParameters();
	    fitparerrors[led] = fit->GetParErrors();
	    
	    fitchisq = fit->GetChisquare()/fit->GetNDF(); // use reduced chi-squared
	    
	    Npar = fit->GetNpar();
	    TString out_fit_string = "";
	    out_fit_string += run;                out_fit_string += "\t"; 
	    out_fit_string += i;                  out_fit_string += "\t"; 
	    out_fit_string += wavelength[led];    out_fit_string += "\t"; 
	    cout << out_fit_string << endl;
	    for (int k = 0; k < Npar; k++){
	      //	      out_fit_string += fitpars[k];       out_fit_string += "\t";
	      out_fit_string += fitpars[led][k];       out_fit_string += "\t";
	      //	      out_fit_string += fitparerrors[k];  out_fit_string += "\t";
	      out_fit_string += fitparerrors[led][k];  out_fit_string += "\t";
	      //	      cout << fitpars[k] << " +/- " << fitparerrors[k] << endl;
	      cout << fitpars[led][k] << " +/- " << fitparerrors[led][k] << endl;
	    }
	    out_fit_string += fitchisq;           out_fit_string += "\t";
	    out_fit_string += fit->GetNDF();      out_fit_string += "\t";
	    out_fit_string += range_min[led][i];  out_fit_string += "\t";
	    out_fit_string += range_max[led][i];  out_fit_string += "\t";
	    //	    out_fit_string += extended_range_max[led][i];     out_fit_string += "\t";
	    out_fit_string += "\n";
	    //if (!constrainfit){
	    fitfile << out_fit_string;
	    cout << "Saving unconstrained fit." << endl;
	    //}
	    /*if (constrainfit){
	      constrainfitfile << out_fit_string;
	      cout << "Saving constrained fit." << endl;
	      }
	    */
	    std::cout << "FITSTRING: " << endl <<  out_fit_string;
#endif	
#endif
	    // } // Don't do constrained fit
	    
	}

#if RELATIVEBETAPLOTS
      TF1 * relativefit = new TF1("relative", 
				  "[4] * ([0] + [1]*(x-[3]) + [2]*(x-[3])**2)",
				  range_min[UP][i], range_max[UP][i]);
      relativefit->FixParameter(0, fitpars[0][0]); // use 405 fitpars for 465 LED
      relativefit->FixParameter(1, fitpars[0][1]);
      relativefit->FixParameter(2, fitpars[0][2]);
      relativefit->FixParameter(3, fitpars[0][3]);
      relativefit->SetParameter(4, 0.3);
      graph[UP][i]->Fit(relativefit, "R");
      
      // Testing
      /*TCanvas *ppp = new TCanvas();
      graph[DOWN][i]->Draw("A*");
      graph[UP][i]->Draw("same*");
      app.Run();
      return 0;
      }*/
#endif	

#if MOVEDUPGRAPH      
      // Rescale 465 to 405 slope
      graph_465_moved_up[i] = new TGraphErrors(pulser_steps);
      /*     cout << "--------------------------------------------" << endl;
      cout << range_min[0]/pd_to_keV_factor[0] << " " << range_max[0]/pd_to_keV_factor[0] << endl;
      cout << range_min[1]/pd_to_keV_factor[1] << " " << range_max[1]/pd_to_keV_factor[1] << endl;
      cout << "PD/keV conversion factor: " << pd_to_keV_factor[0] << "; " << pd_to_keV_factor[1] << endl;
      cout << "--------------------------------------------" << endl;*/
      //TF1 * rescalefit = new TF1("rescalefit","pol1",range_min[0]/pd_to_keV_factor[0], range_max[0]/pd_to_keV_factor[0]);  // use range for 405 LED
      TF1 * rescalefit = new TF1("rescalefit","pol1",10.,100);  // use range to get a decent line
      rescalefit->SetParameters(0.0, 0.75);
      graph[UP][i]->Fit("rescalefit","R");
      graph[UP][i]->Fit("rescalefit","R");
      //      slope_UP = rescalefit->GetParameter(1);
      slope_UP = rescalefit->Eval(20.);
      rescalefit->SetParameters(0, 0.75); 
      graph[DOWN][i]->Fit("rescalefit","R");
      graph[DOWN][i]->Fit("rescalefit","R");
      //      slope_DOWN = rescalefit->GetParameter(1);
      slope_DOWN = rescalefit->Eval(20.);
      cout << "sloperatio = " << slope_DOWN/slope_UP << endl;
      for (int pulse = 0; pulse < pulser_steps; pulse++){
	double xUP,yUP;
	graph[UP][i]->GetPoint(pulse, xUP, yUP);
	double y_scaled = yUP*slope_DOWN/slope_UP;
	graph_465_moved_up[i]->SetPoint(pulse, xUP, y_scaled);
	//	cout << "yrat = " << y_scaled/yUP << endl; 
      }
      //Testing
      /*TCanvas *pop = new TCanvas();
      graph_465_moved_up[i]->Draw("A*");
      graph[DOWN][i]->Draw("*same");
      graph[UP][i]->Draw("*Same");
      rescalefit->Draw("same");
      app.Run();
      return 0;*/
#endif

  /// Plotting 

  printf("Plotting LED intensity...");	
  TString canvas_name = "canvas_";
  canvas_name += detector[i];
  TString canvas_title = "LED Scan ";
  canvas_title += detector[i];
  canvas[i] = new TCanvas(canvas_name, canvas_title);
  canvas[i]->Divide(2,1);
  canvas[i]->GetPad(1)->Divide(1,2);


  for (int led = 0; led < 2; led++)
    {
      canvas[i]->GetPad(1)->cd(led+1);
      
      TString title = "LED Scan ";
      title += detector[i];
      title += " (";
      title += wavelength[led];
      title += " nm)";
      pd_pmt_his2D[led][i]->SetTitle(title);
      pd_pmt_his2D[led][i]->GetZaxis()->SetRangeUser(0, 16);
      pd_pmt_his2D[led][i]->Draw("colz");
      
      // Draw what we fit
      graph[led][i]->SetMarkerColor(1);
      graph[led][i]->SetLineColor(1);
      graph[led][i]->SetMarkerStyle(21);
      graph[led][i]->SetMarkerSize(0.75);
      graph[led][i]->Draw("SameP");
    }
  printf("done.\n");	
  
  
  // split the canvas
  canvas[i]->GetPad(2)->Divide(1,2);
  
      // find the number of photoelectrons and plot
      printf("Building nPE plot...");	
      g[i] = new TGraph(pulser_steps);
      g_PE_PMT[i] = new TGraphErrors(pulser_steps);
      float d_pe = 0, d_pmt = 0;
      for (unsigned pulse = 0; pulse < pulser_steps; pulse++)
	{
	  double x, y;
	  graph[LED_TYPE][i]->GetPoint(pulse,x,y);
	  double sy = graph[LED_TYPE][i]->GetErrorY(pulse)*sqrt(num_pulses[LED_TYPE][pulse]);  //-1); sigma_mean = STDDEV/sqrt(N)
	  //	  std::cout << "PE: " << "y = " << y << " sy = " << sy << "\n";
	  float pe = (sy<1)? 0:pow(y/sy,2);
	  if (not isnan(pe) and y < max_adc_channel and pe < max_npe){
	    g[i]->SetPoint(pulse, x, pe);
	    g_PE_PMT[i]->SetPoint(pulse, pe, y);

	    // propogation of errors of pe, neglecting uncertainty in sy --> d(pe) = pe*2*d(y)/y = 2*d(y)*y/sy^2 =  2*y/(sqrt(N)*sy)
	    d_pe  = 2*y/( sy * sqrt(num_pulses[LED_TYPE][pulse]) );
	    d_pmt = graph[LED_TYPE][i]->GetErrorY(pulse);
	    g_PE_PMT[i]->SetPointError(pulse, d_pe, d_pmt);

	  }
        }
      
      // Find range in PE counts corresponding to ADC range
      printf("Finding PE range... ");
      float _max_pe = 0, _min_pe = 0;
      float range_max_pe = 0, range_min_pe = 0;

      for (int k = 0; k < pulser_steps; k++)
	{
	  double x,y;
	  // if (led){
	  //  g_PE_PMT[i]->GetPoint(k, x, y);
	  //}
	  //if (!led){ // 405 pulser counts backwards from large signal to small
	  
	  g_PE_PMT[i]->GetPoint(pulser_steps - k - 1, x, y);
	  
	  //}
	  //	      cout << x << "  " << y << endl;
	  if (y < ADC_min){_min_pe = x;}
	  if (y < ADC_max){_max_pe = x;}
	  if (y > ADC_max){break;}
	  
	  range_max_pe = _max_pe;
	  range_min_pe = _min_pe;
	}
      if (range_min_pe < 5.0) range_min_pe = 5.0;  // avoid bad clump of below-threshold LED shots
      printf("found PE range (%f,%f)\n", range_min_pe, range_max_pe); 
                 
      printf("Plotting number of photoelectrons...");	
      canvas[i]->GetPad(2)->cd(1);
      g[i]->SetTitle("Number of photoelectrons");
      g[i]->SetMarkerColor(2);
      g[i]->SetLineColor(1);
      g[i]->SetMarkerStyle(21);
      g[i]->SetMarkerSize(0.75);
      g[i]->Draw("AP");

      TString PE_PMT_title = "LED PMT-PE Curve: " + detector[i];
      PE_PMT_canvas->cd(i+1);
      g_PE_PMT[i]->SetTitle(PE_PMT_title.Data());
      g_PE_PMT[i]->SetMarkerColor(2);
      g_PE_PMT[i]->SetLineColor(1);
      g_PE_PMT[i]->SetMarkerStyle(21);
      g_PE_PMT[i]->SetMarkerSize(0.75);
      g_PE_PMT[i]->GetXaxis()->SetTitle("# PE");
      g_PE_PMT[i]->GetYaxis()->SetTitle("ADC Counts");
      g_PE_PMT[i]->Draw("AP");

#if DO_LED_FIT
      TF1 *pd_fit = new TF1("polyfit", "[0]*x", range_min[LED_TYPE][i], range_max[LED_TYPE][i]);
      if (g[i]->Fit(pd_fit, "R"))
	continue;
      printf("done.\n");

      TF1 *pe_fit = new TF1("pefit", "[0] + [1]*x + [2]*x*x", range_min_pe, range_max_pe);
      if (g_PE_PMT[i]->Fit(pe_fit, "R"))
	continue;
      printf("PE fit done.\n");
 #endif
      
      printf("Plotting residuals...");	
      canvas[i]->GetPad(2)->cd(2);
      resg[i] = new TGraphErrors(*graph[LED_TYPE][i]);
      resg[i]->SetTitle("PMT Linearity Residual");
      resg[i]->SetMinimum(-0.1);
      resg[i]->SetMaximum(0.1);
      resg[i]->Draw("AP");
      //g[i]->SetLineColor(4);
      printf("done.\n");
      
      // Write pe_pmt fits to file
      pefitpars = pe_fit->GetParameters();
      pefitparerrors = pe_fit->GetParErrors();
      pefitchisq = pe_fit->GetChisquare()/pe_fit->GetNDF(); // use reduced chi-squared
      
      peNpar = pe_fit->GetNpar();
      TString pe_out_fit_string = "";
      pe_out_fit_string += run;                pe_out_fit_string += "\t"; 
      pe_out_fit_string += i;                  pe_out_fit_string += "\t"; 
      pe_out_fit_string += "405";              pe_out_fit_string += "\t"; 
      //      cout << pe_out_fit_string << endl;
      for (int k = 0; k < Npar; k++){
	pe_out_fit_string += pefitpars[k];       pe_out_fit_string += "\t";
	pe_out_fit_string += pefitparerrors[k];  pe_out_fit_string += "\t";
	cout << pefitpars[k] << " +/- " << pefitparerrors[k] << endl;
      }
      pe_out_fit_string += pefitchisq;         pe_out_fit_string += "\t";
      pe_out_fit_string += pe_fit->GetNDF();      pe_out_fit_string += "\t";
      cout << pefitchisq << endl;
      pe_out_fit_string += "\n";
      pepmtfitfile << pe_out_fit_string;

    }  // end loop over channels
  
  //  const int nvars = 27;
// const int nvars = 28; 
  const int nvars = 19; 

  // Do Minuit stuff 
  TMinuit *gMinuit = new TMinuit(nvars);  //initialize TMinuit with a maximum of nvars params
  gMinuit->SetPrintLevel(0); // 0 = normal, 1 = verbose
  gMinuit->SetFCN(fcn);

  Double_t arglist[10];
  Int_t ierflg = 0;
  
  arglist[0] = 1;
  gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
  
  // Set starting values and step sizes for parameters
  //led = 0
  /*  static Double_t vstart[nvars] = {1000., 10., -0.01,
				   1000., 10., -0.01,
				   1000., 10., -0.01,
				   1000., 10., -0.01,
				   1000., 10., -0.01,
				   1000., 10., -0.01,
				   1000., 10., -0.01,
				   1000., 10., -0.01,
				   0., 1., 0., 
				   3};
  
  static Double_t step[nvars] = {10 , 1  , 0.01,
				 10, 1, 0.01,
				 10 ,1 , 0.01,
				 10 ,1 , 0.01,
				 10 ,1 , 0.01,
				 10 ,1 , 0.01,
				 10 ,1 , 0.01,
				 10 ,1 , 0.01, 
				 1. ,1.  , 0.0001, 
				 1.};*/

  /*static Double_t vstart[nvars] = {100., 1., -0.01,
				   100., 1., -0.01,
				   100., 1., -0.01,
				   100., 1., -0.01,
				   100., 1., -0.01,
				   100., 1., -0.01,
				   100., 1., -0.01,
				   100., 1., -0.01,
				   100., 1., 0., 
				   3};
  
  static Double_t step[nvars] = {1 , 0.1  , 0.1,
				 1, 0.1, 0.1,
				 1 ,0.1 , 0.1,
				 1 ,0.1 , 0.1,
				 1 ,0.1 , 0.1,
				 1 ,0.1 , 0.1,
				 1 ,0.1 , 0.1,
				 1 ,0.1 , 0.1, 
				 1. ,0.1 , 0.1, 
				 0.1};*/

  /*  static Double_t vstart[nvars] = {10., 3., -0.0001,
				   10., 3., -0.0001,
				   10., 3., -0.0001,
				   10., 3., -0.0001,
				   10., 10., -0.0001,
				   10., 10., -0.0001,
				   10., 10., -0.0001,
				   10., 10., -0.0001,
				   10., 3., 0., 
				   3};
 
  static Double_t step[nvars] = {1 , 0.1  , 0.00001,
				 1, 0.1, 0.00001,
				 1 ,0.1 , 0.00001,
				 1 ,0.1 , 0.00001,
				 1 ,0.1 , 0.00001,
				 1 ,0.1 , 0.00001,
				 1 ,0.1 , 0.00001,
				 1 ,0.1 , 0.00001, 
				 1. ,0.1 , 0.00001, 
				 0.1};*/

  static Double_t vstart[nvars] = {10., 1.,      10., 1.,
				   10., 1.,      10., 1.,
				   10., 1.,      10., 1.,
				   10., 1.,      10., 1.,
				   3., -0.001, 5.};
 
  static Double_t step[nvars] = {1., 0.1,        1., 0.1, 
				 1., 0.1,        1., 0.1, 
				 1., 0.1,        1., 0.1, 
				 1., 0.1,        1., 0.1, 
				 0.1, 0.00001, 0.1};

  /*  //led = 1 
  static Double_t vstart_465[nvars] = {1000., 0.5, -0.001,
				       1000., 0.5, -0.001,
				       1000., 0.5, -0.001,
				       1000., 0.5, -0.001,
				       1000., 0.5, -0.001,
				       1000., 0.5, -0.001,
				       1000., 0.5, -0.001,
				       1000., 0.5, -0.001,
				       0., 1., 0., 
				       3};
  
  static Double_t step_465[nvars] = {10 , 0.1  , 0.001,
				     10, 0.1, 0.001,
				     10 ,0.1 , 0.001,
				     10 ,0.1 , 0.001,
				     10 ,0.1 , 0.001,
				     10 ,0.1 , 0.001,
				     10 ,0.1 , 0.001,
				     10 ,0.1 , 0.001, 
				     1. ,0.1  , 0.0001, 
				     1.};*/
  //  int led = COMBINEDLED;  

  double PD_parms[3], PD_errs[3];
  TF1 fittedFunctions[2][NUM_CHANNELS]; // separate functions for each LED
  //  TF1 LEDFreeFuncs[2][NUM_CHANNELS]; // not as interesting 

  //  for (int led = 0; led < 2; led++){
    //  cout << "----------------------------------------------------" << endl;
  //  cout << "-------  MINUIT FOR LED " << led << "  -------------------------" << endl;
  //  cout << "----------------------------------------------------" << endl;

    //    gLED = led; // fill the global parameter telling us which LED we're on

  //    if (led) { // use correct initial parameters.
  //    memcpy(vstart, vstart_465, sizeof(vstart));
  //   memcpy(step, step_465, sizeof(step));
  //  }   
    
    //  gMinuit->FixParameter(27);
    
    /*    gMinuit->mnparm(0, "t0p0", vstart[0], step[0], 0,0,ierflg);
    gMinuit->mnparm(1, "t0p1", vstart[1], step[1], 0,0,ierflg);
    gMinuit->mnparm(2, "t0p2", vstart[2], step[2], 0,0,ierflg);
    gMinuit->mnparm(3, "t1p0", vstart[3], step[3], 0,0,ierflg);
    gMinuit->mnparm(4, "t1p1", vstart[4], step[4], 0,0,ierflg);
    gMinuit->mnparm(5, "t1p2", vstart[5], step[5], 0,0,ierflg);
    gMinuit->mnparm(6, "t2p0", vstart[6], step[6], 0,0,ierflg);
    gMinuit->mnparm(7, "t2p1", vstart[7], step[7], 0,0,ierflg);
    gMinuit->mnparm(8, "t2p2", vstart[8], step[8], 0,0,ierflg);
    gMinuit->mnparm(9, "t3p0", vstart[9], step[9], 0,0,ierflg);
    gMinuit->mnparm(10, "t3p1", vstart[10], step[10], 0,0,ierflg);
    gMinuit->mnparm(11, "t3p2", vstart[11], step[11], 0,0,ierflg);
    gMinuit->mnparm(12, "t4p0", vstart[12], step[12], 0,0,ierflg);
    gMinuit->mnparm(13, "t4p1", vstart[13], step[13], 0,0,ierflg);
    gMinuit->mnparm(14, "t4p2", vstart[14], step[14], 0,0,ierflg);
    gMinuit->mnparm(15, "t5p0", vstart[15], step[15], 0,0,ierflg);
    gMinuit->mnparm(16, "t5p1", vstart[16], step[16], 0,0,ierflg);
    gMinuit->mnparm(17, "t5p2", vstart[17], step[17], 0,0,ierflg);
    gMinuit->mnparm(18, "t6p0", vstart[18], step[18], 0,0,ierflg);
    gMinuit->mnparm(19, "t6p1", vstart[19], step[19], 0,0,ierflg);
    gMinuit->mnparm(20, "t6p2", vstart[20], step[20], 0,0,ierflg);
    gMinuit->mnparm(21, "t7p0", vstart[21], step[21], 0,0,ierflg);
    gMinuit->mnparm(22, "t7p1", vstart[22], step[22], 0,0,ierflg);
    gMinuit->mnparm(23, "t7p2", vstart[23], step[23], 0,0,ierflg);
    gMinuit->mnparm(24, "PDp0", vstart[24], step[24], 0,0,ierflg);
    gMinuit->mnparm(25, "PDp1", vstart[25], step[25], 0,0,ierflg);
    gMinuit->mnparm(26, "PDp2", vstart[26], step[26], 0,0,ierflg);
    gMinuit->mnparm(27, "PDratio", vstart[27], step[27], 0, 0, ierflg);*/
   
    //  gMinuit->mnparm(27, "PD_Beta", vstart[27], step[27], vstart[27], vstart[27], ierflg);
   
    for (int pp = 0; pp < nvars-3; pp++){
      gMinuit->mnparm(pp, Form("p%i", pp%2), vstart[pp], step[pp], 0, 0, ierflg);
      //      gMinuit->mnparm(pp, Form("p%i", pp%2), vstart[pp], step[pp], -100., 1000., ierflg);
    }
    gMinuit->mnparm(nvars-3, "PDp1", vstart[nvars-3], step[nvars-3], 0.,10.,ierflg);
    gMinuit->mnparm(nvars-2, "PDp2", vstart[nvars-2], step[nvars-2], -1.,1.,ierflg);
    gMinuit->mnparm(nvars-1, "PDratio", vstart[nvars-1], step[nvars-1], 0., 10., ierflg);
 
    // Now ready for minimization step
    arglist[0] = 30000;
    arglist[1] = 0.1;
    gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
    
    // Print results
    Double_t amin,edm,errdef;
    Int_t nvpar,nparx,icstat;
    gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
    gMinuit->mnprin(3,amin);
    
    // write to file
    double pp, pperr; 
    TString comb_fit_string = "";
    //    for (int i = 0; i < NUM_CHANNELS + 1; i++){
    // for (int j = 0; j < 3; j++){
    for (int p = 0; p < nvars; p++){
      //	int p = 2*i + j;
      gMinuit->GetParameter(p, pp, pperr);
      comb_fit_string += run;                comb_fit_string += "\t"; 
      int i = p/2; 
      comb_fit_string += i;                  comb_fit_string += "\t"; 
      //	comb_fit_string += led;                comb_fit_string += "\t"; 
      comb_fit_string += gMinuit->fCpnam[p]; comb_fit_string += "\t";
      comb_fit_string += pp;                 comb_fit_string += "\t";
      comb_fit_string += pperr;              comb_fit_string += "\t";
      //	comb_fit_string += pefitchisq;         comb_fit_string += "\t";
      //	comb_fit_string += pe_fit->GetNDF();      comb_fit_string += "\t";
      comb_fit_string += range_min[0][i];  comb_fit_string += "\t";
      comb_fit_string += range_max[0][i];  comb_fit_string += "\t";
      comb_fit_string += range_min[1][i];  comb_fit_string += "\t";
      comb_fit_string += range_max[1][i];  comb_fit_string += "\t";
      comb_fit_string += "\n";
    }
    // find number of points used in fit
    Double_t NFit = 0;
    for (int led = 0; led < 2; led++){
      for (int i = 0; i < NUM_CHANNELS; i++){
	NFit += gPD[led][i].size();
      }
    }
    comb_fit_string += run;                  comb_fit_string += "\t"; 
    comb_fit_string += "10";                 comb_fit_string += "\t"; 
    comb_fit_string += "Chisq/nvars";        comb_fit_string += "\t"; 
    comb_fit_string += amin;                 comb_fit_string +="\t";
    comb_fit_string += nvpar;                comb_fit_string +="\t";
    comb_fit_string += NFit;                 comb_fit_string +="\t";
    comb_fit_string += icstat;               comb_fit_string +="\n";
    combfitfile << comb_fit_string;
	
    //Define functions from fitted parms.
      
    //    double PD_parms[3], PD_errs[3];
    gMinuit->GetParameter(nvars-3, PD_parms[0], PD_errs[0]);
    gMinuit->GetParameter(nvars-2, PD_parms[1], PD_errs[1]);
    gMinuit->GetParameter(nvars-1, PD_parms[2], PD_errs[2]);
    
      
    double p_val[2], p_err[2];
    double PDratio, PDratioErr;
    for (int led = 0; led < 2; led++){
      for (int i = 0; i < NUM_CHANNELS; i++){
	fittedFunctions[led][i] = TF1(Form("f_%i", i), "(1./[1]) * (-[0] + [4]*([2]*x + [3]*x*x) )", 
				      RANGE_MIN, range_max[led][i]); 
	//			      RANGE_MIN, extended_range_max[led][i]); // deprecated

     	for (int j = 0; j < 2; j++){
	  gMinuit->GetParameter(2*i + j, p_val[j], p_err[j]);
	  fittedFunctions[led][i].SetParameter(j, p_val[j]);
	  fittedFunctions[led][i].SetParameter(j+2, PD_parms[j]);
	  if (led == 0) fittedFunctions[led][i].SetParameter(4, PD_parms[2]);
	  if (led == 1) fittedFunctions[led][i].SetParameter(4, 1.0);
	}
      }
    }
    //fittedFunctions[led][i].SetParameter(6, gPDoff[led]);
    //fittedFunctions[led][i].SetParameter(7, gPDBetaEP[led][i]);
    //	fittedFunctions[led][i].SetParameter(6, gPDBetaEP[led][i]);
    //fittedFunctions[led][i].SetParameter(7, gPMTBetaEP[i]);
    
    // if (led == 0) fittedFunctions[led][i].SetParameter(8, PD_parms[2]);
    //if (led == 1) fittedFunctions[led][i].SetParameter(8, 1);
      
      /*     for (int p = 0; p < 9; p++){
	     cout << led << " " << i << " " << p << " " << fittedFunctions[led][i].GetParameter(p) << endl;
	     } */



    //Testing
    /*  TCanvas * c1 = new TCanvas;
	fittedFunctions[0].Draw();
	c1->SaveAs("test.gif");
	
	return 0;
    */
    /*#if USE_ROOT_APPLICATION
  // run the root application
  app.Run();
  #endif
  return 0;
  */

  fitfile.close();
  pepmtfitfile.close();
  convfactorfile.close();
  combfitfile.close();
  //  constrainfitfile.close();

   TString pe_pmt_filename = "pe_pmt";
   pe_pmt_filename = OUTPUT_IMAGE_DIR + pe_pmt_filename;
   pe_pmt_filename += argv[1];
   TString pe_pmt_rootfilename = pe_pmt_filename;
   pe_pmt_filename += ".gif";
   pe_pmt_rootfilename += ".root";
   PE_PMT_canvas->SaveAs(pe_pmt_filename,"9");
   PE_PMT_canvas->SaveAs(pe_pmt_rootfilename,"9");

  
  std::cout << "Fitting Gain Parameters" << std::endl;
  TF1 * mygaus = new TF1("mygaus", "gaus", 0, 1000);
  time_his1D->Fit("mygaus");
  double * projfit_pars = mygaus->GetParameters();
  double * projfit_errs = mygaus->GetParErrors();
  int projfit_Npars = mygaus->GetNpar();
  double projfit_chisq = mygaus->GetChisquare()/mygaus->GetNDF();
  
  std::cout << "Writing gain parameters" << std::endl;
  TString gainfilename = "GainResults.txt";
  gainfilename = OUTPUT_IMAGE_DIR + gainfilename;
  ofstream gainfile;
  gainfile.open(gainfilename, std::ofstream::out | std::ofstream::app);
  TString out_gain_string = "";
  out_gain_string += run;                        out_gain_string += "\t"; 
  for (int q = 0; q < projfit_Npars; q++){
    out_gain_string += projfit_pars[q];          out_gain_string += "\t";
    out_gain_string += projfit_errs[q];          out_gain_string += "\t";
  }
  out_gain_string += projfit_chisq;              out_gain_string += "\t";
  out_gain_string += "\n";
  gainfile << out_gain_string;	
  
  gainfile.close();
  
  // Fit and write out results for PMT response during gain period
  double * pmt_pars;
  double * pmt_errs;
  int pmt_Npars = 0;
  double pmt_chisq = 0.;
  TF1 * pmtgaus = new TF1("pmtgaus", "gaus", 0, 3000);
  
  TString PMTgainfilename = "PMTGainResults.txt";
  PMTgainfilename = OUTPUT_IMAGE_DIR + PMTgainfilename;
  ofstream PMTgainfile;
  PMTgainfile.open(PMTgainfilename, std::ofstream::out | std::ofstream::app);
  
  for (int chan = 0; chan < NUM_CHANNELS; chan++){
    std::cout << "Fitting PMT gain parameters " << chan << std::endl;
    pmt_gain_his1D[chan]->Fit("pmtgaus","R");
    pmt_pars = pmtgaus->GetParameters();
    pmt_errs = pmtgaus->GetParErrors();
    pmt_Npars = pmtgaus->GetNpar();
    pmt_chisq = pmtgaus->GetChisquare()/pmtgaus->GetNDF();
    
    std::cout << "Writing ADC gain parameters " << chan << std::endl;
    
    TString adc_out_gain_string = "";
    adc_out_gain_string += run;                    adc_out_gain_string += "\t"; 
    adc_out_gain_string += chan;                   adc_out_gain_string += "\t"; 
    for (int q = 0; q < pmt_Npars; q++){
      adc_out_gain_string += pmt_pars[q];          adc_out_gain_string += "\t";
      adc_out_gain_string += pmt_errs[q];          adc_out_gain_string += "\t";
    }
    adc_out_gain_string += pmt_chisq;              adc_out_gain_string += "\t";
    adc_out_gain_string += "\n";
    PMTgainfile << adc_out_gain_string;	
  }
  PMTgainfile.close();
  
  // draw time sequence
  TCanvas *time_canvas = new TCanvas("time_canvas", "PD time sequence");
  time_canvas->Divide(2,1);
  time_canvas->cd(1);
  //	true_time_his2D->GetYaxis()->SetRangeUser(-pd_pedestal, 500-pd_pedestal);
  true_time_his2D->Draw("scat");
  time_canvas->cd(2);
  time_his2D[TIME]->GetZaxis()->SetRangeUser(0, 4);
  time_his2D[TIME]->Draw("scat");
  
#if OUTPUT_IMAGE
  TString pd_time_filename = "pd_time_";
  pd_time_filename = OUTPUT_IMAGE_DIR + pd_time_filename;
  pd_time_filename += argv[1];
  //		TString pd_time_rootfilename = pd_time_filename;
  pd_time_filename += ".gif";
  //	pd_time_rootfilename += ".root";
  time_canvas->SaveAs(pd_time_filename,"9");
  //		time_canvas->SaveAs(pd_time_rootfilename,"9");
#endif
  
  TCanvas *colored_time_canvas = new TCanvas("colored_time_canvas", "Cut PD time sequence");
  time_his2D[UP]->SetMarkerColor(1);
  time_his2D[UP]->Draw("scat");
  time_his2D[DOWN]->SetMarkerColor(2);
  time_his2D[DOWN]->Draw("samescat");
  time_his2D[GAIN]->SetMarkerColor(3);
  time_his2D[GAIN]->Draw("samescat");
    
#if OUTPUT_IMAGE
  TString cut_pd_time_filename = "cut_pd_time_";
  cut_pd_time_filename = OUTPUT_IMAGE_DIR + cut_pd_time_filename;
  cut_pd_time_filename += argv[1];
  TString cut_pd_time_rootfilename = cut_pd_time_filename;
  cut_pd_time_filename += ".gif";
  cut_pd_time_rootfilename += ".root";
  colored_time_canvas->SaveAs(cut_pd_time_filename,"9");
  colored_time_canvas->SaveAs(cut_pd_time_rootfilename,"9");
#endif
  
#if NUM_CHANNELS == 8
  TCanvas *ew_canvas = new TCanvas("all_pmt_linearity_canvas", 
				   "PMT linearity scans for all tubes", 1280, 720);
  TCanvas *ew_canvas_scaled = new TCanvas("all_pmt_linearity_canvas_scaled", 
				   "PMT linearity scans for all tubes", 1280, 720);
  ew_canvas->Divide(2,1);
  ew_canvas_scaled->Divide(2,1);
  for (int ew = 0; ew < 2; ew++) {
    ew_canvas->GetPad(ew+1)->Divide(2,2);
    ew_canvas_scaled->GetPad(ew+1)->Divide(2,2);
    for (int tx = 0; tx < 2; tx++) {
      for (int ty = 0; ty < 2; ty++) {
	int i = 4*ew+tx+2*ty;
	ew_canvas->GetPad(ew+1)->cd(tx+2*ty+1);
	TString title = "Run "; 
	title += run;
	title += "    LED Scan: ";
	title += detector[i];
	title += " (";
	title += wavelength[DOWN];
#if PLOTBOTHLEDS
	title += " nm and ";
	title += wavelength[UP];
#endif
	title += " nm)";
	//pd_pmt_his2D[UP][i]->SetTitle(title);
	//pd_pmt_his2D[DOWN][i]->SetTitle(title);
	//pd_pmt_his2D[DOWN][i]->Draw("colz");
	//graph[DOWN][i]->Draw("SameP");
	gStyle->SetOptFit(0);
	gStyle->SetOptStat(0);
	graph[DOWN][i]->SetMarkerColor(6);
	graph[UP][i]->Draw("AP");
	graph[UP][i]->SetName("GraphUP");
	//	graph[DOWN][i]->Draw("AP");
	graph[DOWN][i]->Draw("sameP");
	graph[DOWN][i]->SetName("GraphDOWN");
	graph[UP][i]->SetMarkerColor(4);
	//graph[UP][i]->Draw("SameP");
	//	graph[UP][i]->SetName("GraphUP");
	/*	graph[DOWN][i]->SetTitle(title);
	graph[DOWN][i]->GetXaxis()->SetTitle("PD");
	graph[DOWN][i]->GetYaxis()->SetTitle("PMT (ADC)");*/
	graph[UP][i]->SetTitle(title);
	graph[UP][i]->GetXaxis()->SetTitle("PD");
	graph[UP][i]->GetYaxis()->SetTitle("PMT (ADC)");
	graph[UP][i]->GetYaxis()->SetRangeUser(-50, 4000);

	for (int led = 0; led < 2; led++){
	  fittedFunctions[led][i].SetLineColor(8);
	  fittedFunctions[led][i].Draw("same");
	}

	TPad *subpad = new TPad("subpad","",0.37,0.58,0.7,0.89);
	subpad->Draw();
	subpad->cd();
	TGraph * newgraphUP = (TGraph*)graph[UP][i]->Clone("newUP");
	TGraph * newgraphDOWN = (TGraph*)graph[DOWN][i]->Clone("newDOWN");
	newgraphUP->Draw("AP");
	newgraphDOWN->Draw("sameP");
	newgraphUP->SetTitle("");
	newgraphUP->GetXaxis()->SetTitle("");
	newgraphUP->GetYaxis()->SetRangeUser(-50, 750);
	newgraphUP->GetXaxis()->SetRangeUser(-10, 50);

	for (int led = 0; led < 2; led++){
	  fittedFunctions[led][i].SetLineColor(8);
	  fittedFunctions[led][i].Draw("same");
	}	// redraw on inset

	ew_canvas->Update();
	TPaveStats * st = (TPaveStats*)graph[DOWN][i]->GetListOfFunctions()->FindObject("stats");
	if (st==NULL) cout << "Null TPaveStats - oops?" << endl;
	else {
	  st->SetX1NDC(0.1);
	  st->SetX2NDC(0.5);
	}

#if MOVEDUPGRAPH
	graph_465_moved_up[i]->SetMarkerColor(4);
	graph_465_moved_up[i]->SetMarkerStyle(24);
	graph_465_moved_up[i]->Draw("sameP");
#endif
	ew_canvas_scaled->GetPad(ew+1)->cd(tx+2*ty+1);
	PMT_keV_graph[DOWN][i]->SetTitle(title);
	PMT_keV_graph[DOWN][i]->GetXaxis()->SetTitle("PD (pseudo-keV)");
	PMT_keV_graph[DOWN][i]->GetYaxis()->SetTitle("PMT (ADC)");
	PMT_keV_graph[DOWN][i]->Draw("AP");
		
#if PLOTBOTHLEDS
	//pd_pmt_his2D[UP][i]->Draw("Samecolz");
	PMT_keV_graph[UP][i]->Draw("sameP");
	PMT_keV_graph[UP][i]->SetName("ScaledGraphUP");
#endif
	//	pd_pmt_his2D[DOWN][i]->GetYaxis()->SetRangeUser(0, 2000);
      }
    }
  }

  /* // don't bother, not doing  
  // Repeat for constrained fits
  TCanvas *ew_constr_canvas = new TCanvas("all_pmt_linearity_constrained_canvas", 
					  "PMT linearity scans for all tubes", 1280, 720);
  ew_constr_canvas->Divide(2,1);
  for (int ew = 0; ew < 2; ew++) {
    ew_constr_canvas->GetPad(ew+1)->Divide(2,2);
    for (int tx = 0; tx < 2; tx++) {
      for (int ty = 0; ty < 2; ty++) {
	int i = 4*ew+tx+2*ty;
	ew_constr_canvas->GetPad(ew+1)->cd(tx+2*ty+1);
	TString title = "Run "; 
	title += run;
	title += "    LED Scan: ";
	title += detector[i];
	title += " (";
	title += wavelength[DOWN];
	title += " nm and ";
	title += wavelength[UP];
	title += " nm)";
	pd_pmt_his2D[UP][i]->SetTitle(title);
	pd_pmt_his2D[DOWN][i]->SetTitle(title);
	pd_pmt_his2D[UP][i]->Draw("colz");
	constrained_graph[UP][i]->Draw("SameP");
	constrained_graph[UP][i]->SetName("GraphUP");
	pd_pmt_his2D[DOWN][i]->Draw("Samecolz");
	constrained_graph[DOWN][i]->Draw("SameP");
	constrained_graph[DOWN][i]->SetName("GraphDOWN");
      }
    }
  }
  */
#if OUTPUT_IMAGE
  TString pd_led_pmt_filename = "pd_led_pmt_";
  pd_led_pmt_filename = OUTPUT_IMAGE_DIR + pd_led_pmt_filename;
  //    TString pd_led_pmt_filename = "images/pd_led_pmt_";
  pd_led_pmt_filename += argv[1];
  TString pd_led_pmt_constr_filename = pd_led_pmt_filename + "_constrained.gif";
  TString pd_led_pmt_constr_rootfilename = pd_led_pmt_filename + "_constrained.root";

#if RANGE_MAX_OVERRIDE
  char c[32];
  sprintf(c,"%.0f",RANGE_MAX_VALUE);
  pd_led_pmt_filename += "_";
  pd_led_pmt_filename += std::string(c); 
#endif

  TString pd_led_pmt_rootfilename = pd_led_pmt_filename;
  pd_led_pmt_filename += ".gif";
  pd_led_pmt_rootfilename += ".root";
  TString pd_led_pmt_scaled_filename = OUTPUT_IMAGE_DIR;
  pd_led_pmt_scaled_filename +=  "pd_led_pmt_scaled_";
  pd_led_pmt_scaled_filename +=  argv[1];
  pd_led_pmt_scaled_filename += ".gif";
  ew_canvas->SaveAs(pd_led_pmt_filename,"9");
  // ew_constr_canvas->SaveAs(pd_led_pmt_constr_filename, "9");
  ew_canvas->SaveAs(pd_led_pmt_rootfilename,"9");
  //ew_constr_canvas->SaveAs(pd_led_pmt_constr_rootfilename, "9");
  ew_canvas_scaled->SaveAs(pd_led_pmt_scaled_filename, "9");

#endif
#endif
  
  TCanvas * gain_canvas = new TCanvas("gain_canvas", "Fitted Gain", 1000, 800);
  gain_canvas->Divide(1,2);
  gain_canvas->cd(1);
  time_his2D[GAIN]->Draw("scat");
  //   gain_his_pol0->Draw("scat");
  gain_canvas->cd(2);
  //    gain_his_pol1->Draw("scat");
  //gain_canvas->cd(3); 
  time_his1D->Draw();
  
  
#if OUTPUT_IMAGE
  TString pd_gain_filename = "pd_gain_";
  pd_gain_filename = OUTPUT_IMAGE_DIR + pd_gain_filename;
  pd_gain_filename += argv[1];
  TString pd_gain_rootfilename = pd_gain_filename;
  pd_gain_filename += ".gif";
  //	   pd_gain_rootfilename += ".root";
  gain_canvas->SaveAs(pd_gain_filename, "9");
  //	   gain_canvas->SaveAs(pd_gain_rootfilename, "9");
#endif
  
  TCanvas * pmt_gain_canvas = new TCanvas("pmt_gain_canvas", "Fitted ADC gains", 1000, 800);
  pmt_gain_canvas->Divide(4,2);
  for (int chpmt = 0; chpmt < NUM_CHANNELS; chpmt++){
    pmt_gain_canvas->cd(chpmt+1);
    pmt_gain_his1D[chpmt]->Draw();
  }
  TCanvas * pmt_time_canvas = new TCanvas("pmt_time", "", 1000, 800);
  pmt_time_canvas->Divide(4,2);
  for (int chpmt = 0; chpmt < NUM_CHANNELS; chpmt++){
    pmt_time_canvas->cd(chpmt+1);
    pmt_gain_his2D[chpmt]->SetMaximum(5);
    pmt_gain_his2D[chpmt]->Draw("colz");
  }
  
#if OUTPUT_IMAGE
  TString pmt_gain_filename = "pmt_gain_";
  pmt_gain_filename = OUTPUT_IMAGE_DIR + pmt_gain_filename;
  pmt_gain_filename += argv[1];
  TString pmt_gain_rootfilename = pmt_gain_filename;
  pmt_gain_filename += ".gif";
  //	   pmt_gain_rootfilename += ".root";
  pmt_gain_canvas->SaveAs(pmt_gain_filename, "9");
  //	   pmt_gain_canvas->SaveAs(pmt_gain_rootfilename, "9");
  
  TString pmt_time_filename = "pmt_time_";
  pmt_time_filename = OUTPUT_IMAGE_DIR + pmt_time_filename;
  pmt_time_filename += argv[1];
  TString pmt_time_rootfilename = pmt_time_filename;
  pmt_time_filename += ".gif";
  //	   pmt_time_rootfilename += ".root";
  //  pmt_time_canvas->SaveAs(pmt_time_filename, "9");
  //	   pmt_time_canvas->SaveAs(pmt_time_rootfilename, "9");
  
#endif
  
#if USE_ROOT_APPLICATION
  // run the root application
  app.Run();
#endif
  
  return 0;
}
