// Quick Script to fit the Bi peak in the AnodeCut Spectra
// Simon Slutsky 05/28/14

// Start a root session and run 
// .L QuickBiPeakStudy.cc

// Compile with:
// g++ -o QuickBiPeakStudy QuickBiPeakStudy.cc `root-config --cflags --glibs`

#include <iostream>
#include <TH1.h>
#include <TF1.h>
#include <TFile.h>
#include <TCanvas.h>
#include <math.h>
#include <TStyle.h>
#include <stdio.h>
#include <string>

using namespace std;

//double gaus_fit_routine(TH1F* htemp){
TF1 * gaus_fit_routine(TH1F* htemp){
  TF1 * gausfit = new TF1("gausfit", "gaus(0)",0, 3000);
  gausfit->SetParameters(300, 1000, 100); 
  htemp->Fit("gausfit", "Q"); // Fit 0
  
  double  fitmean = gausfit->GetParameter(1);
  double  fitsig = fabs(gausfit->GetParameter(2));
  // asymmetric fractions push fit away from secondary peak at 600 kev
  htemp->Fit("gausfit", "RQ", "", fitmean - 1.0*fitsig, fitmean + 2.0*fitsig);   // Fit 1
  
  fitmean = gausfit->GetParameter(1);
  fitsig = fabs(gausfit->GetParameter(2));
  htemp->Fit("gausfit", "RQ", "", fitmean - 1.0*fitsig, fitmean + 2.0*fitsig);  // Fit 2
  
  htemp->Draw();

  // return gausfit->GetParameter(1);
  return gausfit;
}

double find_half_height(TF1 * fitgaus){
  // For gaussian f(x) = exp[ -(x-u)^2/2s^2 ], 
  // half_height occurs when x = s*sqrt(2*log(2)) + u
  
  double mu = fitgaus->GetParameter(1);
  double sigma = fitgaus->GetParameter(2);
  double factor = sqrt(2*log(2));

  double x_halfheight = sigma*factor + mu;
  return x_halfheight;
}

double find_hh_peak_ratio(double half_height, double peak){
  ratio = half_height/peak;
  return ratio;
}

vector<double> fit_bi_peaks(char* filename){
  //cout << "Recieved file " << filename << endl;
  string filestring(filename);
  string run_num = filestring.substr(34,5);
  // ------------------------------------------------
  cout << run_num << "\t"; 
  // ------------------------------------------------
  TFile *MyFile = new TFile(filename,"READ");
  TCanvas * c1 = new TCanvas("c1", "c1");
  c1->Divide(2,4);
  gStyle->SetOptFit(1111);
  vector<double> mean_vector;
  vector<double> mean_errs_vector;
  
  int nTubes = 8;
  for (int i = 0; i < nTubes; i++){
    TH1F * htemp  = (TH1F*)MyFile->Get( Form("h190%i", i) );
    c1->cd(i+1);
    //    double mean_i = gaus_fit_routine(htemp);
    TF1* fitfunc = gaus_fit_routine(htemp);
    mean_vector.push_back(fitfunc->GetParameter(1));
    mean_errs_vector.push_back(fitfunc->GetParError(1) );
    // ------------------------------------------------
    cout << mean_vector[i] << "\t";
    cout << mean_errs_vector[i] << "\t";
    // ------------------------------------------------
    if (i < nTubes - 1){
      // ------------------------------------------------
      cout << "\t";
      // ------------------------------------------------
    }
  }
  // ------------------------------------------------
  cout << "\n";
  // ------------------------------------------------

  c1->SaveAs( Form("./Figures/BiSourcePeaks_%s.pdf", run_num.c_str()) ); 
  return mean_vector;

}

int main(int argc, char *argv[]){
  if (argc==1){
    cout << "No filename passed. Exiting" << endl;
    //    return -1;
  }  
  char * MyFileName = argv[1];
  //cout << MyFileName << endl;
  
  //cout << "Processing File" << MyFileName << endl;
  fit_bi_peaks(MyFileName);
 
  return 0;
}

/* Testing
  TF1 * mygaus = new TF1("mygaus", "gaus(0)", -1500, 1500);
  mygaus->SetParameters(15000, 0, 200);
  double hahe = find_half_height(mygaus);
  cout << hahe << endl;
  cout << mygaus->Eval(hahe)/mygaus->Eval(mygaus->GetParameter(1)) << endl;
  */
