#include "KurieStudy.hh"
#include "KurieFitter.hh"
#include <TH2F.h>

using namespace std;

float_err testKurie(){
  float beta_endpoint = 782; // keV
  
  cout << "Kurie Testing" << endl;
  TFile * infile = 
    new TFile ("/data4/saslutsky/OfficialReplayData_06_2014/hists/spec_23171.root");
  TH1F * hRaw = (TH1F*)infile->Get("hRawPMTSpectrum_E2")->Clone("hRaw");
  float_err ep = kurieIterator(hRaw, 800, NULL, beta_endpoint, 250 , 700, 500);
  //  float_err ep = kurieIterator(hRaw, 800, NULL, 1100, 400 , 700, 500);
  
  return ep;
}

float_err makeKurieFit(TH1F* hIn){
  float beta_endpoint = 782; // keV
  //  cout << "makeKurieFit" << endl;
  //  TCanvas * c1 = new TCanvas("c1","c1");
  
  //  TGraphErrors ** tg;
 
  float_err ep = kurieIterator(hIn, 800, NULL, beta_endpoint, 250 , 700, 500);
  //  float_err ep = kurieIterator(hIn, 800, tg, beta_endpoint, 250 , 700, 500);
  //*tg.Draw("A*");
  // string hName = hIn->GetName();
  //hName=hName + "_graph.png";
  //c1->SaveAs(hName.c_str());

  // TGraphErrors *tgout = *tg;
  //tgout->Draw("A*");
  // c1->SaveAs("hi.png");

  return ep;
}

float_err makeKurieFitsforRun(std::string runnum){
  const char * base_char = getenv("UCNAOUTPUTDIR");
  string base_path = base_char;
  string full_path = base_path + "hists/spec_" + runnum + ".root";
  TFile * infile =  new TFile (full_path.c_str());
  
  cout << runnum << "\t";

  float_err ep;
  for (int s = 0; s < 2; s++){
    char * side = "E";
    if (s > 0){
      side = "W";
    }
    for (int i = 0; i < 4; i++){
      char *tubeID = Form("%s%i", side, i);
      //      cout << tubeID << endl;
      TH1F * hRaw = (TH1F*)infile->Get(Form("hRawPMTSpectrum_%s", tubeID))->Clone("hRaw");
      ep = makeKurieFit(hRaw);
      cout << ep.x;
      if (i + 4*s < 7){
	cout << "\t";
      }
    }
  }
  cout << endl;
  return ep;
}

