
{
#include <vector>
#include <cstdlib>

  string year = "2011";
  TCanvas *c0 = new TCanvas("c0", "c0", 1350,700);
  gStyle->SetOptStat(11);
  gStyle->SetOptFit(1);
  
  const char *path = getenv("UCNA_ANA_PLOTS");
  cout << path << endl;

  std::vector <string> srcs;
  std::vector <string> titles;

  if (year=="2011")
    {
      string srcs_hold[] = {"Bi 1", "Bi 2", "Sn", "Cd", "Ce", "In"};
      string titles_hold[] = {"^{207}Bi 1", "^{207}Bi 2", "^{113}Sn", "^{109}Cd", "^{139}Ce", "^{114}In"};
      srcs.resize(6);
      titles.resize(6);
      c0->Divide(3,2);
      for (unsigned int k=0; k<6; k++)
	{
	  srcs[k]= srcs_hold[k];
	  titles[k] = titles_hold[k];
	}
    }
  else if (year=="2012")
    {
      string srcs_hold[] = {"Bi 1", "Bi 2", "Sn", "Cd", "Ce", "In", "Cs"};
      string titles_hold[] = {"^{207}Bi 1", "^{207}Bi 2", "^{113}Sn", "^{109}Cd", "^{139}Ce", "^{114}In", "^{137}Cs"};
      srcs.resize(7);
      titles.resize(7);
      c0->Divide(4,2);
      for (unsigned int k=0; k<7; k++)
	{
	  srcs[k]= srcs_hold[k];
	  titles[k] = titles_hold[k];
	}
    }
      
  std::vector <TH1F*> hists(srcs.size());
  cout << "Size of hists is " << hists.size() << endl;

  for (unsigned int i=0; i<hists.size();i++) {
    c0->cd(i+1);
    hists[i] = new TH1F(srcs[i].c_str(),titles[i].c_str(),40,-25.,25.);
    hists[i]->GetXaxis()->SetTitle("Energy Error (keV)");
    string fname = string(path) + "/" + srcs[i] + "_" + year +".txt";
    ifstream infile;
    infile.open(fname.c_str());
    cout << "Opened file " << fname.c_str() << endl;
    double val=0.;
    while (infile >> val) {hists[i]->Fill(val);}
    infile.close();
    hists[i]->Fit("gaus");
    hists[i]->Draw();
    //fname = "";  
  }
  
}

/*TH1F *Bi1 = new TH1F("Bi1", "^{207}Bi 1", 100, 0., high0);
  east0->GetXaxis()->SetTitle("Energy Error (keV)");
  TH1F *Bi2 = new TH1F("Bi2", "^{207}Bi 2", 100, 0., high0);
  east0->GetXaxis()->SetTitle("Energy Error (keV)");
  TH1F *Sn = new TH1F("Sn", "^{113}Sn", 100, 0., high0);
  east0->GetXaxis()->SetTitle("Energy Error (keV)");
  TH1F *Cd = new TH1F("Cd", "^{109}Cd", 100, 0., high0);
  east0->GetXaxis()->SetTitle("Energy Error (keV)");
  TH1F *Ce = new TH1F("Ce", "^{139}Ce", 100, 0., high0);
  east0->GetXaxis()->SetTitle("Energy Error (keV)");
  TH1F *In = new TH1F("In", "^{114}In", 100, 0., high0);
  east0->GetXaxis()->SetTitle("Energy Error (keV)");*/
