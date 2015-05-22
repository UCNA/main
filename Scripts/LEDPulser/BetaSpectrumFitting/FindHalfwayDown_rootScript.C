{
  TFile f("/data4/saslutsky/OfficialReplayData_06_2014/hists/spec_21805.root");
  TCanvas * rawcan = new TCanvas();
  hRawPMTSpectrum_W0->Draw();
  
  const int NB = 10;
  TVector * linendpoints = new TVector(NB);
  TVector * nbins = new TVector(NB);

  // Identify Maximum height
  double max = hRawPMTSpectrum_W0->GetMaximum();
  int maxbin = hRawPMTSpectrum_W0->GetMaximumBin();

  // Identify location on spectrum where counts fall to half-max
  int halfwaybin;   
  int NB = hRawPMTSpectrum_W0->GetNbinsX();
  // look for max diff of 300, smaller might miss the bin
  hRawPMTSpectrum_W0->GetBinWithContent(max/2., halfwaybin, maxbin, NB, 300); 
  double halfway_energy = hRawPMTSpectrum_W0->GetBinCenter(halfwaybin);
  double halfway_height = hRawPMTSpectrum_W0->GetBinContent(halfwaybin);

  cout << "HALF " << halfwaybin << endl;

  // Do a linear fit around the half-max
  TF1 * linfit = new TF1("linfit", "pol1", -100, 3900);
  //  for (int nb = 2; nb < 10; nb++){ //optimized to be around nb = 6
  int nb = 6;
  double fitmin = hRawPMTSpectrum_W0->GetBinCenter(halfwaybin - nb);
  double fitmax = hRawPMTSpectrum_W0->GetBinCenter(halfwaybin + nb);
  hRawPMTSpectrum_W0->Fit("linfit","R","",fitmin, fitmax);
  
  // Extract "linearized endpoint"
  double * results = linfit->GetParameters();
  double linendpoint = (-1)*results[0]/results[1];
  cout << nb << endl << "linearized endpoint:" << endl << linendpoint << endl;
  //linendpoints(nb) = linendpoint;
  //  nbins(nb) = nb*2;
  }

}
/* TCanvas * graphcan = new TCanvas();
  TGraph graphie = TGraph(*nbins, *linendpoints);
  graphie->Draw("A*");

}
*/
  /*  TF1 * fhalf = new TF1("fhalf", "pol0", -100, 3900);
  fhalf->SetParameter(0, halfway_height);
  fhalf->SetLineStyle(2);
  fhalf->Draw("Same");
  TF1 * ffull= new TF1("ffull", "pol0", -100, 3900);
  ffull->SetParameter(0, max);
  ffull->SetLineStyle(2);
  ffull->SetLineColor(4);
  ffull->Draw("Same");
  
  cout << max << "/" << halfway_height << "=" << max/halfway_height << endl;
  cout << halfway_height << "*2 - " << max << " = " << 2*halfway_height - max << endl;
  */

    
