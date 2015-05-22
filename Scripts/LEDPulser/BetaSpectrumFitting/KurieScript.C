{
  TFile f("/data4/saslutsky/OfficialReplayData_06_2014/hists/spec_21805.root");
  TCanvas * rawcan = new TCanvas();
  hRawPMTSpectrum_E1->Draw();
  int NB = hRawPMTSpectrum_E1->GetNbinsX();

  TF1 * fKurie = new TF1("fKurie", "[1]*(([0]*x + 511)*sqrt([0]*[0]*x*x + 2*x*[0]*511)*(782-[0]*x)*(782-[0]*x))", 0, 4000);
  
  const int NIter = 50;
  TVector * scale = new TVector(NIter);
  TVector * endpoint = new TVector(NIter);
  TVector * halfwaypoint  = new TVector(NIter);
  TH1F * hDivided[NIter]; 

  TCanvas * histcan = new TCanvas();
  histcan->Divide(7,7);
  
  for (int N = 1; N < NIter; N++){
    cout << "Plot " << N << endl;
    //  for (int N = 1; N < 2; N++){
    hDivided[N] = (TH1F*)hRawPMTSpectrum_E1->Clone(Form("hDivided_%i", N));
    hDivided[N]->SetTitle(Form("%i", N));
    
    double scale_N = N*2./NIter;
    fKurie->SetParameters(scale_N, 1.2e-12);
    scale(N) = scale_N;
    
    double rawcontent = 0;
    double bincenter = 0;
    double divcontent = 0;
    double firstmaxcontent = 0;
    double firstmaxbin = 0;
    double halfwaycontent = 0;
    double halfwaybin = 0;
    bool passedhalfway = 0;
    bool passedmax = 0;

    for (int i = 0; i < NB; i++){
      rawcontent = hRawPMTSpectrum_E1->GetBinContent(i);
      bincenter  = hRawPMTSpectrum_E1->GetBinCenter(i);
      if (fKurie(bincenter) == 0){
	hDivided[N]->SetBinContent(i, 0);
	continue;
      }
      
      divcontent = sqrt(rawcontent/fKurie(bincenter));
      hDivided[N]->SetBinContent(i, divcontent);
       
      /*      if (bincenter > 30){ //(stay away from pedestal rise)
	if (divcontent < mincontent){
	  if (passedMin == 0){
	    mincontent = divcontent;
	    minbin = bincenter;
	  }
	}
	if (divcontent > mincontent){
	  passedMin = 1;
	}
      } */
      
      if (divcontent > firstmaxcontent){ // identify max 
	if (passedmax == 0){
	  firstmaxcontent = divcontent;
	  firstmaxbin = bincenter;
	}
      }
      if (divcontent < firstmaxcontent){ // stop looking for max
	passedmax = 1;
      }
      if (passedmax == 1 && passedhalfway == 0){ // identify "half-way down"
	double halfway = 0.5*firstmaxcontent;
	if (divcontent < halfway){
	  if (i > 0){
	    passedhalfway = 1;
	    halfwaycontent = divcontent;
	    
	    //halfwaybin = bincenter; // not robust to bin changes
	    // interpolate betwen the two nearest bins instead: UPDATE: This isn't robust against bin changes either
	    xi_prev = hDivided[N]->GetBinCenter(i-1);
	    yi_prev = hDivided[N]->GetBinContent(i-1);
	    halfwaybin = xi_prev + (halfway - yi_prev)*( (bincenter-xi_prev)/(divcontent - yi_prev) );
	    cout << "(" << xi_prev << "," << yi_prev << ")" << endl;
	    cout << halfway << endl;
	    cout << halfwaybin << endl;
	  }
	}
      }
      halfwaypoint(N) = halfwaybin;
    }
    rawcan->cd();
    //    hDivided[N]->Draw("same");
    hDivided[N]->SetMaximum(500);

    double fit_min = firstmaxbin;   
    double fit_max = 1.6*(halfwaybin - firstmaxbin);  // stay away from upturn
    TF1 * fitpol = new TF1("fitpol", "pol1", fit_min, fit_max);
    histcan->cd(N);
    hDivided[N]->Fit("fitpol","R", "", fit_min, fit_max);
    double endpoint_N = (-1)*fitpol->GetParameter(0)/fitpol->GetParameter(1);
    endpoint(N) = endpoint_N;
    
    TF1 * fhalf = new TF1("fhalf", "pol0", 0, 4000);
    fhalf->SetParameter(0, halfwaycontent);
    fhalf->SetLineColor(4); fhalf->SetLineStyle(2);
    fhalf->Draw("Same");
  }
    
  TCanvas * graphcan = new TCanvas();
  graphcan->Divide(1,2);
  graphcan->cd(1);
  TGraph graphie = TGraph(*scale, *endpoint);
  graphie.Draw("A*");
  TGraph graphhalf = TGraph(*scale, *halfwaypoint);
  graphcan->cd(2);
  graphhalf.Draw("A*");
}
