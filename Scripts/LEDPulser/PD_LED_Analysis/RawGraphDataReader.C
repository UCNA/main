{
  gStyle->SetOptFit(1111);

  ifstream infile;
  infile.open("./images_dmp/RawGraphData_21927.txt");
  float run, tube, LED, PD, PDE;// PMT, PMTE; // nan's screwing things up
  float PMTf;
  string PMT, PMTE;
  vector <float> PD_v, x_v, PDE_v, xE_v, PD_PMT_v, PMT_v;
  int counter = 0;
  tube = -1;
  LED = -1;
  infile >> run >> tube >> LED >> PD >> PDE >> PMT >> PMTE;
  while (1){
    if (LED == 0 && tube == 0){
      cout << run << tube << LED << PD << PDE << PMT << PMTE << endl;
      PD_v.push_back(PD); PDE_v.push_back(PDE);
      x_v.push_back(counter); xE_v.push_back(0);
      if (PMT != "nan"){ 
	const char* PMTc = PMT.c_str();
	PMTf = atof(PMTc);
	PMT_v.push_back(PMTf);
	PD_PMT_v.push_back(PMTf/PD/PD);
	counter++;
	cout << "Count: " << counter << endl;
      }
    }
    if (!infile.good()) break;
    infile >> run >> tube >> LED >> PD >> PDE >> PMT >> PMTE;
  }

  int graphsize = PD_v.size();
  int gratiosize = PD_PMT_v.size();
  if (graphsize){
    //TGraph *gr = new TGraph(graphsize, &x_v[0], &PD_v[0]);
    TGraphErrors *gr = new TGraphErrors(graphsize, &x_v[0], &PD_v[0], &xE_v[0], &PDE_v[0]);
    TGraph *gratio = new TGraph(gratiosize, &x_v[0], &PD_PMT_v[0]);
    TGraph *gPMT = new TGraph(gratiosize, &x_v[0], &PMT_v[0]);
    TCanvas * c1 = new TCanvas("c1", "c1");
    c1->Divide(1,3);
    c1->Update();
    c1->cd(1);
    gr->Draw("A*");
    c1->Update();
    c1->cd(2);
    // gr->Fit("pol2", "R","", 17, 25);
    gratio->Draw("A*");
    c1->cd(3);
    gPMT->Draw("A*");
    c1->Update();

  }
  else cout << "DOH" << endl;
  //  gr->Fit("pol2");
} 
