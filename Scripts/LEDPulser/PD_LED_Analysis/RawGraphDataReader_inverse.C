{
  gStyle->SetOptFit(1111);

  ifstream infile;
  infile.open("./images_dmp/RawGraphData_22767.txt");
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
//	PD_PMT_v.push_back(PMTf*PMTf/PD/PD);
//	PD_PMT_v.push_back(PMTf/PD);
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

//// PD versus x
    TGraphErrors *gr = new TGraphErrors(graphsize, &PD_v[0], &x_v[0], &PDE_v[0], &xE_v[0]);
    TGraph *gratio = new TGraph(gratiosize, &PD_PMT_v[0], &x_v[0]);
    TGraph *gPMT = new TGraph(gratiosize, &PMT_v[0], &x_v[0]);
    TCanvas * c1 = new TCanvas("c1", "c1", 1800, 1000);
    c1->Divide(1,4);
    c1->Update();
    c1->cd(1);
    gr->Draw("A*");

    TF1 * fit = new TF1("fit", "[0]/x/x + [1]/x + [2] + [3]*x + [4]*x*x", 0, 300);
//    TF1 * fit = new TF1("fit", "[0]/x + [1] + [2]*x + [3]*x*x", 0, 300);
//    TF1 * fit = new TF1("fit", "[0] + [1]*x + [2]*x*x + exp([3]+[4]*x)", 0, 300);

//      TF1 * fit = new TF1("fit", "[0] + [1]*x + [2]*x*x + [3]*log([4]+[5]*x)", 0, 300);
//    TF1 * fit = new TF1("fit", "[0] + [1]*x + [2]*exp([3]*x)", 0, 50);
//    TF1 * fit = new TF1("fit", "[0] + [1]*x + [2]*x*x", 0, 50);
 
    fit->SetParameters(0.0, -0.1, 0.0, -3.0, 0.0, 0.0); 
    gr->Fit("fit", "R", "", 10, 200); // LED = 0    
    //gr->Fit("fit", "R", "", 5.5, 35); // LED = 1

 /*   Double_t * fitparms = fit->GetParameters();
    TF1 * fit_pol = new TF1("fit_pol", "[0]*x + [1]*x*x", 0, 300);
    TF1 * fit_exp = new TF1("fit_exp", "[0] + exp([1]+[2]*x)", 0, 300);
    fit_pol->SetParameter(0, fitparms[1]);
    fit_pol->SetParameter(1, fitparms[2]);
    fit_exp->SetParameter(0, fitparms[0]);
    fit_exp->SetParameter(1, fitparms[3]);
    fit_exp->SetParameter(2, fitparms[4]);
//    fit_exp->SetParameter(3, fitparms[5]);

    //for (int p = 0; p < 3; p++) fit_pol->SetParameter(p, fitparms[p]);
    //for (int p = 0; p < 2; p++) fit_exp->SetParameter(p, fitparms[p+3]);
    fit_pol->SetLineColor(4);
    fit_exp->SetLineColor(4);
    fit_exp->Draw("same");
    fit_pol->Draw("same");
*/

    c1->Update();
    c1->cd(2);
    gratio->Draw("A*");

//// PMT/PD^2 versus x
//    TF1 * fitrat = new TF1("fitrat", "[0] + [1]*x + [2]*exp([3]*x)", 0, 50);
    TF1 * fitrat = new TF1("fitrat", "[0] + [1]*x + [2]*x*x + exp([3]+[4]*x)", 0, 50);
//    TF1 * fitrat = new TF1("fitrat", "[0] + [1]*exp([2]*x)", 0, 50);
//    TF1 * fitrat = new TF1("fitrat", "[0] + [1]*x + [2]*x*x", 0, 50);
    fitrat->SetParameters(0, 0.1, 0.01, 1, 0.5);
    gratio->Fit("fitrat", "R", "", 0.5, 24.5); // LED = 0
    //gratio->Fit("fitrat", "R", "", 5.5, 35); // LED = 1
    gratio->GetYaxis()->SetRangeUser(0, 2);

    Double_t * fitratparms = fitrat->GetParameters();
    TF1 * fitrat_pol = new TF1("fitrat_pol", "pol2", 0, 24);
    TF1 * fitrat_exp = new TF1("fitrat_exp", "expo", 0, 24);
    for (int p = 0; p < 3; p++) fitrat_pol->SetParameter(p, fitratparms[p]);
    for (int p = 0; p < 2; p++) fitrat_exp->SetParameter(p, fitratparms[p+3]);
    fitrat_pol->SetLineColor(4);
    fitrat_exp->SetLineColor(4);
    fitrat_exp->Draw("same");
    fitrat_pol->Draw("same");


//// PMT versus x
    c1->cd(3);    
//    TF1 * fitpmt = new TF1("fitpmt", "[0] + [1]*x + [2]*x*x + exp([3] + [4]*x)", 0, 5000);
    TF1 * fitpmt = new TF1("fitpmt", "[0] + [1]*x + [2]*x*x", 0, 50);
    gPMT->Draw("A*");
    gPMT->Fit("fitpmt", "R", "", 5, 2500); // LED = 0    
    //gratio->Fit("fitrat", "R", "", 5.5, 35); // LED = 1
    
/*    Double_t * fitpmtparms = fitpmt->GetParameters();
    TF1 * fitpmt_pol = new TF1("fitpmt_pol", "[0]*x + [1]*x*x", 0, 5000);
    TF1 * fitpmt_exp = new TF1("fitpmt_exp", "[0] + exp([1]+[2]*x)", 0, 5000);
    fitpmt_pol->SetParameter(0, fitpmtparms[1]);
    fitpmt_pol->SetParameter(1, fitpmtparms[2]);
    fitpmt_exp->SetParameter(0, fitpmtparms[0]);
    fitpmt_exp->SetParameter(1, fitpmtparms[3]);
    fitpmt_exp->SetParameter(2, fitpmtparms[4]);
*/
/*    TF1 * fitpmt_pol = new TF1("fitpmt_pol", "pol2", 0, 50);
    TF1 * fitpmt_exp = new TF1("fitpmt_exp", "expo", 0, 50);
    for (int p = 0; p < 3; p++) fitpmt_pol->SetParameter(p, fitpmtparms[p]);
    for (int p = 0; p < 2; p++) fitpmt_exp->SetParameter(p, fitpmtparms[p+3]);
  */  
/*    fitpmt_pol->SetLineColor(4);
    fitpmt_exp->SetLineColor(4);
    fitpmt_pol->Draw("same");
    fitpmt_exp->Draw("same");
*/   
    c1->Update();

//// PD versus PMT
    c1->cd(4);
    TGraphErrors *grpdpmt = new TGraphErrors(graphsize, &PMT_v[0], &PD_v[0]);
    grpdpmt->Draw("A*");
    c1->Update();
  }
  else cout << "DOH" << endl;

} 
