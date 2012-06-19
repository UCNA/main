#include "EfficCurve.hh"
#include "Types.hh"

Double_t poiscdf(const Double_t *x, const Double_t *par) {
	float w0 = par[2]*par[2];
	float x0 = (x[0]-par[0])/par[1]*par[2] + w0 + 1.0/3.0 - 0.015/w0;
	if(x0<0)
		return 0;
	return (1.0-TMath::Gamma(x0,w0))*par[3];
}

Double_t fancyfish(const Double_t *x, const Double_t *par) {
	float adc50 = par[0];						// ADC at 50% efficiency
	float w = par[1];							// width of efficiency transition
	float n50 = par[2]*adc50/w;					// "number of photoelectrons" at 50% efficiency
	float n0 = n50+sqrt(n50)*(x[0]-adc50)/w;	// "number of photoelectrons" at current point
	if(n0 < 0)
		return 0;
	return TMath::Gamma(n50+1/3-0.015/n50,n0)*par[3];
}

void EfficCurve::genEffic(TH1F* hAll, TH1F* hTrig, bool adcChan) {
	
	float xmin = hAll->GetXaxis()->GetXmin();
	float xmax = hAll->GetXaxis()->GetXmax();
	
	// input spectra plots
	if(defaultCanvas) {
		defaultCanvas->cd();
		defaultCanvas->SetLogy(true);
		hAll->SetLineColor(4);
		if(adcChan)
			hAll->GetXaxis()->SetTitle("ADC channels above pedestal");
		else
			hAll->GetXaxis()->SetTitle("Energy [keV]");
		hAll->GetYaxis()->SetTitle("N Events");
		hAll->Draw();
		hTrig->SetLineColor(2);
		hTrig->Draw("Same");
		printCanvas("Trigger_Effic_Input");
		defaultCanvas->SetLogy(false);
	}
	
	// divide histograms
	gEffic = new TGraphAsymmErrors(hTrig->GetNbinsX());
	gEffic->BayesDivide(hTrig,hAll,"w");
	
	// scan for 50% point
	int b = gEffic->GetN();
	double midx,y;
	while(b > 1) {
		gEffic->GetPoint(--b,midx,y);
		if(y < 0.5)
			break;
	}	
	
	// fit efficiency curve
	TF1 efficfit("efficfit",&fancyfish,xmin,xmax,4);
	efficfit.SetParameter(0,midx);
	efficfit.SetParameter(1,10.0);
	efficfit.SetParameter(2,1.4);
	efficfit.SetParameter(3,0.99);
	efficfit.SetParLimits(0,0,100.0);
	efficfit.SetParLimits(1,2,200.0);
	efficfit.SetParLimits(2,0.1,1000.0);
	efficfit.SetParLimits(3,0.75,1.0);
	
	efficfit.SetLineColor(4);
	printf("Pre-fit threshold guess: %.1f\n",midx);
	gEffic->Fit(&efficfit,"Q");
	
	float_err trigef = float_err(efficfit.GetParameter(3),efficfit.GetParError(3));
	float_err trigc = float_err(efficfit.GetParameter(0),efficfit.GetParError(0));
	float_err trigw = float_err(efficfit.GetParameter(1),efficfit.GetParError(1));
	float_err trign_adj = float_err(efficfit.GetParameter(2),efficfit.GetParError(2));
	float trign = trigc.x/trigw.x*trign_adj.x;
	for(unsigned int i=0; i<4; i++) params[i]=efficfit.GetParameter(i);
	
	printf("Poisson CDF Fit: h = %.4f(%.4f), x0 = %.1f(%.1f), dx = %.1f(%.1f), n = %.2f [adjust %.2f(%.2f)]\n",
		   trigef.x, trigef.err, trigc.x, trigc.err, trigw.x, trigw.err, trign, trign_adj.x, trign_adj.err);
	
	// draw fit plot
	if(defaultCanvas) {
		gEffic->SetMinimum(-0.10);
		gEffic->SetMaximum(1.10);
		gEffic->Draw("AP");
		gEffic->SetTitle("Trigger Efficiency");
		if(adcChan)
			gEffic->GetXaxis()->SetTitle("ADC channels above pedestal");
		else
			gEffic->GetXaxis()->SetTitle("Energy [keV]");
		gEffic->GetXaxis()->SetLimits(xmin,xmax);
		gEffic->GetYaxis()->SetTitle("Efficiency");
		gEffic->Draw("AP");
		printCanvas("Trigger_Efficiency");
	}
}

double EfficCurve::effic(double x) const {
	//if(gEffic) return gEffic->Eval(x);
	return fancyfish(&x, params);
}

void EfficCurve::invertEffic(TH1F* hIn, float th) {
	TH1F* hInC = NULL;
	if(defaultCanvas)
		hInC = new TH1F(*hIn);
	for(int i=0; i<hIn->GetNbinsX(); i++) {
		float e = effic(hIn->GetBinCenter(i));
		if(e<th)
			e = th;
		hIn->SetBinContent(i,hIn->GetBinContent(i)/e);
		hIn->SetBinError(i,hIn->GetBinError(i)/e);
	}
	if(defaultCanvas) {
		defaultCanvas->cd();
		hInC->SetLineColor(4);
		hInC->Draw();
		hIn->SetLineColor(2);
		hIn->Draw("Same");
		printCanvas("Corrected_Spectrum");
		delete(hInC);
	}
	
}

