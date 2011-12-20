#include "SRAsym.hh"

#include <cassert>
#include <TPad.h>
#include "PMTGenerator.hh"
#include "BetaSpectrum.hh"

SRAsym::SRAsym(TH1F* eOff, TH1F* wOff, TH1F* eOn, TH1F* wOn, bool bonehead): hAsym(NULL), hSum(NULL), boneit(bonehead), afit(NULL), nafit(NULL) {
	
	// gather data
	hOn[EAST] = eOn;
	hOn[WEST] = wOn;
	hOff[EAST] = eOff;
	hOff[WEST] = wOff;
	// make sure errors are calculated
	hOn[EAST]->Sumw2();
	hOn[WEST]->Sumw2();
	hOff[EAST]->Sumw2();
	hOff[WEST]->Sumw2();
	
	// scale to unit area
	float sfactor = 1.0/hOn[EAST]->Integral("width");
	hOn[EAST]->Scale(sfactor);
	hOn[WEST]->Scale(sfactor);
	hOff[EAST]->Scale(sfactor);
	hOff[WEST]->Scale(sfactor);
	
	if(boneit)
		calcBonehead();
	else
		calcAsym();
	normAsym();
}

void SRAsym::fitAsym(float emin, float emax) {
	assert(hAsym);
	assert(hAnorm);
	if(afit) delete(afit);
	afit = new TF1("AsymFit","pol0",emin,emax);
	hAsym->Fit(afit,"QR");
	if(boneit)
		printf("Bonehead Asymmetry A = %.4f~%.4f\n",afit->GetParameter(0),afit->GetParError(0));
	else
		printf("Super Ratio Asymmetry A = %.4f~%.4f\n",afit->GetParameter(0),afit->GetParError(0));
	if(nafit) delete(nafit);
	nafit = new TF1("nAsymFit","pol0",emin,emax);
	hAnorm->Fit(nafit,"QR");
}

Stringmap SRAsym::toStringmap() const {
	Stringmap m;
	if(afit) {
		m.insert("AsymFit",afit->GetParameter(0));
		m.insert("dAsymFit",afit->GetParError(0));
		m.insert("fitMin",afit->GetXmin());
		m.insert("fitMax",afit->GetXmax());
		m.insert("AsymNormFit",nafit->GetParameter(0));
		m.insert("dAsymNormFit",nafit->GetParError(0));
	}
	return m;
}

void SRAsym::drawSpectra(int lcolor,bool same,std::string outname) {
	printf("Drawing asymmetry spectra for '%s'\n",outname.c_str());
	if(lcolor > 0) {
		hOff[WEST]->SetLineColor(lcolor);
		hOff[EAST]->SetLineColor(lcolor);
		hOn[EAST]->SetLineColor(lcolor);
		hOn[WEST]->SetLineColor(lcolor);
	} else {
		hOff[WEST]->SetLineColor(3);
		hOff[EAST]->SetLineColor(2);
		hOn[EAST]->SetLineColor(4);
		hOn[WEST]->SetLineColor(5);
	}
	if(same)
		hOff[WEST]->Draw("Same");
	else
		hOff[WEST]->Draw();
	hOff[EAST]->Draw("Same");
	hOn[EAST]->Draw("Same");
	hOn[WEST]->Draw("Same");
	if(outname.size())
		gPad->Print(outname.c_str());
}

void SRAsym::calcAsym() {
	// calculate super-ratio S
	if(hAsym) delete(hAsym);
	hAsym = (TH1F*)hOff[EAST]->Clone("hAsym");
	hAsym->Multiply(hOn[WEST]);
	hAsym->Divide(hOff[WEST]);
	hAsym->Divide(hOn[EAST]);
	
	// convert S to asymmetry
	for(int n=1; n<=hAsym->GetNbinsX(); n++) {
		float s = hAsym->GetBinContent(n);
		float ds = hAsym->GetBinError(n);
		hAsym->SetBinContent(n,(1.0-sqrt(s))/(1.0+sqrt(s)));
		hAsym->SetBinError(n,ds/(sqrt(s)*(1.0+sqrt(s))*(1.0+sqrt(s))));
	}
}

void SRAsym::normAsym() {
	hAnorm = (TH1F*)hAsym->Clone("hAnorm");
	for(int i=1; i<=hAnorm->GetNbinsX(); i++) {
		float b = beta(hAnorm->GetBinCenter(i));
		hAnorm->SetBinContent(i,hAnorm->GetBinContent(i)/b);
		hAnorm->SetBinError(i,hAnorm->GetBinError(i)/b);
	}
}

void SRAsym::calcSum() {
	if(hSum)
		delete(hSum);
	hSum = (TH1F*)hOff[EAST]->Clone("hSum");
	hSum->Add(hOff[WEST]);
	hSum->Add(hOn[EAST]);
	hSum->Add(hOn[WEST]);
}

void SRAsym::calcBonehead() {
	
	calcSum();
	
	if(hAsym)
		delete(hAsym);
	hAsym = (TH1F*)hOn[EAST]->Clone("hAsym");
	hAsym->Add(hOff[WEST]);
	hAsym->Add(hOn[WEST],-1.0);
	hAsym->Add(hOff[EAST],-1.0);
	
	hAsym->Divide(hSum);
}

