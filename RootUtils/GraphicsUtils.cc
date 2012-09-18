#include "GraphicsUtils.hh"
#include "strutils.hh"
#include <cassert>
#include <algorithm>

#include <TEllipse.h>
#include <TPad.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLine.h>
#include <TLatex.h>
#include <TBox.h>
#include <TStyle.h>
#include <TColor.h>

bool compareHistosByMax(TH1* i, TH1* j) {
	assert(i && j);
	return i->GetMaximum() < j->GetMaximum();
}

double getXmin(const TH1* h) { assert(h); return h->GetBinLowEdge(1); }
double getXmax(const TH1* h) { assert(h); return h->GetBinLowEdge(h->GetNbinsX()+1); }

bool compareHistosByXmin(TH1* i, TH1* j) {
	assert(i && j);
	return getXmin(i) < getXmin(j);
}
bool compareHistosByXmax(TH1* i, TH1* j) {
	assert(i && j);
	return getXmax(i) < getXmax(j);
}

double drawSimulHistos(std::vector<TH1*>& hists, const std::string& opt, const std::string& newTitle) {
	if(!hists.size())
		return 0;
	printf("Drawing %i histograms together",(int)hists.size()); fflush(stdout);
	double xmin = getXmin(*std::min_element(hists.begin(),hists.end(),compareHistosByXmin));
	printf(" spanning x=(%g,",xmin); fflush(stdout);
	double xmax = getXmax(*std::max_element(hists.begin(),hists.end(),compareHistosByXmax));
	printf("%g), ",xmax); fflush(stdout);
	TH1* maxHist = *std::max_element(hists.begin(),hists.end(),compareHistosByMax);
	assert(maxHist);
	printf("with ymax = %g...",maxHist->GetMaximum()); fflush(stdout);
	//maxHist->SetAxisRange(xmin,xmax,"X");
	std::string oldTitle = maxHist->GetTitle();
	if(newTitle != "DEFAULT")
		maxHist->SetTitle(newTitle.c_str());
	maxHist->Draw(opt.c_str());
	for(std::vector<TH1*>::iterator it = hists.begin(); it != hists.end(); it++) {
		assert(*it);
		if(*it == maxHist)
			continue;
		if(opt.size())
			(*it)->Draw((opt+" SAME").c_str());
		else
			(*it)->Draw("SAME");
	}
	printf(" Done.\n");
	
	maxHist->SetTitle(oldTitle.c_str());
	return maxHist->GetMaximum();
}

void drawHistoPair(TH1* hRed, TH1* hBlue, const std::string& opt, Int_t c1, Int_t c2) {
	assert(hRed && hBlue);
	hRed->SetLineColor(c1);
	hRed->SetMarkerColor(c1);
	hBlue->SetLineColor(c2);
	hBlue->SetMarkerColor(c2);
	std::vector<TH1*> hToPlot;
	hToPlot.push_back(hRed);
	hToPlot.push_back(hBlue);
	drawSimulHistos(hToPlot,opt);
}

void drawFiducialCuts(Int_t color) {
	// decay tube (5")
	TEllipse* e1 = new TEllipse(0,0,decay_tube_radius,decay_tube_radius);
	e1->SetFillStyle(0);
	e1->SetLineColor(color);
	e1->SetLineWidth(1);
	e1->SetLineStyle(2);
	e1->Draw();
	
	// 50mm nominal acceptance
	TEllipse* e2 = new TEllipse(0,0,50,50);
	e2->SetFillStyle(0);
	e2->SetLineColor(color);
	e2->SetLineStyle(2);
	e2->Draw();
}

void drawEllipseCut(Source E, Float_t nSigma, std::string label) {
	
	TEllipse* e = new TEllipse(E.x,E.y,nSigma*E.wx,nSigma*E.wy);
	e->SetFillStyle(0);
	e->SetLineColor(2);
	e->SetLineWidth(2);
	e->Draw();
	
	if(label.size()) {
		TLatex lx;
		lx.SetTextColor(2);
		lx.SetTextSize(0.02);
		lx.DrawLatex(E.x+nSigma*E.wx,E.y+nSigma*E.wy,label.c_str());
	}
	
}


void drawVLine(Float_t x, TVirtualPad* C, Int_t color) {
	Double_t xmin,ymin,xmax,ymax;
	C->Update();
	C->GetRangeAxis(xmin,ymin,xmax,ymax);
	if(C->GetLogy()) {
		ymin = pow(10,ymin);
		ymax = pow(10,ymax);
	}
	TLine* l = new TLine(x,ymin,x,ymax);
	l->SetLineColor(color);
	l->Draw();
}

void drawExcludedRegion(Float_t x0, Float_t x1, TCanvas* C, Int_t color, Int_t fill) {
	Double_t xmin,ymin,xmax,ymax;
	C->Update();
	C->GetRangeAxis(xmin,ymin,xmax,ymax);
	if(C->GetLogy()) {
		ymin = pow(10,ymin);
		ymax = pow(10,ymax);
	}
	TBox* r = new TBox(x0,ymin,x1,ymax);
	r->SetFillColor(color);
	r->SetFillStyle(fill);
	r->Draw();
}

void drawSectors(const SectorCutter& S, int color) {
	
	TEllipse* e = new TEllipse(0,0,S.r,S.r);
	e->SetFillStyle(0);
	e->SetLineColor(color);
	e->Draw();
	
	for(unsigned int i=1; i<S.n; i++) {
		float r1 = S.ringRadius(i-1);
		float r2 = S.ringRadius(i);
		TEllipse* e = new TEllipse(0,0,r1,r1);
		e->SetFillStyle(0);
		e->SetLineColor(color);
		e->Draw();
		for(unsigned int j=0; j<S.ndivs[i]; j++) {
			float ph = (6.28318531*j)/S.ndivs[i];
			TLine* l = new TLine(r1*cos(ph),r1*sin(ph),r2*cos(ph),r2*sin(ph));
			l->SetLineColor(color);
			l->Draw();
		}
	}
}

void labelSectors(const SectorCutter& S, int color) {
	TLatex l;
	float x,y;
	for(unsigned int i=0; i<S.nSectors(); i++) {
		S.sectorCenter(i,x,y);
		l.SetTextColor(color);
		l.SetTextAlign(22);
		l.SetTextSize(0.02);
		l.DrawLatex(x,y,itos(i).c_str());
	}
}

void makeGrayscalepalette() {
	const unsigned int ncol = 256;
	Int_t cnum[ncol];
	for(unsigned int i=0; i<ncol; i++)
		cnum[i] = TColor::GetColor(float(i)/float(ncol-1),float(i)/float(ncol-1),float(i)/float(ncol-1));
	gStyle->SetPalette(ncol,cnum);
	gStyle->SetNumberContours(64);	
}

