#include "StyleSetup.hh"
#include "GraphicsUtils.hh"

void ROOTStyleSetup(bool b2w) {
	gROOT->SetStyle("Plain");
	gStyle->SetPalette(1);
	gStyle->SetNumberContours(255);
	gStyle->SetOptStat("e");
	TCanvas* defaultCanvas = new TCanvas();
#ifdef PUBLICATION_PLOTS
	gStyle->SetOptStat("");
	makeGrayscalepalette(b2w);
	defaultCanvas->SetGrayscale(true);
#endif
	defaultCanvas->SetFillColor(0);
	defaultCanvas->SetCanvasSize(300,300);
}
	