#include "StyleSetup.hh"
#include "GraphicsUtils.hh"

void ROOTStyleSetup() {
	gROOT->SetStyle("Plain");
	gStyle->SetPalette(1);
	gStyle->SetNumberContours(255);
	gStyle->SetOptStat("e");
	TCanvas* defaultCanvas = new TCanvas();
#ifdef PUBLICATION_PLOTS
	gStyle->SetOptStat("");
	makeGrayscalepalette();
	defaultCanvas->SetGrayscale(true);
#endif
	defaultCanvas->SetFillColor(0);
	defaultCanvas->SetCanvasSize(300,300);
}
	