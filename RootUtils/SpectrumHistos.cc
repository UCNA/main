#include "SpectrumHistos.hh"

unsigned int SpectrumHistos::nHistos = 0;
unsigned int SpectrumHistos::nHistos2 = 0;

void SpectrumHistos::setupHistos(const std::string& hName, unsigned int nbins, float xmin, float xmax) {
	printf("Processing FG/BG pair '%s'\n",hName.c_str());
	chName = hName;
	
	// define spectrum histograms
	fg = registeredTH1F(hName+"_FG_"+itos(++nHistos),hName+" Foreground",nbins,xmin,xmax);
	fg->Sumw2();
	bg = registeredTH1F(hName+"_BG_"+itos(++nHistos),hName+" Background",nbins,xmin,xmax);
	bg->Sumw2();
	bg->SetLineColor(2);
}

void SpectrumHistos::finishHistos(float fgNorm, float bgNorm) {
	
	// normalization
	fgCount = fg->Integral();
	fg->Scale(1.0/(fgNorm*fg->GetBinWidth(1)));
	bgCount = bg->Integral();
	if(bgNorm)
		bg->Scale(1.0/(bgNorm*bg->GetBinWidth(1)));
	
	// bg subtraction
	diff = new TH1F(*fg);
	diff->SetLineColor(4);
	diff->Add(bg,-1.0);
	addObject(diff);
	
	Stringmap m;
	m.insert("fgCounts",fgCount);
	m.insert("bgCounts",bgCount);
	m.insert("fgNorm",fgNorm);
	m.insert("bgNorm",bgNorm);
	m.insert("fgRate",fg->Integral("width"));
	m.insert("bgRate",bg->Integral("width"));
	m.insert("diffRate",diff->Integral("width"));
	qOut.insert(std::string("Spectrum_")+chName,m);
	
	printf("** FG: %f Hz, BG: %f Hz, Total: %f Hz [%i/%i/%i events] **\n",
		   fg->Integral("width"),bg->Integral("width"),diff->Integral("width"),
		   (int)fgCount,(int)bgCount,(int)(fgCount-bgCount));
	
}

void SpectrumHistos::fill(float x, bool isFg) {
	if(isFg)
		fg->Fill(x);
	else
		bg->Fill(x);
}

void SpectrumHistos::fill(float x, float y, bool isFg) {
	if(isFg)
		fg2->Fill(x,y);
	else
		bg2->Fill(x,y);
}

void SpectrumHistos::fill(float fgNorm, float bgNorm,
						  const std::vector<float>& fgPoints, const std::vector<float>& bgPoints,
						  const std::string& hName, unsigned int nbins, float xmin, float xmax) {
	
	setupHistos(hName,nbins,xmin,xmax);
	
	// fill FG spectrum histograms
	for(UInt_t e=0; e<fgPoints.size(); e++)
		fill(fgPoints[e],true);
	
	// fill BG spectrum histograms
	for(UInt_t e=0; e<bgPoints.size(); e++)
		fill(bgPoints[e],false);
	
	finishHistos(fgNorm, bgNorm);
}

void SpectrumHistos::setup2D(const std::string& hName,
							 unsigned int nbinsX, float xmin, float xmax,
							 unsigned int nbinsY, float ymin, float ymax) {
	
	printf("Processing FG/BG pair '%s'\n",hName.c_str());
	chName = hName;
	
	// define spectrum histograms
	fg2 = registeredTH2F(hName+"_FG_"+itos(++nHistos2),hName+" Foreground",nbinsX,xmin,xmax,nbinsY,ymin,ymax);
	fg2->Sumw2();
	bg2 = registeredTH2F(hName+"_BG_"+itos(++nHistos2),hName+" Background",nbinsX,xmin,xmax,nbinsY,ymin,ymax);
	bg2->Sumw2();	
}

void SpectrumHistos::finish2D(float fgNorm, float bgNorm) {
	fgCount = fg2->Integral();
	fg2->Scale(1.0/fgNorm);
	bgCount = bg2->Integral();
	if(bgNorm)
		bg2->Scale(1.0/bgNorm);
	diff2 = new TH2F(*fg2);
	diff2->Add(bg2,-1.0);
	addObject(diff2);
}
	
void SpectrumHistos::fill(float fgNorm, float bgNorm,
						  const std::vector<float>& fgPointsX, const std::vector<float>& fgPointsY,
						  const std::vector<float>& bgPointsX, const std::vector<float>& bgPointsY,
						  const std::string& hName,
						  unsigned int nbinsX, float xmin, float xmax,
						  unsigned int nbinsY, float ymin, float ymax) {
	
	setup2D(hName,nbinsX,xmin,xmax,nbinsY,ymin,ymax);
		
	// fill FG spectrum histograms
	for(UInt_t e=0; e<fgPointsX.size(); e++)
		fg2->Fill(fgPointsX[e],fgPointsY[e]);
	
	// fill BG spectrum histograms
	for(UInt_t e=0; e<bgPointsX.size(); e++)
		bg2->Fill(bgPointsX[e],bgPointsY[e]);

	finish2D(fgNorm,bgNorm);
}

void setLabels(TH1* h, const char* xlabel, const char* ylabel, const char* title) {
	if(xlabel)
		h->GetXaxis()->SetTitle(xlabel);
	if(ylabel)
		h->GetYaxis()->SetTitle(ylabel);
	if(title)
		h->SetTitle(title);
}

void SpectrumHistos::draw(bool logScale, const char* xlabel, const char* ylabel, const char* title) {
	
	defaultCanvas->cd();
	
	gPad->SetLogy(logScale);
	if(!logScale)
		fg->SetMinimum(diff->GetMinimum());	
	
	setLabels(fg,xlabel,ylabel,title);
	
	fg->Draw();
	bg->Draw("Same");
	diff->Draw("Same");
	printCanvas(chName);
	
	gPad->SetLogy(false);
}
