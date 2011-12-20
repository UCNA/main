#include "MWPC.hh"
#include "GraphicsUtils.hh"
#include <math.h>
#include <utility>

void MWPC::specialize() {

	if(mySide == EAST)
		setName("MWPC_E");
	else
		setName("MWPC_W");

	// anode
	if(thisRun->RI.runNum < 7000) {
		if(mySide == EAST)
			loadData("Pdc20","MWPCEAnode","anode");
		else
			loadData("Padc0","MWPCWAnode","anode");
	} else {
		if(mySide == EAST)
			loadData("Pdc30","MWPCEAnode","anode");
		else
			loadData("Pdc34","MWPCWAnode","anode");
	}
}

// function for determining wire pedestals
std::vector< std::pair<float,float> > anode_pedestal_finder(Subsystem* S, void*) {
	MWPC* M = (MWPC*)S;
	std::vector< std::pair<float,float> > v;
	for(unsigned int e=0; e<M->nEvents; e++)
		if(!M->BS->Trig->beta2of4(e,M->mySide) && !M->BS->Trig->isBeamnoise(e))
			v.push_back(std::pair<float, float>(M->BS->Trig->eventTime(e),M->getData("anode")[e]));
	return v;
}

float pedEdgeFinder(TH1F* h, Float_t width, TGraph*& g) {
	
	Float_t w,c;
	float n = 0;
	c = h->GetBinCenter(1);
	for(Int_t j = 2; j<=h->GetNbinsX(); j++) {
		w = (h->GetBinCenter(j) - c)/width;
		n += exp(-w*w);
	}
	n = 1.0/(2.0*n + 1.0);
	float* convolved = new float[h->GetNbinsX()];
	for(Int_t i = 1; i<=h->GetNbinsX(); i++) {
		convolved[i-1] = 0;
		c = h->GetBinCenter(i);
		for(Int_t j = 1; j<=h->GetNbinsX(); j++) {
			w = (h->GetBinCenter(j) - c)/width;
			if(h->GetBinContent(j) > 0)
				convolved[i-1] += n*exp(-w*w)*log(h->GetBinContent(j));
		}
	}
	for(int i=0; i<h->GetNbinsX(); i++)
		convolved[i] = exp(convolved[i]);
	
	// first local minimum finder
	int b0 = h->GetMaximumBin();
	while(convolved[b0] > convolved[b0+1] && convolved[b0] >= convolved[b0+1] )
		b0++;
	
	// precision minimum finder
	float maxval = (convolved[b0+1]-convolved[b0-1])/(4*convolved[b0]-2*(convolved[b0+1]+convolved[b0-1]));
	maxval =  (1-maxval)*h->GetBinCenter(b0+1) + maxval*h->GetBinCenter(b0+2);
	
	g = new TGraph(h->GetNbinsX());
	for(int i=0; i<h->GetNbinsX(); i++)
		g->SetPoint(i,h->GetBinCenter(i+1),convolved[i]);
	
	delete(convolved);
	return maxval;
}

MWPC::MWPC(RunManager* T, Side s, Wirechamber* x, Wirechamber* y, BetaScint* bs):
Subsystem(T,"",s), xPlane(x), yPlane(y), BS(bs) {
		
	specialize();

	if(!verifyPedestal(0))
		monitorPeak(anode_pedestal_finder, 60.0, 3000, NULL, sensorNames[0]);
	
	if(!(thisRun->runMode & RUNMODE_POSTRACK))
		return;
	
	// initialize events data
	events.resize(nEvents);
	
	//----------------------------
	// set triggers from wireplanes
	//----------------------------
	
	unsigned int ntriggers = 0;
	for(UInt_t e=0; e<nEvents; e++) {
		triggers[e] = xPlane->triggers[e] && yPlane->triggers[e];
		ntriggers += triggers[e];
	}
	printf("*** Found %i MWPC triggers (%f%%)\n",ntriggers,100.0*float(ntriggers)/T->nEvents);
	
	//----------------------------
	// cathode sum, anode, source ID
	//----------------------------

	float* f_anode = getData("anode");
	for(UInt_t i=0; i<nEvents; i++) {
		events[i].cathodeSum = xPlane->hits[i].cathodeSum + yPlane->hits[i].cathodeSum;
		events[i].anode = f_anode[i];
		events[i].sourceID = 0;
	}
	
	//----------------------------
	// calc anode cut
	//----------------------------
	
	TH1F* hAnode = registeredTH1F("Anode","MWPC Anode",204,-20,4000);
	for(UInt_t e=0; e<nEvents; e++)
		hAnode->Fill(events[e].anode);
	TGraph* gAnode = NULL;
	anodeCut = pedEdgeFinder(hAnode,10,gAnode);
	
	
	//----------------------------
	// calc cathode cut
	//----------------------------
	
	TH1F* hCathode = registeredTH1F("Cathode","MWPC Cathode",200,-500,30000);
	for(UInt_t e=0; e<nEvents; e++)
		hCathode->Fill(events[e].cathodeSum);
	TGraph* gCathode = NULL;
	cathodeCut = pedEdgeFinder(hCathode,100,gCathode);
	
	//----------------------------
	// cuts on anode & cathode histograms
	//----------------------------
	
	TH1F* hAnodeCut = registeredTH1F("Anode_Cut","MWPC Anode Cut",204,-20,4000);
	hAnodeCut->SetLineColor(2);
	TH1F* hCathodeCut = registeredTH1F("Cathode_Cut","MWPC Cathode Cut",200,-500,30000);
	hCathodeCut->SetLineColor(2);
	for(UInt_t e=0; e<nEvents; e++) {
		if(triggers[e]) {
			hAnodeCut->Fill(events[e].anode);
			hCathodeCut->Fill(events[e].cathodeSum);
		}
	}
	hAnode->Scale(1.0/(BS->Trig->runTime()*hAnode->GetBinWidth(1)));
	hAnodeCut->Scale(1.0/(BS->Trig->runTime()*hAnodeCut->GetBinWidth(1)));
	
	defaultCanvas->SetLogy(true);
	hAnode->Draw();
	hAnodeCut->Draw("Same");
	printCanvas("Anode");
	hCathode->Draw();
	hCathodeCut->Draw("Same");
	printCanvas("Cathode");
	defaultCanvas->SetLogy(false);
	
	delete(gAnode);
	delete(gCathode);
	
	genHistograms();
	
}

void MWPC::genHistograms() {
	
	//----------------------------
	// hit positions
	//----------------------------
	TH2F* h1 = blankHisto("Position","Wirechamber Hits",400);
	for(UInt_t e=0; e<nEvents; e++)
		if(triggers[e])
			h1->Fill(xPlane->hits[e].center,yPlane->hits[e].center);
	h1->Draw("COL");
	drawWires();
	drawFiducialCuts();
	printCanvas("Hit_Locations");
}

std::vector<unsigned int> MWPC::cutEllipse(Source p, Float_t nsigma, bool* selected) const {
	std::vector<unsigned int> v;
	for(unsigned int i=0; i<nEvents; i++)
		if( (!selected || selected[i]) && 
		   (xPlane->hits[i].center-p.x)*(xPlane->hits[i].center-p.x)/(nsigma*nsigma*p.wx*p.wx) 
		   + (yPlane->hits[i].center-p.y)*(yPlane->hits[i].center-p.y)/(nsigma*nsigma*p.wy*p.wy) < 1.0)
			v.push_back(i);
	return v;
}

