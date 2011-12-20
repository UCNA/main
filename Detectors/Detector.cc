#include "Detector.hh"
#include <TSpectrum.h>
#include <TSpectrum2.h>
#include "Source.hh"
#include "TSpectrumUtils.hh"
#include "MultiGaus.hh"
#include "ReSource.hh"
#include "PMTGenerator.hh"
#include "KurieFitter.hh"

Detector::Detector(RunManager* T, Trigger* tg, MWPC* mwpce, MWPC* mwpcw, BetaScint* bse, BetaScint* bsw, MuonVeto* mv):
EventParser(NULL,NULL), Subsystem(T,"Main",NONE), muVeto(mv), TG(tg), events(new BetaHit[T->nEvents]) {
	
	if(!(thisRun->runMode & RUNMODE_FULLID && thisRun->runMode & RUNMODE_POSTRACK))
		return;
	
	idflags[EAST].resize(nEvents);
	idflags[WEST].resize(nEvents);
	setFlagsReadPoint(&*idflags[EAST].begin(),&*idflags[WEST].begin());
	MWPCs[EAST] = mwpce;
	MWPCs[WEST] = mwpcw;
	Scints[EAST] = bse;
	Scints[WEST] = bsw;
	
	printf("Combining detector subsystems...\n");
	idEvents();
	calibrationSpectra();
	if(thisRun->runMode & RUNMODE_MAKEPLOTS)
		genHistograms();
}

void Detector::idEvents() {
	printf("Identifying events..."); fflush(stdout);
	//counts.rtime = TG->runTime();
	for(UInt_t e=0; e<nEvents; e++) {
		
		if(e%(nEvents/20)==0) { printf("*"); fflush(stdout); }
		
		events[e].sourceID = 0;
		
		// fill side hit characteristics
		for(Side s = EAST; s <= WEST; s=nextSide(s)) {
			idflags[s][e] = 0;
			float x = MWPCs[s]->xPlane->hits[e].center;
			float y = MWPCs[s]->yPlane->hits[e].center;
			thisRun->PCal->calibrateEnergy(s,x,y,Scints[s]->events[e],TG->eventTime(e));
			
			if(MWPCs[s]->triggers[e])
				idflags[s][e] |= HIT_WIRES;
			if(Scints[s]->triggers[e])
				idflags[s][e] |= HIT_SCINT;
			if(muVeto->muonTrigger(e))
				idflags[s][e] |= IS_MUON;
			if(muVeto->hitBacking(e,s))
				idflags[s][e] |= HIT_MUON_BACK;
			if(TG->primaryScint(e,s))
				idflags[s][e] |= HIT_SIDE_FIRST;
			if(TG->isCrud(e))
				idflags[s][e] |= IS_BEAMCRUD;
			if(TG->isUCNMon(e))
				idflags[s][e] |= IS_UCN_MON;
			if(TG->gmsTrigger(e))
				idflags[s][e] |= IS_GMS;
			//if(Scints[s]->isPulserTrigger(e))
			if(TG->pulserTrigger(e))
				idflags[s][e] |= IS_PMT_PULSER;
		}
		
		// classify events and fill events data
		classifyEventType(e);
		for(Side s = EAST; s <= WEST; s = nextSide(s))
			events[e].typeFlags[s] = idflags[s][e];
		
	}
	
	printf(" Done.\n");
	//counts.print();	// TODO event type counter display
	//parent->qOut.insert("EventCounts",counts.toStringmap());
}



// allocate a 2D float** array
float** allocArray(unsigned int nx, unsigned int ny) {
	float** a = new float*[nx];
	for(unsigned int i=0; i<nx; i++)
		a[i] = new float[ny];
	return a;
}

// delete a 2D float** array
void deleteArray(float** a, unsigned int nx) {
	for(unsigned int i=0; i<nx; i++)
		delete[](a[i]);
	delete[](a);
}

// comparison to sort sources by x position in descending order
bool sortSourceXPos(Source a, Source b) {
	return a.x < b.x;
}

void Detector::calibrationSpectra() {
	
	if(thisRun->RI.type != SOURCE_RUN || !(thisRun->SI.sourceHolder[EAST].size() || thisRun->SI.sourceHolder[WEST].size()))
		return;
	
	// temporary arrays
	bool* correctEvents = new bool[nEvents];
	const unsigned int nbins = 200;
	float** zout = allocArray(nbins,nbins);
	float** historray = allocArray(nbins,nbins);
	
	
	for(Side s = EAST; s <= WEST; s = nextSide(s)) {
		
		if(!thisRun->SI.sourceHolder[s].size())
			continue;
		
		//---------
		// locate sources
		//---------
		
		printf("Searching wirechamber for source positions...\n");
		
		// create and fill hit position histogram
		TH2F* hitpos = MWPCs[s]->blankHisto("HitPositions","HitPositions",nbins);
		for(UInt_t e=0; e<nEvents; e++) {
			if( (correctEvents[e] = (isType0Beta(s,e)
									 && MWPCs[s]->xPlane->hits[e].errflags == WIRES_GOOD 
									 && MWPCs[s]->yPlane->hits[e].errflags == WIRES_GOOD) ) ) {
				float x0 = MWPCs[s]->xPlane->hits[e].center;
				float y0 = MWPCs[s]->yPlane->hits[e].center;
				hitpos->Fill(x0,y0);
			}
		}
				
		// convert hit position histogram to float** array
		for(unsigned int x=0; x<nbins; x++)
			for(unsigned int y=0; y<nbins; y++)
				historray[x][y] = hitpos->GetBinContent(x+1,y+1);
		
		// run TSpectrum 2D peak-finding algorithm
		TSpectrum2 TS = TSpectrum2();
		unsigned int npks = TS.SearchHighRes(
											 historray,		// input spectrum
											 zout,			// output deconvolved spectrum
											 nbins,nbins,	// x-y dimensions
											 5.0,			// sigma
											 0.2,			// threshold in percent of highest peak
											 true,			// remove background?
											 1,				// deconvolution iterations
											 false,			// regenerate input spectrum using Markov chains
											 0				// averaging window for markov chains
											 );
		
		printf("\t%i of %i sources found.\n",npks,(int)thisRun->SI.sourceHolder[s].size());
		
		// extract located peaks into source positions vector
		// drop peaks past expected number of sources
		Float_t* xPos = TS.GetPositionX();
		Float_t* yPos = TS.GetPositionY();
		std::vector<Source> foundSources;
		for(unsigned int i=0; i<npks; i++) {
			if(i>=thisRun->SI.sourceHolder[s].size())
				break;
			Source pk;
			pk.x = binterpolate(hitpos->GetXaxis(),xPos[i]);
			pk.wx = 2.0;
			pk.y = binterpolate(hitpos->GetYaxis(),yPos[i]);
			pk.wy = 2.0;
			pk.mySide = s;
			foundSources.push_back(pk);
		}
		
		
		// sort sources and apply source types from source holder listing
		std::sort(foundSources.begin(),foundSources.end(),&sortSourceXPos);
		for(unsigned int i=0; i < foundSources.size(); i++)
			foundSources[i].t = thisRun->SI.sourceHolder[s][i];
		
		// refine source locations and widths, enter into run sources list
		for(std::vector<Source>::iterator it = foundSources.begin(); it!=foundSources.end(); it++) {
			std::vector<unsigned int> srcPts = MWPCs[s]->cutEllipse(*it,3.0,correctEvents);
			// collect points and fit
			std::vector<float> srcX(srcPts.size());
			std::vector<float> srcY(srcPts.size());
			for(UInt_t e=0; e<srcPts.size(); e++) {
				srcX[e] =  MWPCs[s]->xPlane->hits[srcPts[e]].center;
				srcY[e] =  MWPCs[s]->yPlane->hits[srcPts[e]].center;
			}
			SpectrumPeak px = gausFitter(srcX, it->x, it->wx);
			SpectrumPeak py = gausFitter(srcY, it->y, it->wy);
			it->x = px.center.x;
			it->wx = (px.width.x < 4.0)?px.width.x : 4.0;
			it->y = py.center.x;
			it->wy = (py.width.x < 4.0)?py.width.x : 4.0;			
		}
		
		if(!foundSources.size())
			warn(SEVERE_WARNING,"Missing_All_Sources");
		else if(foundSources.size() < thisRun->SI.sourceHolder[s].size())
			warn(MODERATE_WARNING,"Missing_Some_Sources");
		
		//-------------
		// fit peaks
		//-------------
		
		// redefine "good" events with less stringent positioning quality cuts
		for(UInt_t e=0; e<nEvents; e++)
			correctEvents[e] = isType0Beta(s,e);
		
		for(std::vector<Source>::iterator it = foundSources.begin(); it!=foundSources.end(); it++) {
			// cut source events by position
			std::vector<unsigned int> srcPts = MWPCs[s]->cutEllipse(*it,4.0);
			it->nCounts = srcPts.size();
			Source thisSource = thisRun->SI.addSource(s,*it);
			printf("Found %i points in wirechamber source cut\n",(int)srcPts.size());
			thisSource.nCounts = 0;
			for(std::vector<unsigned int>::const_iterator pit = srcPts.begin(); pit != srcPts.end(); pit++) {
				if(isAnyBeta(EAST,*pit) || isAnyBeta(WEST,*pit)) { // TODO && whichSideFirst(*pit)) {
					events[*pit].sourceID = thisSource.sID;
					thisSource.nCounts++;
				}
				MWPCs[s]->events[*pit].sourceID = thisSource.sID;
			}
			
			thisSource.display();
		}
		
	}
	
	// cleanup and return
	delete[](correctEvents);
	deleteArray(historray,nbins);
	deleteArray(zout,nbins);
	
	
	return;
}




void Detector::draw2Dcuts(Side s) {
	drawFiducialCuts();
	for(std::vector<Source>::const_iterator it = thisRun->SI.sourcesBegin(s); it != thisRun->SI.sourcesEnd(s); it++)
		drawEllipseCut(*it,4.0,it->name());
}

void Detector::genHistograms() {
	
	float sideCounts[2];
	
	for(Side m = EAST; m <= WEST; m = nextSide(m)) {
		
		//----------------------------
		// Beta hits and radial distribution
		//----------------------------
		TH2F* h1 = MWPCs[m]->blankHisto("Beta_Hits","Beta Hits");
		TH1F* hBetaR = registeredTH1F(std::string("Beta_Radius_")+sideNames(m),"Beta Hits Radius^2 Distribution",50,0,90*90*0.6);
		for(UInt_t e=0; e<nEvents; e++) {
			if(isType0Beta(m,e)) {
				float x0 = MWPCs[m]->xPlane->hits[e].center;
				float y0 = MWPCs[m]->yPlane->hits[e].center;
				h1->Fill(x0,y0);
				hBetaR->Fill(MWPCs[m]->radius(e)*MWPCs[m]->radius(e));
			}
		}
		h1->Draw("COL");
		draw2Dcuts(m);
		MWPCs[m]->drawWires();
		printCanvas(std::string("Beta_Hits_")+sideNames(m));
		hBetaR->Draw();
		for(float x=0; x<100; x+=10)
			drawVLine(x*x,defaultCanvas);
		drawVLine(wirechamber_window_radius*wirechamber_window_radius,defaultCanvas,2);
		printCanvas(std::string("Beta_Radial_")+sideNames(m));
		
		//----------------------------
		// muon hits
		//----------------------------
		
		TH2F* h2 =  MWPCs[m]->blankHisto("Muon_Hits","Muon Hits");
		for(UInt_t e=0; e<nEvents; e++) {
			if(isMuon(e) && hitSide(m,e) && muVeto->hitBacking(e,m)) {
				h2->Fill(MWPCs[m]->xPlane->hits[e].center,MWPCs[m]->yPlane->hits[e].center);
			}
		}
		h2->Draw("COL");
		draw2Dcuts(m);
		MWPCs[m]->drawWires();
		printCanvas(std::string("Muon_Hits_")+sideNames(m));
		
		//----------------------------
		// non-backing-veto muon hits
		//----------------------------
		
		TH2F* h3 =  MWPCs[m]->blankHisto("Nonbacking_Muon_Hits","Nonbacking Muon Hits");
		for(UInt_t e=0; e<nEvents; e++) {
			if(isMuon(e) && hitSide(m,e) && !muVeto->hitBacking(e,m)) {
				h3->Fill(MWPCs[m]->xPlane->hits[e].center,MWPCs[m]->yPlane->hits[e].center);
			}
		}
		h3->Draw("COL");
		draw2Dcuts(m);
		MWPCs[m]->drawWires();
		printCanvas(std::string("Muon_Nonbacking_Hits_")+sideNames(m));
		
		//----------------------------
		// QADC, Energy spectra with simple cuts
		//----------------------------
		TH1F* hBeta = registeredTH1F(std::string("Beta_Spectrum_")+sideNames(m),"Beta Events Spectrum",100,0,2000);
		TH1F* hTubens[nBetaTubes];
		TH1F* hTubensCut[nBetaTubes];
		TH1F* hPMTs[nBetaTubes];
		TH1F* hPMTsCut[nBetaTubes];
		for(unsigned int t=0; t<nBetaTubes; t++) {
			hTubens[t] = registeredTH1F(std::string("Tube_Spectrum_")+sideNames(m)+itos(t),"Tube Energy Spectrum",300,-100,2000);
			hTubens[t]->SetLineColor(t+1);
			hTubensCut[t] = registeredTH1F(std::string("Tube_Beta_Spectrum_")+sideNames(m)+itos(t),"Tube Energy Spectrum",300,-100,2000);
			hTubensCut[t]->SetLineColor(t+1);
			hPMTs[t] = registeredTH1F(std::string("PMT_QADC_")+sideNames(m)+itos(t),"Beta PMT QADC",100,-100,4096);
			hPMTsCut[t] = registeredTH1F(std::string("PMT_QADC_Cut_")+sideNames(m)+itos(t),"Beta PMT QADC, cut",100,-100,4096);
			hPMTsCut[t]->SetLineColor(2);
		}
		ScintEvent evt;
		
		hBeta->Sumw2();
		for(UInt_t e=0; e<nEvents; e++) {
			
			float x = MWPCs[m]->xPlane->hits[e].center;
			float y = MWPCs[m]->yPlane->hits[e].center;
			for(unsigned int t=0; t<nBetaTubes; t++)
				evt.adc[t] = Scints[m]->tubedat[t][e];
			thisRun->PCal->calibrateEnergy(m,x,y,evt,TG->eventTime(e));
			
			for(unsigned int t=0; t<nBetaTubes; t++) {
				hPMTs[t]->Fill(evt.adc[t]);
				if(TG->betaTrigger(e,m,t))
					hTubens[t]->Fill(evt.tuben[t].x);
			}
			
			if(muVeto->muonTrigger(e))
				continue;
			if(isType0Beta(m,e) && MWPCs[m]->radius(e)<60.0) {
				hBeta->Fill(evt.energy.x);
				for(unsigned int t=0; t<nBetaTubes; t++) {
					hTubensCut[t]->Fill(evt.tuben[t].x);
					hPMTsCut[t]->Fill(evt.adc[t]);
				}
			}
		}
		hBeta->Draw();
		printCanvas(std::string("Beta_Spectrum_")+sideNames(m));
		sideCounts[m] = hBeta->Integral();
		
		defaultCanvas->SetLogy(true);
		for(unsigned int t=0; t<nBetaTubes; t++) {
			hPMTs[t]->Draw();
			hPMTsCut[t]->Draw("Same");
			printCanvas(std::string("PMT_QADC_")+sideNames(m)+itos(t));
		}
		hTubens[nBetaTubes-1]->Draw();
		for(unsigned int t=0; t<nBetaTubes; t++) {
			if(t!=nBetaTubes-1)
				hTubens[t]->Draw("SAME");
			hTubensCut[t]->Draw("SAME");
		}
		printCanvas(std::string("Tube_Spectra_")+sideNames(m));
		defaultCanvas->SetLogy(false);
		
		// Kurie plot for beta decay runs
		if(thisRun->RI.type == ASYMMETRY_RUN && thisRun->RI.gvState == GV_OPEN) {
			printf("Generating kurie plot from beta spectrum...\n");
			TGraphErrors* tgData = NULL;
			float_err betaEndpoint = kurieIterator(hBeta,neutronBetaEp,&tgData);
			if(tgData) {
				tgData->SetMarkerColor(2);
				tgData->Draw("AP");
				printCanvas(std::string("Kurie_Plot_")+sideNames(m));
				delete(tgData);
			}
			forParent.insert(std::string("Kurie_")+sideNames(m),betaEndpoint.x);
			forParent.insert(std::string("Kurie_")+sideNames(m)+"_err",betaEndpoint.err);
		} else {
			printf("No kurie plotting for %i-%i runs.\n",thisRun->RI.type,thisRun->RI.gvState);
		}
		
		// event spacing
		if(false) {
			TH1F* hDeltaT = registeredTH1F(std::string("Delta_t_")+sideNames(m),"Event Spacing (microseconds)",200,-200,5000);
			TH1F* hDeltaTw = registeredTH1F(std::string("Delta_t_Events_")+sideNames(m),"Event Spacing (microseconds)",200,-200,5000);
			const float* f_delt0 = TG->getData("Delt0");
			for(unsigned int e=0; e<nEvents; e++) {
				if(TG->isDataEvent(e) && TG->beta2of4Side(e) == m) {
					hDeltaT->Fill(f_delt0[e]);
					if(TG->triggers[e])
						hDeltaTw->Fill(f_delt0[e]);
				}
			} 
			defaultCanvas->SetLogy(true);
			hDeltaTw->SetLineColor(2);
			hDeltaT->Draw();
			hDeltaTw->Draw("Same");
			printCanvas(std::string("Event_Spacing_")+sideNames(m));
			defaultCanvas->SetLogy(false);
		}
	}
	
	
	
	//----------------------------
	// event rate histograms
	//----------------------------
	TH1F* hBetaE = registeredTH1F("BetaE","E Type 0 Rate",200,0,TG->totalTime());
	TH1F* hBetaW = registeredTH1F("BetaW","W Type 0 Rate",200,0,TG->totalTime());
	TH1F* hType1E = registeredTH1F("TypeIE","E Type 1 Rate",200,0,TG->totalTime());
	TH1F* hType1W = registeredTH1F("TypeIW","W Type 1 Rate",200,0,TG->totalTime());
	TH1F* hType2E = registeredTH1F("TypeIIE","E Type 2 Rate",200,0,TG->totalTime());
	TH1F* hType2W = registeredTH1F("TypeIIW","W Type 2 Rate",200,0,TG->totalTime());
	TH1F* hMuon = registeredTH1F("Muons","Muon Rate",200,0,TG->totalTime());
	
	hBetaE->SetLineColor(2);		// red
	hBetaW->SetLineColor(4);		// blue
	hType1E->SetLineColor(3);		// green
	hType1W->SetLineColor(5);		// cyan
	hType2E->SetLineColor(6);		// magenta
	hType2W->SetLineColor(7);		// yellow
	hMuon->SetLineColor(1);			// black
	
	for(unsigned int e=0; e<nEvents; e++) {
		if(isType0Beta(EAST,e) && MWPCs[EAST]->radius(e) < 50)
			hBetaE->Fill(TG->eventTime(e));
		if(isType0Beta(WEST,e) && MWPCs[WEST]->radius(e) < 50)
			hBetaW->Fill(TG->eventTime(e));
		if(isTypeIBeta(EAST,e) && MWPCs[EAST]->radius(e) < 50)
			hType1E->Fill(TG->eventTime(e));
		if(isTypeIBeta(WEST,e) && MWPCs[WEST]->radius(e) < 50)
			hType1W->Fill(TG->eventTime(e));
		if(isTypeII_IIIBeta(EAST,e) && MWPCs[EAST]->radius(e) < 50)
			hType2E->Fill(TG->eventTime(e));
		if(isTypeII_IIIBeta(WEST,e) && MWPCs[WEST]->radius(e) < 50)
			hType2W->Fill(TG->eventTime(e));	
		if(isMuon(e))
			hMuon->Fill(TG->eventTime(e));
	}
	hBetaE->Scale(1.0/hBetaE->GetBinWidth(1));
	hBetaW->Scale(1.0/hBetaW->GetBinWidth(1));
	hType1E->Scale(1.0/hType1E->GetBinWidth(1));
	hType1W->Scale(1.0/hType1W->GetBinWidth(1));
	defaultCanvas->SetLogy(true);
	hBetaE->SetMinimum(1e-1);
	hBetaE->SetMaximum(1e2);
	hBetaE->Draw();
	hBetaW->Draw("Same");
	hType1E->Draw("Same");
	hType1W->Draw("Same");
	hType2E->Draw("Same");
	hType2W->Draw("Same");
	hMuon->Draw("Same");
	printCanvas("EventRates");
	defaultCanvas->SetLogy(false);
	
	
	float rate[2];
	rate[EAST] = sideCounts[EAST]/TG->runTime(EAST);
	rate[WEST] = sideCounts[WEST]/TG->runTime(WEST);
	
	/*
	printf("\n\n*************************\n");
	printf("%i E counts, %i W counts\n",int(sideCounts[EAST]),int(sideCounts[WEST]));
	printf("%.2f hours running\n",TG->runTime()/3600.0);
	printf("%.2fHz E, %.2fHz W => %.2fHz total\n",rate[EAST],rate[WEST],rate[EAST]+rate[WEST]);
	printf("(E-W)/(E+W) Asymmetry = %.2f%%\n",100.0*(rate[EAST]-rate[WEST])/(rate[EAST]+rate[WEST]));
	printf(    "*************************\n\n");
	 */
}

void Detector::addOutBranches(TTree* T) { T->Branch(name.c_str(),&anEvent,"flagsE/I:flagsW/I:sID/I"); }
void Detector::fillEvent(UInt_t e) { anEvent = events[e]; }
