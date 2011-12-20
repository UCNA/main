#include "MuonVeto.hh"

void MuonVeto::specialize() {
		
	if(thisRun->RI.runNum < 7000) {
		loadData("Qadc8","ADCBackingE","back_adc");		// Qadc EBack
		loadData("Qadc10","ADCBackingW");				// Qadc WBack
		
		loadData("Pdc217","ADCTopE","Etop_adc");		// Qadc ETop
		loadData("Padc16","ADCSideE","ESide_adc");		// Padc ESide
		loadData("Qadc9","ADCSideW","WSide_adc");		// Padc WTopSide
		
		nchan = 5;
		
		loadData("Tdc018","TDCBackingE","back_tdc");	// Tdc  EBack
		loadData("Tdc020","TDCBackingW");				// Tdc  WBack
		loadData("Tdc019","TDCTopE");					// Tdc ETop
	} else {
		loadData("Qadc8","ADCBackingE","back_adc");		// Qadc EBack
		loadData("Qadc10","ADCBackingW");				// Qadc WBack
		
		loadData("Pdc313","ADCTopE","ETop_adc");		// Qadc ETop
		loadData("Pdc315","ADCTopW","WTop_adc");		// Padc WTop
		
		loadData("Pdc313","TACDriftE","drift_tac");		// E drift TAC
		loadData("Pdc315","TACDriftW");					// W drift TAC
		
		nchan = 4;
		
		loadData("Tdc018","TDCBackingE","back_tdc");	// Tdc  EBack
		loadData("Tdc020","TDCBackingW");				// Tdc  WBack
		loadData("Tdc019","TDCTopE","ETop_tdc");		// Tdc ETop	
	}
}

// function for determining muon pedestals
std::vector< std::pair<float,float> > mv_pedestal_finder(Subsystem* S, void* channel) {
	MuonVeto* M = (MuonVeto*)S;
	unsigned int tn = *(unsigned int*)channel;
	std::vector< std::pair<float,float> > v;
	for(unsigned int e=0; e<M->nEvents; e++)
		if(!M->TG->isBeamnoise(e) && !M->TG->beta2of4(e,(Side)(tn%2)))
			v.push_back(std::pair<float,float>(M->TG->eventTime(e),M->getData("back_adc",tn)[e]));
	return v;
}

bool MuonVeto::hitBacking(UInt_t e, Side s) const {
	if(s==EAST || s==WEST)
		return 500 < getData("back_tdc",s)[e] && getData("back_tdc",s)[e] < 3800;
	return s == sideCombo(hitBacking(e, EAST), hitBacking(e, WEST));
}

MuonVeto::MuonVeto(RunManager* T, Trigger* tg): Subsystem(T,"MuVeto",NONE), TG(tg) {
	
	specialize();
	
	striggers[EAST] = striggers[WEST] = NULL;
	
	// check pedestals (and recreate if necessary)
	for (unsigned int c = 0; c<nchan; c++ )
		if(!verifyPedestal(c))		
			monitorPeak(mv_pedestal_finder, 60.0, 3000, (void*)&c, sensorNames[c]);
	
	if(!(thisRun->runMode & RUNMODE_FULLID))
		return;
	
	// pedestal-subtract data
	for(unsigned int e=0; e<nEvents; e++)
		for(unsigned int c=0; c<nchan; c++)
			getData("back_adc",c)[e] -= PC.getPedestal(sensorNames[c],tg->eventTime(e));
				
	striggers[EAST] = new bool[nEvents];
	striggers[WEST] = new bool[nEvents];
	int ntriggers = 0;
	for(UInt_t e=0; e<nEvents; e++) {
		for(Side s = EAST; s <= WEST; s=nextSide(s)) {
			// backing veto, by TDC
			striggers[s][e] = hitBacking(e,s);
			// drift tubes, by TAC
			striggers[s][e] |= 500 < getData("drift_tac",s)[e] && getData("drift_tac",s)[e] < 3500;
		}
		// top veto
		striggers[EAST][e] |= 200 < getData("ETop_tdc")[e] && getData("ETop_tdc")[e] < 3800;
	}
	if(thisRun->runMode & RUNMODE_MAKEPLOTS)
		genHistograms();
	
	printf("*** Found %i muon events (%.2f %%)\n\n",ntriggers,ntriggers*100.0/nEvents);
}




void MuonVeto::genHistograms() {
	
	defaultCanvas->SetLogy(true);
	
	TH1F* hTAC[2];
	TH1F* hBackTDC[2];
	for(Side s=EAST; s<=WEST; s=nextSide(s)) {
		
		TH1F* h1 = registeredTH1F(std::string("Backing_Veto_Spectrum_")+sideNames(s),"Backing Veto Spectrum",500,-200,4000);
		TH1F* h2 = registeredTH1F(std::string("Backing_Veto_Cut_Spectrum_")+sideNames(s),"Timed Spectrum",500,-200,4000);
		hBackTDC[s] = registeredTH1F(std::string("Backing_Veto_TDC_")+sideNames(s),"Backing Veto Timing",500,10,3900);
		hTAC[s] = registeredTH1F(std::string("Drift_TAC_")+sideNames(s),"Timed Spectrum",500,-200,4000);
		
		h2->SetLineColor(2);
		hBackTDC[s]->SetLineColor(2+2*s);
		hTAC[s]->SetLineColor(2+2*s);
		
		for(UInt_t e=0; e<nEvents; e++) {
			h1->Fill(getData("back_adc",s)[e]);
			hBackTDC[s]->Fill(getData("back_tdc",s)[e]);
			hTAC[s]->Fill(getData("drift_tac",s)[e]);
			if(hitBacking(e, s))
				h2->Fill(getData("back_adc",s)[e]);
		}
		
		hTAC[s]->Scale(1.0/(TG->runTime()*hTAC[s]->GetBinWidth(1)));
		hBackTDC[s]->Scale(1.0/(TG->runTime()*hBackTDC[s]->GetBinWidth(1)));
		h1->Scale(1.0/(TG->runTime()*h1->GetBinWidth(1)));
		h2->Scale(1.0/(TG->runTime()*h2->GetBinWidth(1)));
		
		h1->Draw();
		h2->Draw("Same");
		printCanvas(std::string("Backing_Spectrum_")+sideNames(s));
	}
	
	hTAC[EAST]->Draw();
	hTAC[WEST]->Draw("Same");
	printCanvas("Drift_TAC");
	
	hBackTDC[EAST]->Draw();
	hBackTDC[WEST]->Draw("Same");
	printCanvas("Backing_TDC");
	
	defaultCanvas->SetLogy(false);
}
