#include "BetaScint.hh"
#include "TSpectrumUtils.hh"
#include "GraphicsUtils.hh"
#include <TPolyMarker.h>
#include <TGraphAsymmErrors.h>
#include "MultiGaus.hh"

void BetaScint::specialize() {
	
	if(mySide == EAST)
		setName("BetaScE");
	else
		setName("BetaScW");
	
	// PMTs ordered by quadrant
	if(mySide == EAST) {
		loadData("Qadc2","ADCE3Beta","tube");	// E0 (=JE3)
		loadData("Qadc3","ADCE4Beta");			// E1 (=JE4)
		loadData("Qadc0","ADCE1Beta");			// E2 (=JE1)
		loadData("Qadc1","ADCE2Beta");			// E3 (=JE2)
		
		if(thisRun->RI.runNum < 7000)
			loadData("Pdc21","ADCRefELED","reftube");		// 4 GMS reference phototube
		else if(thisRun->RI.runNum < 13000)
			loadData("Pdc32","ADCRefELED","reftube");		// 4 GMS reference phototube
	} else {
		loadData("Qadc4","ADCW1Beta","tube");	// W0
		loadData("Qadc5","ADCW2Beta");			// W1
		loadData("Qadc6","ADCW3Beta");			// W2
		loadData("Qadc7","ADCW4Beta");			// W3
		
		if(thisRun->RI.runNum < 7000)
			loadData("Pdc216","ADCRefWLED","reftube");	// 4 GMS reference phototube
		else if(thisRun->RI.runNum < 13000)		
			loadData("Pdc36","ADCRefWLED","reftube");	// 4 GMS reference phototube		
	}
	
	loadData("Pdc36","","refpd");				// reference LED photodiode for late runs
}

// function for determining beta scintillator pedestals
std::vector< std::pair<float,float> > beta_pedestal_finder(Subsystem* S, void* tb) {
	BetaScint* BS = (BetaScint*)S;
	unsigned int t = *(unsigned int*)tb;
	std::vector< std::pair<float,float> > v;
	// look for pedestals in data events not on this side of detector or firing PMT (TODO reconsider this)
	for(unsigned int e=0; e<BS->nEvents; e++)
		if(BS->Trig->isDataEvent(e) && !BS->Trig->isBeamnoise(e) && !BS->Trig->pulserTrigger(e)
		   && !BS->Trig->beta2of4(e,BS->mySide) && !BS->Trig->betaTrigger(e,BS->mySide,t) )
			v.push_back(std::pair<float,float>(BS->Trig->eventTime(e),BS->tubedat[t][e]));
	return v;
}

// function for determining GMS LED events
std::vector< std::pair<float,float> > beta_LED_finder(Subsystem* S, void* tdat) {
	BetaScint* BS = (BetaScint*)S;
	float* tubedat = (float*)tdat;
	std::vector< std::pair<float,float> > v;
	for(unsigned int e=0; e<BS->nEvents; e++)
		if(BS->Trig->gmsLED(e) && !BS->Trig->isBeamnoise(e))
			v.push_back(std::pair<float,float>(BS->Trig->eventTime(e),tubedat[e]));
	return v;
}

BetaScint::BetaScint(RunManager* T, Trigger* TG, Side s): Subsystem(T,"",s), events(new ScintEvent[T->nEvents]), Trig(TG) {
	
	specialize();
	for(unsigned int t=0; t<nBetaTubes; t++)
		tubedat.push_back(getData("tube",t));
	
	// determine triggers
	int ntriggers = 0;
	for(UInt_t e=0; e<nEvents; e++)
		if( (triggers[e] = Trig->beta2of4(e, mySide)) )
			ntriggers++;
	
	for(unsigned int t=0; t<nBetaTubes+(thisRun->RI.runNum<13000?1:0); t++)
		if(!verifyPedestal(t))
			monitorPeak(beta_pedestal_finder, 60.0, 3000, (void*)&t, sensorNames[t]);
	
	if(!(thisRun->runMode & RUNMODE_FULLID || thisRun->runMode & RUNMODE_LEDTRACK))
		return;
	
	// pre-subtract pedestals
	for(UInt_t t=0; t<nBetaTubes+(thisRun->RI.runNum<13000?1:0); t++)
		for(UInt_t e=0; e<nEvents; e++)
			tubedat[t][e] -= PC.getPedestal(sensorNames[t],TG->eventTime(e));
	
	// check for LED peaks
	/*
	 //TODO
	 for(unsigned int t=0; t<5; t++) {
	 TGraph* g = thisRun->PCal->CDB->getLED(thisRun->RI.runNum,sensorNames[t]);
	 if(!g) {
	 if(!thisRun->manualAnalyze) {
	 warn(SEVERE_WARNING,sensorNames[t]+"_LED_MISSING_ABORT");
	 thisRun->runMode = RUNMODE_LEDTRACK;
	 } else {
	 warn(BENIGN_WARNING,sensorNames[t]+"_LED_Missing");
	 }
	 monitorPeak(beta_LED_finder, 60.0, 1000, (void*)tubedat[t], sensorNames[t]+"_led", false);
	 } else
	 delete(g);
	 }
	 */
	
	if(!(thisRun->runMode & RUNMODE_FULLID))
		return;
	
	for(UInt_t t=0; t<nBetaTubes; t++)
		qMax[t] = 3999.0 - PC.getPedestal(sensorNames[t],0);
	
	// fill ped-subtracted tube event data
	for(unsigned int e=0; e<nEvents; e++)
		for(unsigned int t=0; t<nBetaTubes; t++)
			events[e].adc[t] = tubedat[t][e];
	
	triggerThresholds();
	gmsPlots();
	
	printf("*** Found %i Beta scintillator events (%.2f %%)\n\n",ntriggers,ntriggers*100.0/nEvents);
}

void BetaScint::addOutBranches(TTree* T) { 
	T->Branch(name.c_str(),&anEvent,"q1/F:q2/F:q3/F:q4/F:e1/F:de1/F:e2/F:de2/F:e3/F:de3/F:e4/F:de4/F:energy/F:denergy/F:nPE1/F:nPE2/F:nPE3/F:nPE4/F");
	T->Branch((name+"_led_pd").c_str(),&led_pd,(name+"_led_pd/F").c_str());
}

void BetaScint::fillEvent(UInt_t e) { 
	anEvent = events[e];
	led_pd = getData("refpd")[e];
}

Stringmap BetaScint::finalWords() {
	for(unsigned int t=0; t<nBetaTubes; t++) {
		forParent.insert("gms0_"+itos(t),thisRun->PCal->getGMS0(mySide,t));
		forParent.insert("gmsRel_"+itos(t),thisRun->PCal->gmsFactor(mySide,t,0)/thisRun->PCal->getGMS0(mySide,t));
	}
	return forParent;
}

bool BetaScint::isPulserTrigger(unsigned int e) {
	if(Trig->pulserTrigger(e))
		return true;
	unsigned int nthresh = 0;
	unsigned int nhigh = 0;
	for(unsigned int t=0; t<nBetaTubes; t++) {
		nthresh += tubedat[t][e] > 200;
		nhigh += tubedat[t][e] > 1500;
	}
	return nhigh == 1 && nthresh == 1;
}


void BetaScint::gmsPlots() {
	
	// approximate PMT pedestals
	float tubePeds[nBetaTubes];
	for(unsigned int t=0; t<nBetaTubes; t++)
		tubePeds[t] = PC.getPedestal(sensorNames[t],0);
	
	//----------------
	// Draw GMS reference tube histograms
	//----------------
	
	if(thisRun->RI.runNum < 13000) {
		TH1F* hRef = registeredTH1F("gmsRef","GMS Reference PMT",450,-500,4000);
		hRef->Sumw2();
		TH1F* hRefCoFull = registeredTH1F("gmsRefCoFull","GMS Co60 Reference",450,-500,4000);
		hRefCoFull->Sumw2();
		TH1F* hRefCo = registeredTH1F("gmsRefCo","GMS Co60 Reference",200,1200,3500);
		hRefCo->Sumw2();
		TH1F* hRefLED = registeredTH1F("gmsRefLED","GMS LED Reference",450,-500,4000);
		hRefLED->Sumw2();
		hRefCoFull->SetLineColor(4);
		hRefLED->SetLineColor(2);
		for(UInt_t e=0; e<nEvents; e++) {
			hRef->Fill(tubedat[nBetaTubes][e]);
			if( Trig->gmsCo(e,mySide) ) {
				hRefCo->Fill(tubedat[nBetaTubes][e]);
				hRefCoFull->Fill(tubedat[nBetaTubes][e]);
			}
			if(Trig->gmsLED(e))
				hRefLED->Fill(tubedat[nBetaTubes][e]);
		}
		hRef->Scale(1.0/(Trig->runTime(mySide)*hRef->GetBinWidth(1)));
		hRefCo->Scale(1.0/(Trig->runTime(mySide)*hRefCo->GetBinWidth(1)));
		hRefCoFull->Scale(1.0/(Trig->runTime(mySide)*hRefCoFull->GetBinWidth(1)));
		hRefLED->Scale(1.0/(Trig->runTime(mySide)*hRefLED->GetBinWidth(1)));
		
		defaultCanvas->SetLogy(true);
		hRef->Draw();
		hRefCoFull->Draw("Same");
		hRefLED->Draw("Same");
		printCanvas("GMS_Reference_Spectrum");
		defaultCanvas->SetLogy(false);
	}
	
	
	
	//----------------
	// GMS LED peak fitting (for LED scans)
	//----------------
	if(false) {	
		for(unsigned int t=0; t<nBetaTubes+1; t++) {
			
			// collect GMS LED data
			TH1F* hLED = registeredTH1F(std::string("gmsLED_")+itos(t),"GMS LED",501,-500,4000);
			hLED->SetLineColor(2);
			hLED->Sumw2();
			
			for(UInt_t e=0; e<nEvents; e++)
				if(Trig->gmsLED(e)) // && tubedat[t][e]<3999-PC.getPedestal(sensorNames[t],Trig->eventTime(e)))
					hLED->Fill(tubedat[t][e]);
			printf("%c%i: %i GMS LED events\n",sideNames(mySide),t,(int)hLED->Integral());
			hLED->Scale(1.0/(Trig->runTime(mySide)*hLED->GetBinWidth(1)));
			
			// Draw histograms
			if(hLED->GetEntries()) {
				//defaultCanvas->SetLogy(true);
				hLED->Draw();
				printCanvas(std::string("GMS_LED_Spectrum_")+itos(t));
				//defaultCanvas->SetLogy(false);
			}
			
		}
	}
	
	//----------------
	// Chris Pulser peak fits
	//----------------
	unsigned int nDivisions = (unsigned int)(Trig->totalTime()/300.0);
	nDivisions = nDivisions?nDivisions:1;
	printf("Scanning Chris Pulser peaks (%i divisions)...",nDivisions); fflush(stdout);
	TH1F* pulspectra[nBetaTubes];
	std::vector<TH1*> allSpectra[nBetaTubes];
	unsigned int divn = 0;
	float tavg=0;
	float twt=0;
	TF1 gausFit("gasufit","gaus",1000,4000);
	QFile pulseLocation(dataPath+"ChrisPulser.txt",false);
	for(unsigned int e=0; e<nEvents; e++) {
		if(e>=(nEvents*divn)/nDivisions || e==nEvents-1) {
			printf("*"); fflush(stdout);
			// fit old histograms
			if(divn>0) {
				for(unsigned int t=0; t<nBetaTubes; t++) {
					float bcenter = pulspectra[t]->GetBinCenter(pulspectra[t]->GetMaximumBin());
					Stringmap m;
					m.insert("side",ctos(sideNames(mySide)));
					m.insert("tube",t);
					m.insert("time",tavg/twt);
					m.insert("abstime",tavg/twt+thisRun->RI.startTime);
					m.insert("counts",pulspectra[t]->GetEntries());
					if(!iterGaus(pulspectra[t],&gausFit,3,bcenter,200,1.5)) {
						m.insert("height",gausFit.GetParameter(0));
						m.insert("center",gausFit.GetParameter(1));
						m.insert("width",gausFit.GetParameter(2));
						m.insert("dheight",gausFit.GetParError(0));
						m.insert("dcenter",gausFit.GetParError(1));
						m.insert("dwidth",gausFit.GetParError(2));
					}
					pulseLocation.insert("pulserpeak",m);
				}
			}
			divn++;
			if(divn>nDivisions)
				break;			
			// setup new histograms
			for(unsigned int t=0; t<nBetaTubes; t++) {
				tavg = twt = 0;
				pulspectra[t] = registeredTH1F(std::string("Pulspectrum_")+itos(t)+"_"+itos(divn),"Pulser Spectrum",100,1700-tubePeds[t],4000-tubePeds[t]);
				pulspectra[t]->SetLineColor(t+1);
				allSpectra[t].push_back(pulspectra[t]);
			}
		}
		if(!isPulserTrigger(e))
			continue;
		for(unsigned int t=0; t<nBetaTubes; t++)
			pulspectra[t]->Fill(tubedat[t][e]);
		tavg += Trig->eventTime(e);
		twt += 1.0;
	}
	printf(" Done\n");
	for(unsigned int t=0; t<nBetaTubes; t++) {
		drawSimulHistos(allSpectra[t]);
		printCanvas(std::string("Chris_Pulser_Spectra_")+itos(t));
	}
	pulseLocation.commit();
}


void BetaScint::triggerThresholds(bool allPlots, const std::vector<bool>& otherCuts) {
	for(unsigned int t=0; t<nBetaTubes; t++) {
		char tmp[1024];
		sprintf(tmp,"%c%i Trigger Efficiency",sideNames(mySide),t+1);
		TH1F hAll("hAll",tmp,125,-50,200);
		TH1F hTrig("hTrig",tmp,125,-50,200);
		for(unsigned int e=0; e<nEvents; e++) {
			if(otherCuts.size() == nEvents && !otherCuts[e])
				continue;
			if(!Trig->isDataEvent(e) || !Trig->beta2of4(e,mySide) || Trig->isBeamout(e))
				continue;
			bool tfired = Trig->betaTrigger(e,mySide,t);
			if(Trig->nFiring(e,mySide)-tfired > 1) {
				hAll.Fill(tubedat[t][e]);
				if(tfired)
					hTrig.Fill(tubedat[t][e]);		
			}
		}
		
		// input spectra plots
		if(allPlots) {
			defaultCanvas->SetLogy(true);
			hAll.SetLineColor(4);
			hAll.GetXaxis()->SetTitle("ADC channels above pedestal");
			hAll.GetYaxis()->SetTitle("N Events");
			hAll.Draw();
			hTrig.SetLineColor(2);
			hTrig.Draw("Same");
			printCanvas(std::string("Trigger_Effic_Input_")+itos(t));
			defaultCanvas->SetLogy(false);
		}
		
		TGraphAsymmErrors gEffic(hTrig.GetNbinsX());
		gEffic.BayesDivide(&hTrig,&hAll,"w");
		
		//hTrig.Divide(&hTrig,&hAll,1,1,"B");
		
		// scan for 50% point
		int b = gEffic.GetN();
		double midx,y;
		while(b > 1) {
			gEffic.GetPoint(--b,midx,y);
			if(y < 0.5)
				break;
		}	
		
		/*TF1 efficfit("efficfit",&poiscdf,lowx-50,lowx+250,4);
		 efficfit.SetParameter(0,midx);
		 efficfit.SetParameter(1,50);
		 efficfit.SetParameter(2,2.0);
		 efficfit.SetParameter(3,1.0);
		 efficfit.SetParLimits(1,5.0,200.0);
		 efficfit.SetParLimits(2,0.5,10.0); */
		
		TF1 efficfit("efficfit",&fancyfish,-50,200,4);
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
		gEffic.Fit(&efficfit,"Q");
		
		float_err trigef = float_err(efficfit.GetParameter(3),efficfit.GetParError(3));
		float_err trigc = float_err(efficfit.GetParameter(0),efficfit.GetParError(0));
		float_err trigw = float_err(efficfit.GetParameter(1),efficfit.GetParError(1));
		float_err trign_adj = float_err(efficfit.GetParameter(2),efficfit.GetParError(2));
		float trign = trigc.x/trigw.x*trign_adj.x;
		
		printf("Poisson CDF Fit: h = %.4f(%.4f), x0 = %.1f(%.1f), dx = %.1f(%.1f), n = %.2f [adjust %.2f(%.2f)]\n",
			   trigef.x, trigef.err, trigc.x, trigc.err, trigw.x, trigw.err, trign, trign_adj.x, trign_adj.err);
		
		forParent.insert(std::string("trig_effic_")+itos(t),vtos(efficfit.GetParameters(),efficfit.GetParameters()+4));
		forParent.insert(std::string("trig_effic_err_")+itos(t),vtos(efficfit.GetParErrors(),efficfit.GetParErrors()+4));
		
		if(trigef.x < 0.98) {
			Stringmap m;
			sprintf(tmp,"%c",sideNames(mySide));
			m.insert("side",tmp);
			m.insert("tube",t);
			m.insert("efficiency",trigef.x);
			m.insert("effic_delta",trigef.err);
			warn(MODERATE_WARNING,"Low_Trigger_Efficiency",m);
		}
		
		gEffic.SetMinimum(-0.10);
		gEffic.SetMaximum(1.10);
		gEffic.Draw("AP");
		sprintf(tmp,"%c%i PMT Trigger Efficiency",sideNames(mySide),t+1);
		gEffic.SetTitle(tmp);
		gEffic.GetXaxis()->SetTitle("ADC channels above pedestal");
		gEffic.GetXaxis()->SetLimits(-50,200);
		gEffic.GetYaxis()->SetTitle("Efficiency");
		gEffic.Draw("AP");
		
		printCanvas(std::string("Trigger_Efficiency_")+itos(t));
	}
}

