#include "LEDTracker.hh"
#include "GraphicsUtils.hh"

void LEDTracker::specialize() {
	
	loadData("Qadc2","ADCE3Beta","Etubes");	// 0
	loadData("Qadc3","ADCE4Beta");			// 1
	loadData("Qadc0","ADCE1Beta");			// 2
	loadData("Qadc1","ADCE2Beta");			// 3
	loadData("Pdc32","ADCRefELED","reftube");		// 4 GMS reference phototube
	
	loadData("Qadc4","ADCW1Beta","Wtubes");	// 0
	loadData("Qadc5","ADCW2Beta");			// 1
	loadData("Qadc6","ADCW3Beta");			// 2
	loadData("Qadc7","ADCW4Beta");			// 3
	loadData("Pdc36","ADCRefWLED","reftube");	// 4 GMS reference phototube		

	for(unsigned int t=0; t<=nBetaTubes; t++) {
		tubedat[EAST].push_back(getData("Etubes",t));
		tubedat[WEST].push_back(getData("Wtubes",t));
	}
}

// function for determining beta scintillator pedestals
std::vector< std::pair<float,float> > LED_pedestal_finder(Subsystem* S, void* tdat) {
	LEDTracker* LT = (LEDTracker*)S;
	float* tubedat = (float*)tdat;
	std::vector< std::pair<float,float> > v;
	for(unsigned int e=0; e<LT->nEvents; e++)
		if(!LT->Trig->gmsCo(e,LT->s) && !LT->Trig->beta2of4(e,LT->s) && !LT->Trig->gmsLED(e) && !LT->Trig->isCrud(e))
			v.push_back(std::make_pair(LT->Trig->eventTime(e),tubedat[e]));
	return v;
}

LEDTracker::LEDTracker(RunManager* T,Trigger* TG,unsigned int scanSteps, float stepSize, float stepTime): Subsystem(T,"LEDTracker",BOTH), Trig(TG) {

	specialize();
	char tmp[1024];
	
	float scanTime = stepTime*scanSteps;	// nominal time of scan, seconds
	float scanRange = (scanSteps-1)*stepSize;	// nominal range of scan, DB
	float transTime = 3.0;					// transition time to cut out
	float tdelta = 2.0;						// offset from transition start
	
	// pedestal subtraction
	printf("LEDTracker pedestal subtraction...\n");
	for(s = EAST; s <= WEST; s = nextSide(s)) {
		for(unsigned int t=0; t<=nBetaTubes; t++) {
			if(true || !verifyPedestal(t+5*s))
				monitorPeak(LED_pedestal_finder, 60.0, 3000, (void*)tubedat[s][t], sensorNames[t+5*s]);
			for(unsigned int e=0; e<nEvents; e++)
				tubedat[s][t][e] -= PC.getPedestal(sensorNames[t+5*s],TG->eventTime(e));
		}
	}
	printf("\tPedestal subtraction complete.\n");	
	
	// LED peak tracking TProfiles
	printf("LEDTracker LED peak tracking...\n");
	float binWidth = 2.0;	// data binning width, seconds
	unsigned int ndivs = (unsigned int)(TG->totalTime()/binWidth);	// number of bins
	if(!ndivs)
		ndivs = 1;
	printf("Dividing into %i segments.\n",ndivs);
	TProfile* tracks[2][nBetaTubes+1];
	TH1F* hists[2][nBetaTubes+1];
	for(s = EAST; s <= WEST; s = nextSide(s)) {
		for(unsigned int t=0; t<=nBetaTubes; t++) {
			sprintf(tmp,"LED_%c%i",sideNames(s),t);
			tracks[s][t] = new TProfile(tmp,"LED Peak",ndivs,0,Trig->totalTime());
			sprintf(tmp,"hLED_%c%i",sideNames(s),t);
			hists[s][t] = new TH1F(tmp,"LED Peak",200,-100,4000);
		}
	}
	printf("Filling histograms...\n");
	for(unsigned int e=0; e<nEvents; e++)
		if(Trig->gmsLED(e) && !Trig->isCrud(e))
			for(s = EAST; s <= WEST; s = nextSide(s))
				for(unsigned int t=0; t<=nBetaTubes; t++)
					tracks[s][t]->Fill(Trig->eventTime(e),tubedat[s][t][e]);
	printf("\tLED tracking complete.\n");	
	
	
	// identify tube min/max ranges
	printf("Identifying data range...\n");
	float ymin[2][nBetaTubes+1];
	float ymax[2][nBetaTubes+1];
	for(s = EAST; s <= WEST; s = nextSide(s)) {
		for(unsigned int t=0; t<=nBetaTubes; t++) {
			for(unsigned int n=0; n<ndivs; n++)
				hists[s][t]->Fill(tracks[s][t]->GetBinContent(n+1));
			ymin[s][t] = 10000.0;
			ymax[s][t] = -10000.0;
			for(int n=1; n<=hists[s][t]->GetNbinsX(); n++) {
				if(hists[s][t]->GetBinContent(n)) {
					if(hists[s][t]->GetBinCenter(n) < ymin[s][t])
						ymin[s][t] = hists[s][t]->GetBinCenter(n);
					if(hists[s][t]->GetBinCenter(n) > ymax[s][t])
						ymax[s][t] = hists[s][t]->GetBinCenter(n);
				}
			}
			printf("\t%c%i %.1f-%.1f\n",sideNames(s),t,ymin[s][t],ymax[s][t]);
		}
	}
	
	// identify scan start points at max->min transitions
	printf("Identifying scan ranges...\n");
	std::vector<float> scanStarts;
	int lastmax = -1;
	for(unsigned int n=0; n<ndivs; n++) {
		int nMax = 0;
		int nMin = 0;
		for(s = EAST; s <= WEST; s = nextSide(s)) {
			for(unsigned int t=0; t<=nBetaTubes; t++) {
				nMax += (tracks[s][t]->GetBinContent(n+1) > ymax[s][t]-0.40*(ymax[s][t]-ymin[s][t]));
				nMin += (tracks[s][t]->GetBinContent(n+1) < ymin[s][t]+0.20*(ymax[s][t]-ymin[s][t]));
			}
		}
		if(nMax > 5)
			lastmax = n;
		if(nMin > 5) {
			if(lastmax >= 0 && tracks[EAST][0]->GetBinCenter(lastmax+1)+scanTime-20.0 < TG->totalTime())
				scanStarts.push_back(tracks[EAST][0]->GetBinCenter(lastmax+1));
			lastmax = -1;
		}
	}
	printf("\tFound %i scan starts.\n",(int)scanStarts.size());
	
	// fit and record scan steps
	Stringmap mm;
	mm.insert("scanSteps",scanSteps);
	mm.insert("scanRange",scanRange);
	mm.insert("scanTime",scanTime);
	mm.insert("transTime",transTime);
	mm.insert("nScans",scanStarts.size());	
	TF1 pol0fit("pol0","pol0",0,10);
	pol0fit.SetLineColor(4);
	parent->qOut.insert("LED_Scan",mm);
	unsigned int e0 = 0;
	
	for(int i=0; i<int(scanStarts.size()); i++) {
		
		std::vector<TH1F*> hPeaks[2][nBetaTubes+1];
		
		for(unsigned int n=0; n<scanSteps; n++) {
			
			float t0 = scanStarts[i]+n*scanTime/scanSteps+transTime+tdelta;
			float t1 = scanStarts[i]+(n+1)*scanTime/scanSteps-transTime+tdelta;
			pol0fit.SetRange(t0,t1);
			
			// fit steps			
			for(s = EAST; s <= WEST; s = nextSide(s)) {
				for(unsigned int t=0; t<=nBetaTubes; t++) {
					
					tracks[s][t]->Fit(&pol0fit,"QR+");
					Stringmap m;
					m.insert("side",ctos(sideNames(s)));
					m.insert("tube",itos(t));
					m.insert("scan",itos(i));
					m.insert("point",itos(n));
					m.insert("y",pol0fit.GetParameter(0));
					m.insert("dy",pol0fit.GetParError(0));
					m.insert("x",pow(10,(n*scanRange/float(scanSteps-1.0)-scanRange)/10.0));
					m.insert("dB",n*scanRange/float(scanSteps-1.0)-scanRange);
					parent->qOut.insert("LED_Scan_Point",m);
					
					sprintf(tmp,"hPeaks_%c%i_%i",sideNames(s),t,n);
					hPeaks[s][t].push_back(new TH1F(tmp,"LED Peaks",1000,log(1.0),log(4000.0)));
					hPeaks[s][t].back()->SetLineColor(2+(n%6));
				}
			}
			
			// fill peaks histograms
			while(Trig->eventTime(e0) < t0 && e0 < nEvents)
				e0++;
			while(Trig->eventTime(e0) < t1 && e0 < nEvents) {
				for(s = EAST; s <= WEST; s = nextSide(s))
					for(unsigned int t=0; t<=nBetaTubes; t++)
						if(tubedat[s][t][e0] > 0)
							hPeaks[s][t][n]->Fill(log(tubedat[s][t][e0]));
				e0++;
			}
		}
		
		// peaks histograms plots
		for(s = EAST; s <= WEST; s = nextSide(s)) {
			for(unsigned int t=0; t<=nBetaTubes; t++) {
				
				defaultCanvas->SetLogy(true);
				hPeaks[s][t][0]->GetXaxis()->SetTitle("log ADC");
				hPeaks[s][t][0]->Draw();
				for(unsigned int n=1; n<scanSteps; n++)
					hPeaks[s][t][n]->Draw("Same");
				sprintf(tmp,"Peaks_%c%i_Scan_%i",sideNames(s),t,i);
				printCanvas(tmp);
				defaultCanvas->SetLogy(false);
				
				// cleanup
				for(unsigned int n=0; n<scanSteps; n++)
					delete(hPeaks[s][t][n]);
				
			}
		}

	}
	
	// make plots
	for(s = EAST; s <= WEST; s = nextSide(s)) {
		for(unsigned int t=0; t<=nBetaTubes; t++) {
			
			tracks[s][t]->Draw();
			for(unsigned int i=0; i<scanStarts.size(); i++)
				drawVLine(scanStarts[i], defaultCanvas, 2);
			printCanvas(std::string("LED_Peak_History_")+ctos(sideNames(s))+itos(t));
			hists[s][t]->Draw();
			drawVLine(ymin[s][t], defaultCanvas, 4);
			drawVLine(ymax[s][t], defaultCanvas, 2);
			printCanvas(std::string("LED_Peak_Dist_")+ctos(sideNames(s))+itos(t));
			
		}
	}
	
	// cleanup
	for(s = EAST; s <= WEST; s = nextSide(s)) {
		for(unsigned int t=0; t<=nBetaTubes; t++) {
			delete(tracks[s][t]);
			delete(hists[s][t]);
		}
	}
	
}
