#include "ucnaDataAnalyzer11b.hh"
#include "GraphicsUtils.hh"
#include "SourceDBSQL.hh"
#include <TSpectrum2.h>
void ucnaDataAnalyzer11b::setupHistograms() {
	
	unsigned int nTimeBins = wallTime/10.;
	nTimeBins = nTimeBins?nTimeBins:1;
	float tpad = 10.0;
	
	for(Side s = EAST; s <= WEST; ++s) {
		hCathMax[s][0] = registeredTH1F(sideSubst("hCathMax_%c",s),sideSubst("%s Max Cathode",s),200,-200,4200);
		hCathMax[s][1] = registeredTH1F(sideSubst("hCathMaxCut_%c",s),sideSubst("%s Max Cathode, MWPC Cut",s),200,-200,4200);
		hCathMax[s][1]->SetLineColor(2);
		hCathMax[s][0]->GetXaxis()->SetTitle("ADC Channels");
		hCathMaxSum[s][0] = registeredTH1F(sideSubst("hCathMaxSum_%c",s),sideSubst("%s Max Cathode Sum",s),200,-200,8200);
		hCathMaxSum[s][1] = registeredTH1F(sideSubst("hCathMaxSumCut_%c",s),sideSubst("%s Max Cathode Sum, MWPC Cut",s),200,-200,8200);
		hCathMaxSum[s][1]->SetLineColor(2);
		hCathMaxSum[s][0]->GetXaxis()->SetTitle("ADC Channels");
		hAnode[s][0] = registeredTH1F(sideSubst("hAnode_%c",s),sideSubst("%s Anode",s),200,-200,4000);
		hAnode[s][1] = registeredTH1F(sideSubst("hAnodeCut_%c",s),sideSubst("%s Anode, MWPC Cut",s),200,-200,4000);
		hAnode[s][1]->SetLineColor(2);
		hAnode[s][0]->GetXaxis()->SetTitle("ADC Channels");
		hCathSum[s][0] = registeredTH1F(sideSubst("hCathSum_%c",s),sideSubst("%s Cathode Sum",s),200,-1000,40000);
		hCathSum[s][1] = registeredTH1F(sideSubst("hCathSumCut_%c",s),sideSubst("%s Cathode Sum, MWPC Cut",s),200,-1000,40000);
		hCathSum[s][1]->SetLineColor(2);
		hCathSum[s][0]->GetXaxis()->SetTitle("ADC Channels");
		
		hBackTDC[s] = registeredTH1F(sideSubst("hBackTDC_%c",s),"Backing Veto TDC",200,0,4000);
		hBackTDC[s]->SetLineColor(2+2*s);
		hBackTDC[s]->GetXaxis()->SetTitle("TDC Channels");
		hBackADC[s][0] = registeredTH1F(sideSubst("hBackADC_%c",s),sideSubst("%s Backing Veto ADC",s),200,-200,4000);
		hBackADC[s][1] = registeredTH1F(sideSubst("hBackADCCut_%c",s),sideSubst("%s Backing Veto ADC, TDC Cut",s),200,-200,4000);
		hBackADC[s][1]->SetLineColor(2);
		hDriftTAC[s] = registeredTH1F(sideSubst("hDriftTAC_%c",s),"Drift Tubes TAC",200,-200,4000);
		hDriftTAC[s]->SetLineColor(2+2*s);
		if(s==EAST) {
			hTopTDC[s] = registeredTH1F(sideSubst("hTopTDC_%c",s),"Top Veto TDC",200,0,4000);
			hTopTDC[s]->SetLineColor(2+2*s);
			hBackTDC[s]->GetXaxis()->SetTitle("TDC Channels");
			hTopADC[s][0] = registeredTH1F(sideSubst("hTopADC_%c",s),sideSubst("%s Top Veto ADC",s),200,-200,4000);
			hTopADC[s][1] = registeredTH1F(sideSubst("hTopADCCut_%c",s),sideSubst("%s Top Veto ADC, with TDC Cut",s),200,-200,4000);
			hTopADC[s][1]->SetLineColor(2);
			hTopADC[s][0]->GetXaxis()->SetTitle("ADC Channels");
		}
		for(unsigned int t=0; t<=nBetaTubes; t++) {
			hScintTDC[s][t] = registeredTH1F(sideSubst("hScintTDC_%c",s)+(t<nBetaTubes?itos(t):""),"2-of-4 TDC",200,0,5000);
			hScintTDC[s][t]->SetLineColor(t<nBetaTubes?t+2:2+2*s);
			hScintTDC[s][t]->GetXaxis()->SetTitle("TDC Channels");
		}
		
		for(unsigned int t=TYPE_0_EVENT; t<=TYPE_IV_EVENT; t++) {
			hEtrue[s][t] = registeredTH1F(sideSubst("hErecon_%c_",s)+itos(t),"Reconstructed Energy, Type "+itos(t),200,0,1500);
			hEtrue[s][t]->SetLineColor(2+2*s);
			hEtrue[s][t]->GetXaxis()->SetTitle("Erecon [keV]");
		}
		
		for(unsigned int t=0; t<nBetaTubes; t++) {
			hTuben[s][t] = registeredTH1F(sideSubst("hTuben_%c_",s)+itos(t),"Visible Energy",200,-100,1500);
			hTuben[s][t]->SetLineColor(2+t);
			hTuben[s][t]->GetXaxis()->SetTitle("Erecon [keV]");
		}
		
		hSideRate[s][1] = registeredTH1F(sideSubst("hBetaRate_%c",s),"Beta Event Rate",nTimeBins,-tpad,wallTime+tpad);
		hSideRate[s][1]->SetLineColor(2+2*s);
		hSideRate[s][1]->GetXaxis()->SetTitle("Time [s]");
		
		hSideRate[s][0] = registeredTH1F(sideSubst("hMuonRate_%c",s),"Muon Event Rate",nTimeBins,-tpad,wallTime+tpad);
		hSideRate[s][0]->SetLineColor(1+2*s);
		hSideRate[s][0]->GetXaxis()->SetTitle("Time [s]");
		
		hHitPos[s] = registeredTH2F(sideSubst("HitPos_%c",s),sideSubst("%s Hit Positions",s),400,-65,65,400,-65,65);
		hHitPos[s]->GetXaxis()->SetTitle("x position [mm]");
		hHitPos[s]->GetYaxis()->SetTitle("y position [mm]");
		for(AxisDirection d = X_DIRECTION; d <= Y_DIRECTION; ++d) {
			hHitsProfile[s][d] = registeredTH1F(sideSubst("HitPos_%c",s)+(d==X_DIRECTION?"x":"y"),
												std::string(d==X_DIRECTION?"x":"y")+" Hit Positions",200,-65,65);
			hHitsProfile[s][d]->SetLineColor(2+2*s);
			hHitsProfile[s][d]->GetXaxis()->SetTitle((std::string(d==X_DIRECTION?"x":"y")+" position [mm]").c_str());
		}
		
		// trigger efficiency, pulser spectrum
		unsigned int nBiDivs = int(wallTime/300)?wallTime/300:1;
		for(unsigned int t=0; t<nBetaTubes; t++)
			sevt[s].adc[t]=3950;
		PCal.pedSubtract(s,sevt[s].adc,0);
		for(unsigned int t=0; t<nBetaTubes; t++) {
			hTrigEffic[s][t][0] = registeredTH1F(sideSubst("hTrigEfficAll_%c",s)+itos(t),"Trigger Efficiency Events",100,-50,150);
			hTrigEffic[s][t][0]->SetLineColor(4);
			hTrigEffic[s][t][1] = registeredTH1F(sideSubst("hTrigEfficTrig_%c",s)+itos(t),"Trigger Efficiency Events",100,-50,150);
			hTrigEffic[s][t][1]->SetLineColor(2);
			for(unsigned int i=0; i<nBiDivs; i++) {
				hBiPulser[s][t].push_back(registeredTH1F(sideSubst("hPulserSpectrum_%c",s)+itos(t)+"_"+itos(i),"Bi Pulser Spectrum",
														 200,100,sevt[s].adc[t]));
				hBiPulser[s][t].back()->SetLineColor(t+1);
			}
		}
	}
	
	for(unsigned int t=TYPE_0_EVENT; t<=TYPE_IV_EVENT; t++) {
		hTypeRate[t] = registeredTH1F("hTypeRate_"+itos(t),"Type "+itos(t)+" Event Rate",nTimeBins,-tpad,wallTime+tpad);
		hTypeRate[t]->SetLineColor(5+t);
		hTypeRate[t]->GetXaxis()->SetTitle("Time [s]");
	}
	
	for(unsigned int n=0; n<kNumUCNMons; n++) {
		hMonADC[n] = registeredTH1F("UCN_Mon_"+itos(n+1)+"_ADC","UCN Mon "+itos(n+1)+" ADC",200,0,4000);
		hMonADC[n]->SetLineColor(n+1);
		hMonADC[n]->GetXaxis()->SetTitle("ADC Channels");
		hMonRate[n] = registeredTH1F("UCN_Mon_"+itos(n+1)+"_Rate","UCN Mon "+itos(n+1)+" Rate",nTimeBins,-tpad,wallTime+tpad);
		hMonRate[n]->SetLineColor(n+1);
		hMonRate[n]->GetXaxis()->SetTitle("Time [s]");
	}
	
	hEvnbFailRate = registeredTH1F("EvnbFail","Evnb Fail Rate",nTimeBins,-tpad,wallTime+tpad);
	hEvnbFailRate->SetLineColor(kOrange+10);
	hBkhfFailRate = registeredTH1F("BkhfFail","Bkhf Fail Rate",nTimeBins,-tpad,wallTime+tpad);
	hBkhfFailRate->SetLineColor(2);
	
	hClusterTiming[0] = registeredTH1F("ClusterTimingAll","Events in time window",100,1.0,5.0);
	hClusterTiming[0]->GetXaxis()->SetTitle("log microseconds");
	hClusterTiming[1] = registeredTH1F("ClusterTimingTrig","Events in time window",100,1.0,5.0);
	hClusterTiming[1]->SetLineColor(4);
	hClusterTiming[1]->GetXaxis()->SetTitle("log microseconds");
}

void ucnaDataAnalyzer11b::fillEarlyHistograms() {
	for(unsigned int n=0; n<kNumUCNMons; n++) {
		if(isUCNMon(n)) {
			hMonADC[n]->Fill(fMonADC[n].val);
			hMonRate[n]->Fill(fTimeScaler[BOTH]);
		}
	}
	if(isPulserTrigger()) {
		unsigned int tbin = fTimeScaler[BOTH]*hBiPulser[EAST][0].size()/wallTime;
		if(tbin<hBiPulser[EAST][0].size())
			for(Side s = EAST; s <= WEST; ++s)
				for(unsigned int t=0; t<nBetaTubes; t++)
					hBiPulser[s][t][tbin]->Fill(sevt[s].adc[t]);
	}
	if(!fEvnbGood)
		hEvnbFailRate->Fill(fTimeScaler[BOTH]);
	if(!fBkhfGood)
		hBkhfFailRate->Fill(fTimeScaler[BOTH]);
}

void ucnaDataAnalyzer11b::fillHistograms() {
	
	// already not LED, non-Scint triggers
	
	// gamma events
	if(fType == TYPE_IV_EVENT && fSide <= WEST && fPID != PID_PULSER) {
		hTypeRate[fType]->Fill(fTimeScaler[BOTH]);
		hEtrue[fSide][fType]->Fill(fEtrue);
	}
	
	for(Side s = EAST; s <= WEST; ++s) {
		
		// wirechambers
		hCathMax[s][0]->Fill(fCathMax[s].val);
		hCathMaxSum[s][0]->Fill(fCathMaxSum[s].val);
		hCathSum[s][0]->Fill(fCathSum[s].val);
		hAnode[s][0]->Fill(fMWPC_anode[s].val);
		if(passedMWPC(s)) {
			hCathMax[s][1]->Fill(fCathMax[s].val);
			hCathMaxSum[s][1]->Fill(fCathMaxSum[s].val);
			hCathSum[s][1]->Fill(fCathSum[s].val);
			hAnode[s][1]->Fill(fMWPC_anode[s].val);
			if(fPID==PID_BETA && fSide==s) {
				hHitPos[s]->Fill(wirePos[s][X_DIRECTION].center,wirePos[s][Y_DIRECTION].center);
				for(AxisDirection d = X_DIRECTION; d <= Y_DIRECTION; ++d)
					hHitsProfile[s][d]->Fill(wirePos[s][d].center);
			}
		}
		
		
		// muon vetos
		hBackTDC[s]->Fill(fBacking_tdc[s].val);
		hBackADC[s][0]->Fill(fBacking_adc[s]);
		if(fTaggedBack[s])
			hBackADC[s][1]->Fill(fBacking_adc[s]);
		hDriftTAC[s]->Fill(fDrift_tac[s].val);
		if(s==EAST) {
			hTopTDC[s]->Fill(fTop_tdc[s].val);
			hTopADC[s][0]->Fill(fTop_adc[s]);
			if(fTaggedTop[s])
				hTopADC[s][1]->Fill(fTop_adc[s]);
		}
		
		// trigger efficiency
		for(unsigned int t=0; t<nBetaTubes; t++) {
			int nf = nFiring(s);
			bool tfired = pmtFired(s,t);
			if(nf-tfired<2 || iSis00!=1+s || !fPassedGlobal)
				continue;
			hTrigEffic[s][t][0]->Fill(sevt[s].adc[t]);
			if(tfired)
				hTrigEffic[s][t][1]->Fill(sevt[s].adc[t]);	
		}
		
		if(fPID==PID_BETA)
			for(unsigned int t=0; t<=nBetaTubes; t++)
				hScintTDC[s][t]->Fill(fScint_tdc[s][t].val);
		
		if(fSide != s)
			continue;
		
		if(fPID==PID_BETA) {
			hSideRate[s][1]->Fill(fTimeScaler[BOTH]);
			if(fType <= TYPE_III_EVENT && fPassedGlobal)
				hEtrue[s][fType]->Fill(fEtrue);
			if(fType == TYPE_0_EVENT && fPassedGlobal)
				for(unsigned int t=0; t<nBetaTubes; t++)
					hTuben[s][t]->Fill(sevt[s].tuben[t].x);
		} else if(fPID==PID_MUON) {
			hSideRate[s][0]->Fill(fTimeScaler[BOTH]);
		}
	}
	
	if(fPID==PID_BETA)
		hClusterTiming[0]->Fill(log10(fWindow.val*1.e6));
	if(fPID==PID_BETA && fType<=TYPE_III_EVENT) {
		hTypeRate[fType]->Fill(fTimeScaler[BOTH]);
		if(fPassedGlobal)
			hClusterTiming[1]->Fill(log10(fWindow.val*1.e6));
	}
}

void ucnaDataAnalyzer11b::drawCutRange(const RangeCut& r, Int_t c) {
	drawVLine(r.start, defaultCanvas, c);
	drawVLine(r.end, defaultCanvas, c);
}

void ucnaDataAnalyzer11b::drawExclusionBlips(Int_t c) {
	for(std::vector<Blip>::const_iterator it = cutBlips.begin(); it != cutBlips.end(); it++)
		if(it->end[BOTH]-it->start[BOTH] > 5)
			drawExcludedRegion(it->start[BOTH], it->end[BOTH], defaultCanvas,c,3354);
		else
			drawExcludedRegion(it->start[BOTH], it->end[BOTH], defaultCanvas,c,1001);
}

void ucnaDataAnalyzer11b::plotHistos() {
	printf("\nMaking output plots...\n");
	defaultCanvas->cd();
	defaultCanvas->SetLogy(true);
	std::vector<TH1*> hToPlot;
	
	for(Side s = EAST; s <= WEST; ++s) {
		// wirechambers
		hCathMax[s][0]->Draw();
		hCathMax[s][1]->Draw("Same");
		drawCutRange(fCathMax[s].R);
		printCanvas(sideSubst("Wirechamber/CathMax_%c",s));
		
		hCathMaxSum[s][0]->Draw();
		hCathMaxSum[s][1]->Draw("Same");
		drawCutRange(fCathMaxSum[s].R);
		printCanvas(sideSubst("Wirechamber/CathMaxSum_%c",s));
		
		hAnode[s][0]->Draw();
		hAnode[s][1]->Draw("Same");
		drawCutRange(fMWPC_anode[s].R);
		printCanvas(sideSubst("Wirechamber/Anode_%c",s));
		
		hCathSum[s][0]->Draw();
		hCathSum[s][1]->Draw("Same");
		drawCutRange(fCathSum[s].R);
		printCanvas(sideSubst("Wirechamber/CathSum_%c",s));
		
		defaultCanvas->SetLogy(false);
		hHitPos[s]->Draw("Col");
		printCanvas(sideSubst("Wirechamber/HitPos_%c",s));
		defaultCanvas->SetLogy(true);
		
		// muon vetos
		hBackADC[s][0]->Draw();
		hBackADC[s][1]->Draw("Same");
		printCanvas(sideSubst("MuVeto/BackADC_%c",s));
		if(s==EAST) {
			hTopADC[s][0]->Draw();
			hTopADC[s][1]->Draw("Same");
			printCanvas(sideSubst("MuVeto/TopADC_%c",s));
		}
		
		// trigger efficiency
		for(unsigned int t=0; t<nBetaTubes; t++) {
			hTrigEffic[s][t][0]->Draw();
			hTrigEffic[s][t][1]->Draw("Same");
			printCanvas(sideSubst("PMTs/TrigEffic_Input_%c",s)+itos(t));
		}
		
		// PMT TDCs
		hToPlot.clear();
		for(unsigned int t=0; t<nBetaTubes; t++)
			hToPlot.push_back(hScintTDC[s][t]);
		drawSimulHistos(hToPlot);
		for(unsigned int t=0; t<nBetaTubes; t++)
			drawCutRange(fScint_tdc[s][t].R,2+t);
		printCanvas(sideSubst("PMTs/TDC_2of4_%c",s));
	}
	
	// 1-D hit positions
	defaultCanvas->SetLogy(false);
	for(AxisDirection d = X_DIRECTION; d <= Y_DIRECTION; ++d) {
		hToPlot.clear();
		for(Side s = EAST; s <= WEST; ++s)
			hToPlot.push_back(hHitsProfile[s][d]);
		drawSimulHistos(hToPlot);
		printCanvas("Wirechamber/HitPos_"+std::string(d==X_DIRECTION?"x":"y"));
	}
	defaultCanvas->SetLogy(true);
	
	// muon vetos
	hBackTDC[EAST]->Draw();
	hBackTDC[WEST]->Draw("Same");
	drawCutRange(fBacking_tdc[EAST].R,2);
	drawCutRange(fBacking_tdc[WEST].R,4);
	printCanvas("MuVeto/Backing_TDC");
	
	hTopTDC[EAST]->Draw();
	drawCutRange(fTop_tdc[EAST].R,2);
	printCanvas("MuVeto/Top_TDC");
	
	hDriftTAC[EAST]->Draw();
	hDriftTAC[WEST]->Draw("Same");
	drawCutRange(fDrift_tac[EAST].R,2);
	drawCutRange(fDrift_tac[WEST].R,4);
	printCanvas("MuVeto/Drift_TAC");
	
	// PMT TDCs
	hScintTDC[EAST][nBetaTubes]->Draw();
	hScintTDC[WEST][nBetaTubes]->Draw("Same");
	drawCutRange(fScint_tdc[EAST][nBetaTubes].R,2);
	drawCutRange(fScint_tdc[WEST][nBetaTubes].R,4);
	drawCutRange(ScintSelftrig[EAST],2);
	drawCutRange(ScintSelftrig[WEST],4);
	printCanvas("PMTs/TDC_2of4");
	
	// GV monitors
	hToPlot.clear();
	for(unsigned int n=0; n<kNumUCNMons; n++) {
		hMonRate[n]->Scale(1.0/hMonRate[n]->GetBinWidth(1));
		hToPlot.push_back(hMonRate[n]);
		hMonADC[n]->Draw();
		drawCutRange(fMonADC[n].R,4);
		printCanvas("UCN_Mon/Mon_"+itos(n)+"_ADC");
	}
	drawSimulHistos(hToPlot);
	printCanvas("UCN_Mon/Mon_Rates");
	
	// rates
	hToPlot.clear();
	for(unsigned int t=0; t<=TYPE_II_EVENT; t++)
		hToPlot.push_back(hTypeRate[t]);
	hToPlot.push_back(hTypeRate[TYPE_IV_EVENT]);
	for(Side s = EAST; s <= WEST; ++s) {
		hToPlot.push_back(hSideRate[s][0]);
		hToPlot.push_back(hSideRate[s][1]);
	}
	for(std::vector<TH1*>::iterator it = hToPlot.begin(); it  != hToPlot.end(); it++) {
		(*it)->Scale(1.0/(*it)->GetBinWidth(1));
		(*it)->SetMinimum(0.02);
		(*it)->SetMaximum(200.0);
	}
	drawSimulHistos(hToPlot);
	drawExclusionBlips(4);
	printCanvas("Rates/EventRates");
	
	hToPlot.clear();
	hToPlot.push_back(hEvnbFailRate);
	hEvnbFailRate->Scale(1.0/hEvnbFailRate->GetBinWidth(1.0));
	hToPlot.push_back(hBkhfFailRate);
	hBkhfFailRate->Scale(1.0/hBkhfFailRate->GetBinWidth(1.0));
	hToPlot.push_back(hMonRate[UCN_MON_GV]);
	hToPlot.push_back(hTypeRate[TYPE_0_EVENT]);
	for(std::vector<TH1*>::iterator it = hToPlot.begin(); it  != hToPlot.end(); it++) {
		(*it)->SetMinimum(0.01);
		(*it)->SetMaximum(100.0);
	}	
	drawSimulHistos(hToPlot);
	drawExclusionBlips(4);
	printCanvas("Rates/GlobalRates");
	
	hClusterTiming[0]->Draw();
	hClusterTiming[1]->Draw("Same");
	drawCutRange(RangeCut(log10(fWindow.R.start*1.e6),log10(fWindow.R.end*1.e6)),2);
	printCanvas("Rates/ClusterTiming");
	
	// energy by event type
	defaultCanvas->SetLogy(false);
	for(unsigned int t=TYPE_0_EVENT; t<=TYPE_IV_EVENT; t++) {
		hToPlot.clear();
		for(Side s = EAST; s <= WEST; ++s)
			hToPlot.push_back(hEtrue[s][t]);
		drawSimulHistos(hToPlot);
		printCanvas("PMTs/Erecon_Type_"+itos(t));
	}
	// energy by PMT
	for(Side s = EAST; s <= WEST; ++s) {
		hToPlot.clear();
		for(unsigned int t=0; t<nBetaTubes; t++)
			hToPlot.push_back(hTuben[s][t]);
		drawSimulHistos(hToPlot);
		printCanvas(sideSubst("PMTs/Evis_PMTs_%c",s));
	}
}

void ucnaDataAnalyzer11b::muonVetoAccidentals() {
	for(Side s = EAST; s <= WEST; ++s) {
		int b0 = hDriftTAC[s]->FindBin(300);
		int b1 = hDriftTAC[s]->FindBin(fDrift_tac[s].R.start);
		int b2 = hDriftTAC[s]->FindBin(fDrift_tac[s].R.end);
		int b3 = hDriftTAC[s]->FindBin(s==EAST?4000:2900);
		double nBg = hDriftTAC[s]->Integral(b0,b1)+hDriftTAC[s]->Integral(b2,b3);
		double tBg = (hDriftTAC[s]->GetBinLowEdge(b1+1) - hDriftTAC[s]->GetBinLowEdge(b0)
					  + hDriftTAC[s]->GetBinLowEdge(b3+1) - hDriftTAC[s]->GetBinLowEdge(b2));
		double tCut = hDriftTAC[s]->GetBinLowEdge(b2) - hDriftTAC[s]->GetBinLowEdge(b1+1);
		double nAcc = nBg*tCut/tBg;
		double nVetod = hDriftTAC[s]->Integral(b1+1,b2-1);
		double bgTot = hDriftTAC[s]->Integral()+hDriftTAC[s]->GetBinContent(hDriftTAC[s]->GetNbinsX()+1)-nVetod+nAcc;
		Stringmap m;
		m.insert("side",sideSubst("%c",s));
		m.insert("bgRate",nBg/tBg);
		m.insert("tCut",tCut);
		m.insert("nAcc",nAcc);
		m.insert("nVeto",nVetod-nAcc);
		m.insert("nBG",bgTot);
		m.insert("accFrac",nAcc/bgTot);
		qOut.insert("driftTAC_accidentals",m);
	}
}

// comparison to sort sources by x position in descending order
bool sortSourceXPos(Source a, Source b) { return a.x < b.x; }

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

void ucnaDataAnalyzer11b::locateSourcePositions() {
	
	if(!CDBout) {
		printf("Skipping source location without DB write permission.\n");
		return;
	}
		
	for(Side s = EAST; s <= WEST; ++s) {
		
		unsigned int nbins = hHitPos[s]->GetNbinsX()-2;
		
		std::vector<Source> expectedSources = SourceDBSQL::getSourceDBSQL()->runSources(rn,s);
		if(!expectedSources.size())
			continue;
		
		// convert hit position histogram to float** array
		float** zout = allocArray(nbins,nbins);
		float** historray = allocArray(nbins,nbins);
		for(unsigned int x=0; x<nbins; x++)
			for(unsigned int y=0; y<nbins; y++)
				historray[x][y] = hHitPos[s]->GetBinContent(x+1,y+1);
		
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
		
		printf("\t%i of %i sources found.\n",npks,(int)expectedSources.size());
		
		// extract located peaks into source positions vector
		// drop peaks past expected number of sources
		Float_t* xPos = TS.GetPositionX();
		Float_t* yPos = TS.GetPositionY();
		std::vector<Source> foundSources;
		for(unsigned int i=0; i<npks; i++) {
			if(i>=expectedSources.size())
				break;
			Source pk;
			pk.x = binterpolate(hHitPos[s]->GetXaxis(),xPos[i]);
			pk.wx = 2.0;
			pk.y = binterpolate(hHitPos[s]->GetYaxis(),yPos[i]);
			pk.wy = 2.0;
			pk.mySide = s;
			pk.myRun = rn;
			foundSources.push_back(pk);
		}
		
		
		// sort sources and apply source types from expected sources
		std::sort(foundSources.begin(),foundSources.end(),&sortSourceXPos);
		std::sort(expectedSources.begin(),expectedSources.end(),&sortSourceXPos);
		for(unsigned int i=0; i < foundSources.size(); i++) {
			foundSources[i].t = expectedSources[i].t;
			foundSources[i].sID = expectedSources[i].sID;
		}
		
		// refine source locations and widths; upload to DB
		for(std::vector<Source>::iterator it = foundSources.begin(); it!=foundSources.end(); it++) {
			TH1D* px = hHitPos[s]->ProjectionX("_px",hHitPos[s]->GetYaxis()->FindBin(it->y-3*it->wy),hHitPos[s]->GetYaxis()->FindBin(it->y+3*it->wy));
			TH1D* py = hHitPos[s]->ProjectionY("_py",hHitPos[s]->GetXaxis()->FindBin(it->x-3*it->wx),hHitPos[s]->GetXaxis()->FindBin(it->x+3*it->wx));
			TF1 g1("g1","gaus"); 
			g1.SetLineColor(2);
			px->Fit(&g1,"Q","",it->x-3*it->wx,it->x+3*it->wx);
			it->x = g1.GetParameter(1);
			it->wx = fabs(g1.GetParameter(2));
			it->nCounts = px->Integral(px->FindBin(it->x-3*it->wx),px->FindBin(it->x+3*it->wx));
			py->Fit(&g1,"Q","",it->y-3*it->wy,it->y+3*it->wy);
			it->y = g1.GetParameter(1);
			it->wy = g1.GetParameter(2);
			
			px->Draw();
			printCanvas("Wirechamber/Source_"+it->name()+"_x");
			py->Draw();
			printCanvas("Wirechamber/Source_"+it->name()+"_y");
			delete(px);
			delete(py);
			
			it->display();
			SourceDBSQL::getSourceDBSQL()->addSource(*it);
		}
		
		deleteArray(historray,nbins);
		deleteArray(zout,nbins);
	}
}

