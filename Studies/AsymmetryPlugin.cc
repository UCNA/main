#include "AsymmetryPlugin.hh"
#include "BetaSpectrum.hh"
#include "KurieFitter.hh"
#include "GraphicsUtils.hh"
#include "HighEnergyExcessPlugin.hh"
#include <TH2F.h>

/// nominal asymmetry for fit
Double_t asymmetryFitFunc(const Double_t* x, const Double_t* par) {
	return par[0]*plainAsymmetry(x[0],0.5)/A0_PDG;
}

TF1 AsymmetryPlugin::asymmetryFit = TF1("asymFit",&asymmetryFitFunc,0,neutronBetaEp,1);
TF1 AsymmetryPlugin::averagerFit = TF1("averagerFit","pol0",0,neutronBetaEp);

AsymmetryPlugin::AsymmetryPlugin(OctetAnalyzer* OA):
OctetAnalyzerPlugin(OA,"asymmetry"), nEnergyBins(150), energyMax(1500), hAsym(NULL), hCxn(NULL), anChoice(ANCHOICE_C) {
	for(Side s = EAST; s <= WEST; ++s) {
		myA->ignoreMissingHistos = true;
		qPassesWC[s] = registerCoreHist("hPassedWC","Events Passing Wirechamber Energy",nEnergyBins, 0, energyMax, s);
		myA->ignoreMissingHistos = false;
		for(EventType t=TYPE_0_EVENT; t<=TYPE_IV_EVENT; ++t) {
			if(t==TYPE_II_EVENT || t==TYPE_III_EVENT) {
				q23ProbCut[s][t] = registerCoreHist("h23ProbCut_Tp_"+itos(t),
													"Type "+itos(t)+" Events Energy",
													nEnergyBins, 0, energyMax, s);
				q23ProbCut[s][t]->setAxisTitle(X_DIRECTION,"Energy [keV]");
			}
			for(unsigned int p=0; p<=nBetaTubes; p++) {
				qEnergySpectra[s][p][t] = registerCoreHist("hEnergy_"+(p<nBetaTubes?itos(p)+"_":"")+"Type_"+itos(t),
														   "Type "+itos(t)+" Events Energy",
														   nEnergyBins, 0, energyMax, s);
				qEnergySpectra[s][p][t]->setAxisTitle(X_DIRECTION,"Energy [keV]");
				qEnergySpectra[s][p][t]->setRangeUser(0,800);
			}
		}
	}
}

void AsymmetryPlugin::fillCoreHists(ProcessedDataScanner& PDS, double weight) {
	Side s = PDS.fSide;
	if(!(s==EAST || s==WEST)) return;
	if(PDS.fPID == PID_SINGLE && PDS.fType == TYPE_IV_EVENT)
		qEnergySpectra[s][nBetaTubes][TYPE_IV_EVENT]->fillPoint->Fill(PDS.getErecon(),weight);
	if(PDS.fPID == PID_BETA || PDS.fPID == PID_MUON)
		qPassesWC[s]->fillPoint->Fill(PDS.getErecon(),weight);
	if(PDS.fPID != PID_BETA) return;
	if(PDS.passesPositionCut(s) && PDS.fType <= TYPE_III_EVENT) {
		qEnergySpectra[s][nBetaTubes][PDS.fType]->fillPoint->Fill(PDS.getErecon(),weight);
		if(PDS.fType >= TYPE_II_EVENT) {
			q23ProbCut[s][TYPE_II_EVENT]->fillPoint->Fill(PDS.getErecon(),weight*(1.0-PDS.fProbIII));
			q23ProbCut[s][TYPE_III_EVENT]->fillPoint->Fill(PDS.getErecon(),weight*PDS.fProbIII);
		}
		for(unsigned int t=0; t<nBetaTubes; t++)
			qEnergySpectra[s][t][PDS.fType]->fillPoint->Fill(PDS.scints[s].tuben[t].x,weight);
	}
}

void AsymmetryPlugin::fitAsym(float fmin, float fmax, unsigned int color, bool avg) {
	Stringmap m;
	TF1* fitter = avg?&averagerFit:&asymmetryFit;
	fitter->SetParameter(0,A0_PDG);
	fitter->SetLineColor(color);
	fitter->SetLineStyle(2);
	fitter->SetLineWidth(2.0);
	hAsym->Fit(fitter,(fmin==0 && myA->isSimulated?"Q":"QN"),"",fmin,fmax);
	
	AnaNumber AN("asym_"+ctos(choiceLetter(anChoice))+"_"+itos(fmin)+"-"+itos(fmax));
	if(avg) AN.name += "_avg";
	AN.value = fitter->GetParameter(0);	// asymmetry fit value
	AN.err = fitter->GetParError(0);	// asymmetry fit error
	myA->uploadAnaNumber(AN, GV_OPEN, AFP_OTHER);

	AN.name += "_chisq";
	AN.value = fitter->GetChisquare();	// asymmetry fit Chi^2
	AN.err = fitter->GetNDF();			// asymmetry fit number of degrees of freedom
	myA->uploadAnaNumber(AN, GV_OPEN, AFP_OTHER);
}

void AsymmetryPlugin::fitInstAsym(float fmin, float fmax, unsigned int color) {
	Stringmap m;
	averagerFit.SetLineColor(color);
	hInstAsym->Fit(&averagerFit,"Q+","",fmin,fmax);
	
	AnaNumber AN("instasym_"+ctos(choiceLetter(anChoice))+"_"+itos(fmin)+"-"+itos(fmax));
	AN.value = averagerFit.GetParameter(0);	// instrumental asymmetry average value
	AN.err = averagerFit.GetParError(0);	// uncertainty on instrumental asymmetry average value
	myA->uploadAnaNumber(AN, GV_OPEN, AFP_OTHER);
	
	AN.name += "_chisq";
	AN.value = averagerFit.GetChisquare();	// instrumental asymmetry fit Chi^2
	AN.err = averagerFit.GetNDF();			// instrumental asymmetry fit ndf
	myA->uploadAnaNumber(AN, GV_OPEN, AFP_OTHER);
}


void AsymmetryPlugin::endpointFits() {
	const float fitStart = 150;
	const float fitEnd = 635;
	for(Side s = EAST; s <= WEST; ++s) {
		for(AFPState afp = AFP_OFF; afp <= AFP_ON; ++afp) {
			for(unsigned int t=0; t<=nBetaTubes; t++) {
				float_err ep = kurieIterator((TH1F*)qEnergySpectra[s][t][TYPE_0_EVENT]->fgbg[afp]->h[1],
											 800., NULL, neutronBetaEp, fitStart, fitEnd);
											 
				AnaNumber AN("kurie_"+itos(fitStart)+"-"+itos(fitEnd));
				AN.etypes.insert(TYPE_0_EVENT);
				AN.s = s;			// side
				AN.n = t;			// PMT number (4 = combined energy)
				AN.value = ep.x;	// Kurie endpoint
				AN.err = ep.err;	// estimated endpoint uncertainty (may not be accurate)
				myA->uploadAnaNumber(AN, GV_OPEN, afp);
			}
		}
	}
}

void AsymmetryPlugin::calcSuperCombos() {
	// build total spectra based on analysis choice
	etypes.clear();
	for(Side s = EAST; s <= WEST; ++s) {
		qTotalSpectrum[s] = myA->cloneQuadHist(qEnergySpectra[s][nBetaTubes][TYPE_0_EVENT], "hTotalEvents", "All Events Energy");
		if(!(ANCHOICE_A <= anChoice && anChoice <= ANCHOICE_E))
			*qTotalSpectrum[s] *= 0;	// analysis choices without Type 0 events
		else etypes.insert(TYPE_0_EVENT);
		if(ANCHOICE_A <= anChoice && anChoice <= ANCHOICE_F && anChoice != ANCHOICE_D) {
			*qTotalSpectrum[s] += *qEnergySpectra[s][nBetaTubes][TYPE_I_EVENT];
			etypes.insert(TYPE_I_EVENT);
		}
		if(anChoice == ANCHOICE_A || anChoice == ANCHOICE_G) {
			*qTotalSpectrum[s] += *qEnergySpectra[s][nBetaTubes][TYPE_II_EVENT];
			*qTotalSpectrum[s] += *qEnergySpectra[s][nBetaTubes][TYPE_III_EVENT];
			etypes.insert(TYPE_II_EVENT);
			etypes.insert(TYPE_III_EVENT);
		}
		if(anChoice == ANCHOICE_C || anChoice == ANCHOICE_H) {
			*qTotalSpectrum[s] += *qEnergySpectra[otherSide(s)][nBetaTubes][TYPE_II_EVENT];
			*qTotalSpectrum[s] += *qEnergySpectra[s][nBetaTubes][TYPE_III_EVENT];
			etypes.insert(TYPE_II_EVENT);
			etypes.insert(TYPE_III_EVENT);
		}
		if(anChoice == ANCHOICE_E || anChoice == ANCHOICE_I) {
			*qTotalSpectrum[s] += *q23ProbCut[otherSide(s)][TYPE_II_EVENT];
			*qTotalSpectrum[s] += *q23ProbCut[s][TYPE_III_EVENT];
			etypes.insert(TYPE_II_EVENT);
			etypes.insert(TYPE_III_EVENT);
		}
		if(anChoice == ANCHOICE_J) {
			*qTotalSpectrum[s] += *qEnergySpectra[otherSide(s)][nBetaTubes][TYPE_II_EVENT];
			etypes.insert(TYPE_II_EVENT);
		}
		if(anChoice == ANCHOICE_K) {
			*qTotalSpectrum[s] += *qEnergySpectra[s][nBetaTubes][TYPE_III_EVENT];
			etypes.insert(TYPE_III_EVENT);
		}
	}
}

void AsymmetryPlugin::calculateResults() {
	
	calcSuperCombos();
	
	// calculate SR and SS
	hAsym = calculateSR("Total_Events_SR",qTotalSpectrum[EAST],qTotalSpectrum[WEST]);
	hInstAsym = calculateSR("Total_Instrumental_Asym",qTotalSpectrum[EAST],qTotalSpectrum[WEST],true,true);
	for(GVState gv = GV_CLOSED; gv <= GV_OPEN; ++gv)
		hSuperSum[gv] = calculateSuperSum("Total_Events_SuperSum",qTotalSpectrum[EAST],qTotalSpectrum[WEST],gv);
	for(EventType tp = TYPE_0_EVENT; tp <= TYPE_III_EVENT; ++tp) {
		hTpAsym[tp] = calculateSR("Asymmetry_Type_"+itos(tp),
									qEnergySpectra[tp==TYPE_II_EVENT?WEST:EAST][nBetaTubes][tp],
									qEnergySpectra[tp==TYPE_II_EVENT?EAST:WEST][nBetaTubes][tp]);
		hTpAsym[tp]->SetMinimum(-0.10);
		hTpAsym[tp]->SetMaximum(0.0);
		hEvtSS[tp] = calculateSuperSum("SuperSum_Type_"+itos(tp),
										qEnergySpectra[EAST][nBetaTubes][tp],
										qEnergySpectra[WEST][nBetaTubes][tp]);
	}
	
	// perform data fits
	fitAsym(200,675,1,true);	// match Robby's analysis
	fitAsym(50,800,7);
	fitAsym(225,675,6);
	fitAsym(0,1000,1);
	fitInstAsym();
	endpointFits();
	if(myA->grouping >= GROUP_OCTET) {
		for(Side s = EAST; s <= WEST; ++s) {
			for(EventType tp = TYPE_0_EVENT; tp <= TYPE_IV_EVENT; ++tp)
				fitHighEnergyExcess(myA,qEnergySpectra[s][nBetaTubes][tp],1000,7000);
		}
	}
	
	//------------------------------
	// counts and counts asymmetries
	//------------------------------
	
	double emin = 230;
	double emax = 660;
	
	AnaNumber ANr("asym_energy_window"); // energy window used for asymmetries/counts
	ANr.value = emin;	// lower bound
	ANr.err = emax;		// upper bound
	myA->uploadAnaNumber(ANr, GV_OTHER, AFP_OTHER);
	
	// raw counts asymmetry
	double cts[2][2];
	double cterr[2][2];
	for(Side s = EAST; s <= WEST; ++s) {
		for(AFPState afp = AFP_OFF; afp <= AFP_ON; ++afp) {
			TH1* h = (TH1F*)qTotalSpectrum[s]->fgbg[afp]->h[GV_OPEN];
			int b0 = h->FindBin(emin+0.5);
			int b1 = h->FindBin(emax-0.5);
			cts[s][afp] = h->IntegralAndError(b0,b1,cterr[s][afp])/myA->totalTime[afp][GV_OPEN][s];
			cterr[s][afp] /= myA->totalTime[afp][GV_OPEN][s];
		}
	}
	
	double S = (cts[EAST][AFP_OFF]*cts[WEST][AFP_ON])/(cts[EAST][AFP_ON]*cts[WEST][AFP_OFF]);
	AnaNumber ANa("raw_count_asym");
	ANa.value = (1-sqrt(S))/(1+sqrt(S));	// asymmetry of raw counts in energy window
	ANa.err = 0;
	for(Side s = EAST; s <= WEST; ++s)
		for(AFPState afp = AFP_OFF; afp <= AFP_ON; ++afp)
			ANa.err += pow(cterr[s][afp]/cts[s][afp],2);
	ANa.err = sqrt(ANa.err);
	ANa.err *= sqrt(S)/pow(1+sqrt(S),2);	// uncertainty on asymmetry
	myA->uploadAnaNumber(ANa, GV_OTHER, AFP_OTHER);
	
	// event count results
	AnaNumber ANc("ecut_evt_counts");
	for(Side s = EAST; s <= WEST; ++s) {
		for(EventType tp = TYPE_0_EVENT; tp <= TYPE_III_EVENT; ++tp) {
			for(AFPState afp = AFP_OFF; afp <= AFP_ON; ++afp) {
				for(GVState gv = GV_CLOSED; gv <= GV_OPEN; ++gv) {
					
					if(myA->isSimulated && gv==GV_CLOSED) continue;
					
					ANc.s = s;				// event count side
					ANc.etypes.clear();
					ANc.etypes.insert(tp);	// event backscatter type
					
					TH1* h = (TH1F*)qEnergySpectra[s][nBetaTubes][tp]->fgbg[afp]->h[gv];
					int b0 = h->FindBin(emin+0.5);
					int b1 = h->FindBin(emax-0.5);
					
					Double_t ierr;
					double gvscale = gv==GV_OPEN? 1.0 : myA->totalTime[afp][GV_OPEN][BOTH]/myA->totalTime[afp][GV_CLOSED][BOTH];
					ANc.value = h->IntegralAndError(b0,b1,ierr)*gvscale;	// number of events (time-scaled for background)
					ANc.err = ierr*gvscale;									// uncertainty on event count
					myA->uploadAnaNumber(ANc, gv, afp);
				}
			}
		}
	}
}

void AsymmetryPlugin::makePlots() {
	
	myA->defaultCanvas->SetRightMargin(0.04);
	myA->defaultCanvas->SetLeftMargin(0.12);
	
	// overall asymmetry
	hAsym->SetMinimum(-0.10);
	hAsym->SetMaximum(0.0);
	hAsym->GetXaxis()->SetRangeUser(0,800);
	hAsym->Draw();
	printCanvas("Asymmetry");
	
	// instrumental asymmetry
	hInstAsym->SetMinimum(-0.10);
	hInstAsym->SetMaximum(0.10);
	hInstAsym->GetXaxis()->SetRangeUser(0,800);
	hInstAsym->Draw();
	printCanvas("InstAsym");
	
	// super sum spectrum
	hSuperSum[GV_OPEN]->Draw();
	printCanvas("SuperSum");
	
	// energy spectra event counts
	for(EventType tp=TYPE_0_EVENT; tp<=TYPE_IV_EVENT; ++tp) {
		drawQuadSides(qEnergySpectra[EAST][nBetaTubes][tp], qEnergySpectra[WEST][nBetaTubes][tp], true, "Energy");
		for(GVState gv = GV_CLOSED; gv <= GV_OPEN; ++gv) {
			std::vector<TH1*> hToPlot;
			for(Side s = EAST; s <= WEST; ++s) {
				TH1* hTpRt = myA->flipperSummedRate(qEnergySpectra[s][nBetaTubes][tp], gv);
				hTpRt->GetYaxis()->SetTitleOffset(1.45);
				hTpRt->Scale(1000);
				hTpRt->GetYaxis()->SetTitle("event rate [mHz/keV]");
				if(gv==GV_CLOSED && tp<TYPE_IV_EVENT) {
					hTpRt->Scale(1000);
					hTpRt->GetYaxis()->SetTitle("event rate [uHz/keV]");
				}
				hTpRt->SetTitle(("Type "+itosRN(tp)+" events"+(gv==GV_OPEN?"":" background")).c_str());
				if(tp==TYPE_IV_EVENT) hTpRt->SetTitle("Gamma events");
				hTpRt->SetMinimum(0);
				hTpRt->SetLineStyle(1+s);
#ifdef PUBLICATION_PLOTS
				hTpRt->SetLineColor(1);
#endif
				if(gv==GV_CLOSED && tp>TYPE_0_EVENT) hTpRt->Rebin(2);
				if(gv==GV_CLOSED) hTpRt->GetXaxis()->SetRange(1,hTpRt->GetXaxis()->GetNbins());
				hToPlot.push_back(hTpRt);
				
				// non-WC-vetoed for comparison
				if(tp==TYPE_IV_EVENT) {
					TH1* hTpRtWC = myA->flipperSummedRate(qPassesWC[s], gv);
					hTpRtWC->Scale(1000);
					hTpRtWC->SetLineStyle(1+s);
					hToPlot.push_back(hTpRtWC);
				}
				
			}
			drawSimulHistos(hToPlot,"HIST");
			printCanvas("Energy/EnergyRate_"+itos(tp)+(gv==GV_OPEN?"":"_BG"));
			if(tp==TYPE_IV_EVENT && gv==GV_CLOSED) {
				myA->defaultCanvas->SetLogy(true);
				for(Side s = EAST; s <= WEST; ++s) { hToPlot[2*s]->SetMinimum(0.3); hToPlot[2*s+1]->SetMinimum(0.3); }
				drawSimulHistos(hToPlot,"HIST");
				printCanvas("Energy/EnergyRate_"+itos(tp)+"_BG_Log");
				myA->defaultCanvas->SetLogy(false);
			}
		}
	}
		
}

void AsymmetryPlugin::compareMCtoData(AnalyzerPlugin* AP) {
	// re-cast to correct type
	AsymmetryPlugin& dat = *(AsymmetryPlugin*)AP;
	
	myA->defaultCanvas->SetRightMargin(0.04);
	myA->defaultCanvas->SetLeftMargin(0.14);
	
	//hAsym->SetLineColor(4);
	//dat.hAsym->SetLineColor(2);
	dat.hAsym->SetTitle("Super-ratio Asymmetry");
	dat.hAsym->GetXaxis()->SetRangeUser(0,800);
	dat.hAsym->GetYaxis()->SetRangeUser(-0.07,0);
	
	hAsym->SetMarkerStyle(33);
	hAsym->SetMarkerSize(0.4);
	
	dat.hAsym->Draw("E0");
	hAsym->Draw("P SAME");
	printCanvas("DataComparison/Asymmetry");
	
	// determine scale to total events over analysis window
	TH1* simRef = hSuperSum[GV_OPEN];
	TH1* datRef = dat.hSuperSum[GV_OPEN];
	int b1 = datRef->FindBin(220);
	int b2 = datRef->FindBin(670);
	double evtscale = datRef->Integral(b1,b2)/simRef->Integral(b1,b2);
	
	// total events super sum
	dat.hSuperSum[GV_OPEN]->Scale(1000);
	dat.hSuperSum[GV_OPEN]->SetTitle("beta decay energy spectrum");
	dat.hSuperSum[GV_OPEN]->GetYaxis()->SetTitle("event rate [mHz/keV]");
	dat.hSuperSum[GV_OPEN]->SetMinimum(0);
	dat.hSuperSum[GV_OPEN]->GetYaxis()->SetTitleOffset(1.5);
	hSuperSum[GV_OPEN]->Scale(1000*evtscale);
	hSuperSum[GV_OPEN]->SetMarkerSize(0.5);
	drawDataMCPair(dat.hSuperSum[GV_OPEN],hSuperSum[GV_OPEN]);
	printCanvas("DataComparison/SuperSum");
	
	// residuals
	hSuperSum[GV_OPEN]->Add(dat.hSuperSum[GV_OPEN],-1.);
	hSuperSum[GV_OPEN]->SetTitle("MC - Data residuals");
	hSuperSum[GV_OPEN]->GetYaxis()->SetTitleOffset(1.5);
	hSuperSum[GV_OPEN]->GetYaxis()->SetTitle("event rate [mHz/keV]");
	hSuperSum[GV_OPEN]->Draw();
	drawHLine(0, myA->defaultCanvas);
	printCanvas("DataComparison/SuperSumResid");
	
	for(EventType tp = TYPE_0_EVENT; tp <= TYPE_III_EVENT; ++tp) {
		// data vs. MC supersums by event type
		dat.hEvtSS[tp]->GetXaxis()->SetRangeUser(0,800);
		hEvtSS[tp]->GetXaxis()->SetRangeUser(0,800);
		dat.hEvtSS[tp]->SetTitle(("Type "+itosRN(tp)+" event rate").c_str());
		hEvtSS[tp]->SetTitle(("Type "+itosRN(tp)+" event rate").c_str());
		dat.hEvtSS[tp]->Scale(1000);
		dat.hEvtSS[tp]->GetYaxis()->SetTitle("event rate [mHz/keV]");
		dat.hEvtSS[tp]->SetMinimum(0);
		dat.hEvtSS[tp]->GetYaxis()->SetTitleOffset(1.5);
		hEvtSS[tp]->Scale(1000*evtscale);
		hEvtSS[tp]->GetYaxis()->SetTitle("event rate [mHz/keV]");
		hEvtSS[tp]->SetMinimum(0);
		hEvtSS[tp]->GetYaxis()->SetTitleOffset(1.5);
		hEvtSS[tp]->SetMarkerSize(0.5);
		drawDataMCPair(dat.hEvtSS[tp],hEvtSS[tp]);
		printCanvas("DataComparison/SuperSum_Type_"+itos(tp));
		
		// match data, MC scales to compare shapes
		hEvtSS[tp]->Scale(dat.hEvtSS[tp]->Integral()/hEvtSS[tp]->Integral());
		hEvtSS[tp]->SetMinimum(0);
		drawDataMCPair(dat.hEvtSS[tp],hEvtSS[tp]);
		printCanvas("DataComparison/SuperSumScaled_Type_"+itos(tp));
		
		if(tp>TYPE_0_EVENT) {
			dat.hTpAsym[tp]->Rebin(4);
			dat.hTpAsym[tp]->Scale(1/4.);
			dat.hTpAsym[tp]->SetMinimum(-0.10);
			dat.hTpAsym[tp]->SetMaximum(0.0);
			
			hTpAsym[tp]->Rebin(4);
			hTpAsym[tp]->Scale(1/4.);
			hTpAsym[tp]->SetMinimum(-0.10);
			hTpAsym[tp]->SetMaximum(0.0);
		}
		dat.hTpAsym[tp]->GetXaxis()->SetRangeUser(0,800);
		hTpAsym[tp]->GetXaxis()->SetRangeUser(0,800);
		dat.hTpAsym[tp]->SetMarkerStyle(1);
		hTpAsym[tp]->SetMarkerStyle(1);
		drawHistoPair(dat.hTpAsym[tp],hTpAsym[tp]);
		printCanvas("DataComparison/Asymmetry_Type_"+itos(tp));
	}
	
	for(unsigned int t=TYPE_0_EVENT; t<=TYPE_III_EVENT; t++) {
		std::vector<TH1*> hToPlot;
		for(Side s = EAST; s <= WEST; ++s) {
			for(AFPState afp = AFP_OFF; afp <= AFP_ON; ++afp) {
				qEnergySpectra[s][nBetaTubes][t]->fgbg[afp]->h[GV_OPEN]->SetMarkerColor(2+2*s);
				qEnergySpectra[s][nBetaTubes][t]->fgbg[afp]->h[GV_OPEN]->SetMarkerStyle(22+4*afp);
				qEnergySpectra[s][nBetaTubes][t]->fgbg[afp]->h[GV_OPEN]->SetMarkerSize(0.25);
				hToPlot.push_back(qEnergySpectra[s][nBetaTubes][t]->fgbg[afp]->h[GV_OPEN]);
				dat.qEnergySpectra[s][nBetaTubes][t]->fgbg[afp]->h[GV_OPEN]->SetMarkerColor(2+2*s);
				dat.qEnergySpectra[s][nBetaTubes][t]->fgbg[afp]->h[GV_OPEN]->SetMarkerStyle(20+4*afp);
				dat.qEnergySpectra[s][nBetaTubes][t]->fgbg[afp]->h[GV_OPEN]->SetMarkerSize(0.25);
				hToPlot.push_back(dat.qEnergySpectra[s][nBetaTubes][t]->fgbg[afp]->h[GV_OPEN]);
			}
		}
		drawSimulHistos(hToPlot,"HIST P");
		printCanvas("DataComparison/Type_"+itos(t));
	}
}

