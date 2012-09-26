#include "AsymmetryAnalyzer.hh"
#include "HighEnergyExcess.hh"
#include "BetaSpectrum.hh"
#include "KurieFitter.hh"
#include "GraphicsUtils.hh"
#include <TH2F.h>

/// nominal asymmetry for fit
Double_t asymmetryFitFunc(const Double_t* x, const Double_t* par) {
	return par[0]*plainAsymmetry(x[0],0.5)/A0_PDG;
}

TF1 AsymmetryAnalyzer::asymmetryFit = TF1("asymFit",&asymmetryFitFunc,0,neutronBetaEp,1);
TF1 AsymmetryAnalyzer::averagerFit = TF1("averagerFit","pol0",0,neutronBetaEp);

AsymmetryAnalyzer::AsymmetryAnalyzer(OctetAnalyzer* OA):
OctetAnalyzerPlugin(OA,"asymmetry"), nEnergyBins(150), energyMax(1500), anChoice(ANCHOICE_C) {
	for(Side s = EAST; s <= WEST; ++s) {
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

void AsymmetryAnalyzer::fillCoreHists(ProcessedDataScanner& PDS, double weight) {
	Side s = PDS.fSide;
	if(!(s==EAST || s==WEST)) return;
	if(PDS.fPID == PID_SINGLE && PDS.fType == TYPE_IV_EVENT)
		qEnergySpectra[s][nBetaTubes][TYPE_IV_EVENT]->fillPoint->Fill(PDS.getEtrue(),weight);
	if(PDS.fPID != PID_BETA) return;
	if(PDS.passesPositionCut(s) && PDS.fType <= TYPE_III_EVENT) {
		qEnergySpectra[s][nBetaTubes][PDS.fType]->fillPoint->Fill(PDS.getEtrue(),weight);
		if(PDS.fType >= TYPE_II_EVENT) {
			q23ProbCut[s][TYPE_II_EVENT]->fillPoint->Fill(PDS.getEtrue(),weight*(1.0-PDS.fProbIII));
			q23ProbCut[s][TYPE_III_EVENT]->fillPoint->Fill(PDS.getEtrue(),weight*PDS.fProbIII);
		}
		for(unsigned int t=0; t<nBetaTubes; t++)
			qEnergySpectra[s][t][PDS.fType]->fillPoint->Fill(PDS.scints[s].tuben[t].x,weight);
	}
}

AnaResult AsymmetryAnalyzer::getResultBase() const {
	AnaResult AR;
	AR.datp = myA->isSimulated?AnaResult::G4_DATA:AnaResult::REAL_DATA;
	AR.startRun = myA->runCounts.counts.begin()->first;
	AR.endRun = myA->runCounts.counts.rbegin()->first;
	AR.anach = anChoice;
	AR.grouping = AnaResult::RunGrouping(AnaResult::GROUP_OCTET-myA->depth-myA->isSimulated);
	return AR;
}

void AsymmetryAnalyzer::fitAsym(float fmin, float fmax, unsigned int color, AnaResult AR, bool avg) {
	Stringmap m;
	TF1* fitter = avg?&averagerFit:&asymmetryFit;
	fitter->SetParameter(0,A0_PDG);
	fitter->SetLineColor(color);
	hAsym->Fit(fitter,"Q+","",fmin,fmax);
	m.insert("A0_fit",fitter->GetParameter(0));
	m.insert("dA0",fitter->GetParError(0));
	m.insert("A0_chi2",fitter->GetChisquare());
	m.insert("A0_NDF",fitter->GetNDF());
	m.insert("fitMin",fmin);
	m.insert("fitMax",fmax);
	m.insert("anChoice",itos(anChoice));
	m.insert("method",avg?"average":"fit");
	myA->qOut.insert("asymmetry",m);
	
	// save results
	AR.anatp = AnaResult::ANA_ASYM;
	AR.value = fitter->GetParameter(0);
	AR.err = fitter->GetParError(0);
	AnaCutSpec c;
	c.emin = fmin;
	c.emax = fmax;
	c.radius = 50.;
	asymFits.push_back(AR);
	asymCuts.push_back(c);
}

void AsymmetryAnalyzer::fitInstAsym(float fmin, float fmax, unsigned int color) {
	Stringmap m;
	averagerFit.SetLineColor(color);
	hInstAsym->Fit(&averagerFit,"Q+","",fmin,fmax);
	m.insert("IA",averagerFit.GetParameter(0));
	m.insert("dIA",averagerFit.GetParError(0));
	m.insert("IA_chi2",averagerFit.GetChisquare());
	m.insert("IA_NDF",averagerFit.GetNDF());
	m.insert("fitMin",fmin);
	m.insert("fitMax",fmax);
	m.insert("anChoice",itos(anChoice));
	m.insert("method","average");
	myA->qOut.insert("instasym",m);
}


void AsymmetryAnalyzer::endpointFits() {
	const float fitStart = 150;
	const float fitEnd = 700;
	for(Side s = EAST; s <= WEST; ++s) {
		for(AFPState afp = AFP_OFF; afp <= AFP_ON; ++afp) {
			for(unsigned int t=0; t<=nBetaTubes; t++) {
				float_err ep = kurieIterator((TH1F*)qEnergySpectra[s][t][TYPE_0_EVENT]->fgbg[afp]->h[1],
											 800., NULL, neutronBetaEp, fitStart, fitEnd);
				Stringmap m;
				m.insert("fitStart",fitStart);
				m.insert("fitEnd",fitEnd);
				m.insert("afp",afpWords(afp));
				m.insert("side",ctos(sideNames(s)));
				m.insert("tube",t);
				m.insert("type",TYPE_0_EVENT);
				m.insert("endpoint",ep.x);
				m.insert("dendpoint",ep.err);
				m.insert("counts",qEnergySpectra[s][t][TYPE_0_EVENT]->fgbg[afp]->h[1]->Integral());
				m.display("--- ");
				myA->qOut.insert("kurieFit",m);
			}
		}
	}
}

void AsymmetryAnalyzer::calcSuperCombos() {
	// build total spectra based on analysis choice
	ARtot = getResultBase();
	for(Side s = EAST; s <= WEST; ++s) {
		qTotalSpectrum[s] = myA->cloneQuadHist(qEnergySpectra[s][nBetaTubes][TYPE_0_EVENT], "hTotalEvents", "All Events Energy");
		if(!(ANCHOICE_A <= anChoice && anChoice <= ANCHOICE_E))
			*qTotalSpectrum[s] *= 0;	// analysis choices without Type 0 events
		else ARtot.etypes.insert(TYPE_0_EVENT);
		if(ANCHOICE_A <= anChoice && anChoice <= ANCHOICE_F && anChoice != ANCHOICE_D) {
			*qTotalSpectrum[s] += *qEnergySpectra[s][nBetaTubes][TYPE_I_EVENT];
			ARtot.etypes.insert(TYPE_I_EVENT);
		}
		if(anChoice == ANCHOICE_A || anChoice == ANCHOICE_G) {
			*qTotalSpectrum[s] += *qEnergySpectra[s][nBetaTubes][TYPE_II_EVENT];
			*qTotalSpectrum[s] += *qEnergySpectra[s][nBetaTubes][TYPE_III_EVENT];
			ARtot.etypes.insert(TYPE_II_EVENT);
			ARtot.etypes.insert(TYPE_III_EVENT);
		}
		if(anChoice == ANCHOICE_C || anChoice == ANCHOICE_H) {
			*qTotalSpectrum[s] += *qEnergySpectra[otherSide(s)][nBetaTubes][TYPE_II_EVENT];
			*qTotalSpectrum[s] += *qEnergySpectra[s][nBetaTubes][TYPE_III_EVENT];
			ARtot.etypes.insert(TYPE_II_EVENT);
			ARtot.etypes.insert(TYPE_III_EVENT);
		}
		if(anChoice == ANCHOICE_E || anChoice == ANCHOICE_I) {
			*qTotalSpectrum[s] += *q23ProbCut[otherSide(s)][TYPE_II_EVENT];
			*qTotalSpectrum[s] += *q23ProbCut[s][TYPE_III_EVENT];
			ARtot.etypes.insert(TYPE_II_EVENT);
			ARtot.etypes.insert(TYPE_III_EVENT);
		}
		if(anChoice == ANCHOICE_J) {
			*qTotalSpectrum[s] += *qEnergySpectra[otherSide(s)][nBetaTubes][TYPE_II_EVENT];
			ARtot.etypes.insert(TYPE_II_EVENT);
		}
		if(anChoice == ANCHOICE_K) {
			*qTotalSpectrum[s] += *qEnergySpectra[s][nBetaTubes][TYPE_III_EVENT];
			ARtot.etypes.insert(TYPE_III_EVENT);
		}
	}
}

void AsymmetryAnalyzer::calculateResults() {
	
	calcSuperCombos();
	
	// calculate SR and SS
	hAsym = (TH1F*)calculateSR("Total_Events_SR",qTotalSpectrum[EAST],qTotalSpectrum[WEST]);
	hInstAsym = (TH1F*)calculateSR("Total_Instrumental_Asym",qTotalSpectrum[EAST],qTotalSpectrum[WEST],true,true);
	hSuperSum = (TH1F*)calculateSuperSum("Total_Events_SuperSum",qTotalSpectrum[EAST],qTotalSpectrum[WEST]);
	for(EventType tp = TYPE_0_EVENT; tp <= TYPE_III_EVENT; ++tp) {
		hTpAsym[tp] = (TH1F*)calculateSR("Asymmetry_Type_"+itos(tp),
										 qEnergySpectra[tp==TYPE_II_EVENT?WEST:EAST][nBetaTubes][tp],
										 qEnergySpectra[tp==TYPE_II_EVENT?EAST:WEST][nBetaTubes][tp]);
		hTpAsym[tp]->SetMinimum(-0.10);
		hTpAsym[tp]->SetMaximum(0.0);
		hEvtSS[tp] = (TH1F*)calculateSuperSum("SuperSum_Type_"+itos(tp),
											  qEnergySpectra[EAST][nBetaTubes][tp],
											  qEnergySpectra[WEST][nBetaTubes][tp]);
	}
	
	// perform data fits
	fitAsym(200,675,1,ARtot,true);	// match Robby's analysis
	asymFits.pop_back(); asymCuts.pop_back();
	fitAsym(50,800,7,ARtot);
	fitAsym(225,675,6,ARtot);
	fitInstAsym();
	endpointFits();
	for(Side s = EAST; s <= WEST; ++s) {
		//highEnergyExcess(qMuonSpectra[s][true],1000,2000);
		//highEnergyExcess(qBackMuons[s][true],1000,2000);
		for(EventType tp = TYPE_0_EVENT; tp <= TYPE_IV_EVENT; ++tp)
			fitHighEnergyExcess(myA->qOut,qEnergySpectra[s][nBetaTubes][tp],1000,7000);
	}
	myA->makeRatesSummary();
}

void AsymmetryAnalyzer::uploadAnaResults() {
	
	AnalysisDB* ADB = AnalysisDB::getADB();
	
	// delete old fit results
	for(unsigned int n=0; n<asymFits.size(); n++) {
		std::vector<AnaResult> oldr = ADB->findMatching(asymFits[n]);
		for(unsigned int i=0; i<oldr.size(); i++)
			ADB->deleteAnaResult(oldr[i].arid);
	}
	// upload new fit results
	for(unsigned int n=0; n<asymFits.size(); n++) {
		asymFits[n].csid = ADB->uploadCutSpec(asymCuts[n]);
		ADB->uploadAnaResult(asymFits[n]);
	}
	
	AnaCutSpec c;
	c.emin = 220;
	c.emax = 670;
	c.radius = 50.;
	ARtot.csid = ADB->uploadCutSpec(c);

	// raw counts asymmetry
	ARtot.anatp = AnaResult::ANA_ASYM;
	double cts[2][2];
	double cterr[2][2];
	for(Side s = EAST; s <= WEST; ++s) {
		for(AFPState afp = AFP_OFF; afp <= AFP_ON; ++afp) {
			TH1* h = (TH1F*)qTotalSpectrum[s]->fgbg[afp]->h[GV_OPEN];
			int b0 = h->FindBin(c.emin+0.5);
			int b1 = h->FindBin(c.emax-0.5);
			cts[s][afp] = h->IntegralAndError(b0,b1,cterr[s][afp])/myA->totalTime[afp][GV_OPEN][s];
			cterr[s][afp] /= myA->totalTime[afp][GV_OPEN][s];
		}
	}
	
	// RWP:
	// Ar = (1-sqrt(S))/(1+sqrt(S))
	// S = (<R1><R2>)/(<R3><R4>) with the <R> being the average of the rates then the statistical uncertainty is :
	// dA = (sqrt(S)) / (1+sqrt(S))^2 sqrt(sum (dR/R)^2)
	double S = (cts[EAST][AFP_OFF]*cts[WEST][AFP_ON])/(cts[EAST][AFP_ON]*cts[WEST][AFP_OFF]);
	ARtot.value = (1-sqrt(S))/(1+sqrt(S));
	ARtot.err = 0;
	for(Side s = EAST; s <= WEST; ++s)
		for(AFPState afp = AFP_OFF; afp <= AFP_ON; ++afp)
			ARtot.err += pow(cterr[s][afp]/cts[s][afp],2);
	ARtot.err = sqrt(ARtot.err);
	ARtot.err *= sqrt(S)/pow(1+sqrt(S),2);
	ADB->uploadAnaResult(ARtot);
	
	
	// event count results
	ARtot.anatp = AnaResult::ANA_COUNTS;
	for(Side s = EAST; s <= WEST; ++s) {
		for(EventType tp = TYPE_0_EVENT; tp <= TYPE_III_EVENT; ++tp) {
			for(AFPState afp = AFP_OFF; afp <= AFP_ON; ++afp) {
				ARtot.s = s;
				ARtot.afp = afp;
				ARtot.etypes.clear();
				ARtot.etypes.insert(tp);
				
				// delete old
				std::vector<AnaResult> oldr = ADB->findMatching(ARtot);
				for(unsigned int i=0; i<oldr.size(); i++)
					ADB->deleteAnaResult(oldr[i].arid);
				
				TH1* h = (TH1F*)qEnergySpectra[s][nBetaTubes][tp]->fgbg[afp]->h[GV_OPEN];
				int b0 = h->FindBin(c.emin+0.5);
				int b1 = h->FindBin(c.emax-0.5);
				
				Double_t ierr;
				ARtot.value = h->IntegralAndError(b0,b1,ierr);
				ARtot.err = ierr;
				ADB->uploadAnaResult(ARtot);
			}
		}
	}
	
	
	
}

void AsymmetryAnalyzer::makePlots() {
	
	hAsym->SetMinimum(-0.10);
	hAsym->SetMaximum(0.0);
	hAsym->GetXaxis()->SetRangeUser(0,800);
	hAsym->Draw();
	printCanvas("Asymmetry");
	
	hInstAsym->SetMinimum(-0.10);
	hInstAsym->SetMaximum(0.10);
	hInstAsym->GetXaxis()->SetRangeUser(0,800);
	hInstAsym->Draw();
	printCanvas("InstAsym");
	
	hSuperSum->Draw();
	printCanvas("SuperSum");
		
	for(unsigned int t=TYPE_0_EVENT; t<=TYPE_IV_EVENT; t++)
		drawQuadSides(qEnergySpectra[EAST][nBetaTubes][t], qEnergySpectra[WEST][nBetaTubes][t], true, "Energy");
}

void AsymmetryAnalyzer::compareMCtoData(AnalyzerPlugin* AP) {
	// re-cast to correct type
	AsymmetryAnalyzer& dat = *(AsymmetryAnalyzer*)AP;
	
	hAsym->GetXaxis()->SetRangeUser(0,800);
	hAsym->SetLineColor(4);
	dat.hAsym->SetLineColor(2);
	hAsym->Draw("HIST E1");
	dat.hAsym->Draw("SAME HIST E1");
	printCanvas("DataComparison/Asymmetry");
	
	drawHistoPair(dat.hSuperSum,hSuperSum);
	printCanvas("DataComparison/SuperSum");
	double evtscale = dat.hEvtSS[TYPE_0_EVENT]->Integral()/hEvtSS[TYPE_0_EVENT]->Integral();
	for(EventType tp = TYPE_0_EVENT; tp <= TYPE_III_EVENT; ++tp) {
		// data vs. MC supersums by event type
		dat.hEvtSS[tp]->GetXaxis()->SetRangeUser(0,800);
		hEvtSS[tp]->GetXaxis()->SetRangeUser(0,800);
		dat.hEvtSS[tp]->SetMinimum(0);
		hEvtSS[tp]->SetMinimum(0);
		dat.hEvtSS[tp]->SetMarkerStyle(1);
		hEvtSS[tp]->SetMarkerStyle(1);
		hEvtSS[tp]->Scale(evtscale);
		drawHistoPair(dat.hEvtSS[tp],hEvtSS[tp]);
		printCanvas("DataComparison/SuperSum_Type_"+itos(tp));
		// match data, MC scales to compare shapes
		hEvtSS[tp]->Scale(dat.hEvtSS[tp]->Integral()/hEvtSS[tp]->Integral());
		hEvtSS[tp]->SetMinimum(0);
		drawHistoPair(dat.hEvtSS[tp],hEvtSS[tp]);
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

