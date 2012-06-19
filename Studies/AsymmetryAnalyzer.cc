#include "AsymmetryAnalyzer.hh"
#include "BetaSpectrum.hh"
#include "KurieFitter.hh"
#include "GraphicsUtils.hh"
#include <TH2F.h>

std::string AsymmetryAnalyzer::processedLocation = "";	// set this later depending on situtation

/// nominal asymmetry for fit
Double_t asymmetryFitFunc(const Double_t* x, const Double_t* par) {
	return par[0]*plainAsymmetry(x[0],0.5)/A0_PDG;
}

TF1 AsymmetryAnalyzer::asymmetryFit = TF1("asymFit",&asymmetryFitFunc,0,neutronBetaEp,1);
TF1 AsymmetryAnalyzer::averagerFit = TF1("averagerFir","pol0",0,neutronBetaEp);
AnalysisChoice  AsymmetryAnalyzer::anChoice = ANCHOICE_A;

AsymmetryAnalyzer::AsymmetryAnalyzer(OutputManager* pnt, const std::string& nm, const std::string& inflname): OctetAnalyzer(pnt,nm,inflname) {
	for(Side s = EAST; s <= WEST; ++s) {
		qAnodeCal[s] = registerCoreHist("AnodeCal","Anode Calibration Events",50, 0, 8, s, &hAnodeCal[s]);
		
		for(EventType t=TYPE_0_EVENT; t<=TYPE_IV_EVENT; ++t) {
			for(unsigned int p=0; p<=nBetaTubes; p++) {
				qEnergySpectra[s][p][t] = registerCoreHist(std::string("hEnergy_")+(p<nBetaTubes?itos(p)+"_":"")+"Type_"+itos(t),
														   std::string("Type ")+itos(t)+" Events Energy",
														   100, 0, 1000, s, &hEnergySpectra[s][p][t]);
				qEnergySpectra[s][p][t].setAxisTitle(X_DIRECTION,"Energy [keV]");
			}
			if(t>TYPE_III_EVENT) continue;
			TH2F* hPositionsTemplate = new TH2F((std::string("hPos_Type_")+itos(t)).c_str(),
												(std::string("Type ")+itos(t)+" Positions").c_str(),
												200,-65,65,200,-65,65);
			qPositions[s][t] = registerCoreHist(*hPositionsTemplate,s,(TH1**)&hPositions[s][t]);
			qPositions[s][t].setAxisTitle(X_DIRECTION,"x Position [mm]");
			qPositions[s][t].setAxisTitle(Y_DIRECTION,"y Position [mm]");
			delete(hPositionsTemplate);
		}
	}
}

void AsymmetryAnalyzer::fillCoreHists(ProcessedDataScanner& PDS, double weight) {
	Side s = PDS.fSide;
	if(!(s==EAST || s==WEST)) return;
	if(PDS.fPID != PID_BETA) return;
	if(PDS.passesPositionCut(s) && PDS.fType == TYPE_0_EVENT && PDS.getEtrue()>225)
		hAnodeCal[s]->Fill(PDS.mwpcEnergy[s]/PDS.ActiveCal->wirechamberGainCorr(s,PDS.runClock.t[BOTH]),weight);
	if(PDS.passesPositionCut(s) && PDS.fType <= TYPE_IV_EVENT) {
		hEnergySpectra[s][nBetaTubes][PDS.fType]->Fill(PDS.getEtrue(),weight);
		for(unsigned int t=0; t<nBetaTubes; t++)
			hEnergySpectra[s][t][PDS.fType]->Fill(PDS.scints[s].tuben[t].x,weight);
	}
	if(PDS.fType <= TYPE_III_EVENT)
		hPositions[s][PDS.fType]->Fill(PDS.wires[s][X_DIRECTION].center,PDS.wires[s][Y_DIRECTION].center,weight);
}

AnaResult AsymmetryAnalyzer::getResultBase() const {
	AnaResult AR;
	AR.datp = isSimulated?AnaResult::G4_DATA:AnaResult::REAL_DATA;
	AR.startRun = runCounts.counts.begin()->first;
	AR.endRun = runCounts.counts.rbegin()->first;
	AR.anach = anChoice;
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
	qOut.insert("asymmetry",m);
	
	// save results
	AR.anatp = AnaResult::ANA_ASYM;
	AR.value = fitter->GetParameter(0);
	AR.err = fitter->GetParError(0);
	AnaCutSpec c;
	c.emin = fmin;
	c.emax = fmax;
	//c.radius = fiducialR;
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
	qOut.insert("instasym",m);
}


void AsymmetryAnalyzer::endpointFits() {
	const float fitStart = 250;
	const float fitEnd = 750;	
	for(Side s = EAST; s <= WEST; ++s) {		
		for(AFPState afp = AFP_OFF; afp <= AFP_ON; ++afp) {
			for(unsigned int t=0; t<=nBetaTubes; t++) {
				float_err ep = kurieIterator((TH1F*)qEnergySpectra[s][t][TYPE_0_EVENT].fgbg[afp].h[1],
											 800., NULL, neutronBetaEp, fitStart, fitEnd);
				Stringmap m;
				m.insert("fitStart",fitStart);
				m.insert("fitEnd",fitEnd);
				m.insert("afp",afp);
				m.insert("side",ctos(sideNames(s)));
				m.insert("tube",t);
				m.insert("type",TYPE_0_EVENT);
				m.insert("endpoint",ep.x);
				m.insert("dendpoint",ep.err);
				m.insert("counts",qEnergySpectra[s][t][TYPE_0_EVENT].fgbg[afp].h[1]->Integral());
				m.display("--- ");
				qOut.insert("kurieFit",m);
			}
		}
	}
}

void AsymmetryAnalyzer::anodeCalFits() {
	for(Side s = EAST; s <= WEST; ++s) {		
		for(AFPState afp = AFP_OFF; afp <= AFP_ON; ++afp) {
			TF1 fLandau("landauFit","landau",0,15);
			fLandau.SetLineColor(2+2*s);
			int fiterr = qAnodeCal[s].fgbg[afp].h[1]->Fit(&fLandau,"Q");
			Stringmap m;
			m.insert("afp",afp);
			m.insert("side",ctos(sideNames(s)));
			m.insert("fiterr",itos(fiterr));
			m.insert("height",fLandau.GetParameter(0));
			m.insert("d_height",fLandau.GetParError(0));
			m.insert("mpv",fLandau.GetParameter(1));
			m.insert("d_mpv",fLandau.GetParError(1));
			m.insert("sigma",fLandau.GetParameter(2));
			m.insert("d_sigma",fLandau.GetParError(2));
			qOut.insert("anodeCalFit",m);
		}
	}	
}

double integralError(const TH1* h, int b0, int b1) {
	double e = 0;
	for(int b = b0; b <= b1; b++) e += pow(h->GetBinError(b),2);
	return sqrt(e);
}

void AsymmetryAnalyzer::calculateResults() {
	// build total spectra based on analysis choice
	quadHists qTotalSpectrum[2];
	AnaResult ARtot = getResultBase();
	for(Side s = EAST; s <= WEST; ++s) {
		qTotalSpectrum[s] = cloneQuadHist(qEnergySpectra[s][nBetaTubes][TYPE_0_EVENT], "hTotalEvents", "All Events Energy");
		if(!(anChoice == ANCHOICE_A || anChoice == ANCHOICE_B || anChoice == ANCHOICE_C || anChoice == ANCHOICE_D))
			qTotalSpectrum[s] *= 0;	// analysis choices without Type 0 events
		else ARtot.etypes.insert(TYPE_0_EVENT);
		if(anChoice == ANCHOICE_A || anChoice == ANCHOICE_B || anChoice == ANCHOICE_C) {
			qTotalSpectrum[s] += qEnergySpectra[s][nBetaTubes][TYPE_I_EVENT];
			ARtot.etypes.insert(TYPE_I_EVENT);
		}
		if(anChoice == ANCHOICE_A || anChoice == ANCHOICE_C) {
			qTotalSpectrum[s] += qEnergySpectra[s][nBetaTubes][TYPE_II_EVENT];
			ARtot.etypes.insert(TYPE_II_EVENT);
		}
		if(anChoice == ANCHOICE_A || anChoice == ANCHOICE_C) {
			qTotalSpectrum[s] += qEnergySpectra[s][nBetaTubes][TYPE_III_EVENT];
			ARtot.etypes.insert(TYPE_III_EVENT);
		}
	}
	// calculate SR and SS
	hAsym = (TH1F*)calculateSR("Total_Events_SR",qTotalSpectrum[EAST],qTotalSpectrum[WEST]);
	hInstAsym = (TH1F*)calculateSR("Total_Instrumental_Asym",qTotalSpectrum[EAST],qTotalSpectrum[WEST],true,true);
	hSuperSum = (TH1F*)calculateSuperSum("Total_Events_SuperSum",qTotalSpectrum[EAST],qTotalSpectrum[WEST]);
	for(EventType tp = TYPE_0_EVENT; tp <= TYPE_II_EVENT; ++tp) {
		hTpAsym[tp] = (TH1F*)calculateSR(std::string("Asymmetry_Type_")+itos(tp),
										 qEnergySpectra[EAST][nBetaTubes][tp],
										 qEnergySpectra[WEST][nBetaTubes][tp]);
		hTpAsym[tp]->SetMinimum(-0.10);
		hTpAsym[tp]->SetMaximum(0.0);
		hEvtSS[tp] = (TH1F*)calculateSuperSum(std::string("SuperSum_Type_")+itos(tp),
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
	anodeCalFits();
	makeRatesSummary();
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

	// event count results
	AnaResult ARtot = getResultBase();
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
				
				AnaCutSpec c;
				c.emin = 225;
				c.emax = 675;
				//c.radius = fiducialR;
				TH1* h = (TH1F*)qEnergySpectra[s][nBetaTubes][tp].fgbg[afp].h[GV_OPEN];
				int b0 = h->FindBin(c.emin);
				int b1 = h->FindBin(c.emax);
				
				ARtot.value = h->Integral(b0,b1);
				ARtot.err = integralError(h,b0,b1);
				ARtot.csid = ADB->uploadCutSpec(c);
				ADB->uploadAnaResult(ARtot);
			}
		}
	}
}

void AsymmetryAnalyzer::makePlots() {
	hAsym->SetMinimum(-0.10);
	hAsym->SetMaximum(0.0);
	hAsym->Draw();
	printCanvas("Asymmetry");
	
	hInstAsym->SetMinimum(-0.10);
	hInstAsym->SetMaximum(0.10);
	hInstAsym->Draw();
	printCanvas("InstAsym");
	
	hSuperSum->Draw();
	printCanvas("SuperSum");
	
	drawQuadSides(qAnodeCal[EAST], qAnodeCal[WEST], true, "AnodeCal");
	
	for(unsigned int t=TYPE_0_EVENT; t<=TYPE_II_EVENT; t++)
		drawQuadSides(qEnergySpectra[EAST][nBetaTubes][t], qEnergySpectra[WEST][nBetaTubes][t], true, "Energy");
	
	// positions
	if(depth <= 0) {
		for(Side s = EAST; s <= WEST; ++s) {
			for(unsigned int t=TYPE_0_EVENT; t<=TYPE_II_EVENT; t++) {
				qPositions[s][t].setDrawMinimum(0);
				drawQuad(qPositions[s][t],"Positions","COL");
			}
		}
	}
}
