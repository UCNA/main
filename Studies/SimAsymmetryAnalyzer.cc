#include "SimAsymmetryAnalyzer.hh"
#include "BetaSpectrum.hh"
#include "GraphicsUtils.hh"
#include <TProfile.h>
#include <TGraph.h>
#include <iostream>

SimAsymmetryAnalyzer::SimAsymmetryAnalyzer(OctetAnalyzer* OA): OctetAnalyzerPlugin(OA,"simasymmetry") {
	myA->isSimulated = true;
	int nEnergyBins = 150;
	double energyMax = 1500;
	
	myA->ignoreMissingHistos = true;
	
	qMissedSpectrum = registerCoreHist("MissedSpectrum","Missing Events Energy Spectrum",
									   nEnergyBins, 0, energyMax, BOTH);
	qMissedSpectrum->setAxisTitle(X_DIRECTION,"Energy [keV]");
	for(Side s = EAST; s <= WEST; ++s) {
		for(EventType t=TYPE_0_EVENT; t<=TYPE_IV_EVENT; ++t) {
			
			TProfile pBCTTemplate(("pBcT_Type_"+itos(t)).c_str(),
								  ("Type "+itos(t)+" Beta Cos Theta").c_str(),
								  nEnergyBins,0,energyMax);
			qBCT[s][t] = registerCoreHist(pBCTTemplate,s);
			qBCT[s][t]->setAxisTitle(X_DIRECTION,"Reconstructed Energy [keV]");
			qBCT[s][t]->setAxisTitle(Y_DIRECTION,"<#beta cos #theta>");
			qBCT[s][t]->setDrawRange(-1.0,false);
			qBCT[s][t]->setDrawRange(1.0,true);
			qBCT[s][t]->setRangeUser(0,800);
			
			TProfile pBetaTemplate(("pBeta_Type_"+itos(t)).c_str(),
								   ("Type "+itos(t)+" #beta").c_str(),
								   nEnergyBins,0,energyMax);
			qBeta[s][t] = registerCoreHist(pBetaTemplate,s);
			qBeta[s][t]->setAxisTitle(X_DIRECTION,"Reconstructed Energy [keV]");
			qBeta[s][t]->setAxisTitle(Y_DIRECTION,"<#beta>");
			qBeta[s][t]->setDrawRange(0.0,false);
			qBeta[s][t]->setDrawRange(1.0,true);
			qBeta[s][t]->setRangeUser(0,800);
			
			TProfile pCosthTemplate(("pCosth_Type_"+itos(t)).c_str(),
									("Type "+itos(t)+" cos #theta").c_str(),
									nEnergyBins,0,energyMax);
			qCosth[s][t] = registerCoreHist(pCosthTemplate,s);
			qCosth[s][t]->setAxisTitle(X_DIRECTION,"Reconstructed Energy [keV]");
			qCosth[s][t]->setAxisTitle(Y_DIRECTION,"<cos #theta>");
			qCosth[s][t]->setDrawRange(-1.0,false);
			qCosth[s][t]->setDrawRange(1.0,true);
			qCosth[s][t]->setRangeUser(0,800);
			
			TProfile pAsymAccTemplate(("pAsymAcc_Type_"+itos(t)).c_str(),
									  ("Type "+itos(t)+" Asymmetry Acceptance").c_str(),
									  nEnergyBins,0,energyMax);
			qAsymAcc[s][t] = registerCoreHist(pAsymAccTemplate,s);
			qAsymAcc[s][t]->setAxisTitle(X_DIRECTION,"Reconstructed Energy [keV]");
			qAsymAcc[s][t]->setAxisTitle(Y_DIRECTION,"<A(E,#theta)>");
			qAsymAcc[s][t]->setDrawRange(-0.1,false);
			qAsymAcc[s][t]->setDrawRange(0.1,true);
			qAsymAcc[s][t]->setRangeUser(0,800);
			
			
			qWrongSide[s][t] = registerCoreHist("hWrongSide_"+itos(t),"Wrong-side ID events",
												nEnergyBins, 0, energyMax, s);
			qWrongSide[s][t]->setAxisTitle(X_DIRECTION,"Reconstructed Energy [keV]");
			qWrongSide[s][t]->setDrawRange(0,false);
			qWrongSide[s][t]->setRangeUser(0,800);
		}
	}
	
	myA->ignoreMissingHistos = false;
}

void SimAsymmetryAnalyzer::fillCoreHists(ProcessedDataScanner& PDS, double weight) {
	assert(PDS.isSimulated());
	Sim2PMT& S2P = (Sim2PMT&)PDS;
	Side s = S2P.fSide;
	if(!(s==EAST || s==WEST)) return;
	if(S2P.fPID != PID_BETA) return;
	if(S2P.passesPositionCut(s) && S2P.fType <= TYPE_III_EVENT) {
		// NOTE: these profiles filled WITHOUT weight to capture isotropic acceptance
		((TProfile*)(qBCT[S2P.primSide][S2P.fType]->fillPoint))->Fill(S2P.getEtrue(),beta(S2P.ePrim)*S2P.costheta);
		((TProfile*)(qBeta[S2P.primSide][S2P.fType]->fillPoint))->Fill(S2P.getEtrue(),beta(S2P.ePrim));
		((TProfile*)(qCosth[S2P.primSide][S2P.fType]->fillPoint))->Fill(S2P.getEtrue(),S2P.costheta);
		((TProfile*)(qAsymAcc[S2P.primSide][S2P.fType]->fillPoint))->Fill(S2P.getEtrue(),correctedAsymmetry(S2P.ePrim,S2P.costheta));
		if( (S2P.fType == TYPE_II_EVENT?otherSide(s):s) != S2P.primSide)
			qWrongSide[s][S2P.fType]->fillPoint->Fill(S2P.getEtrue(),weight);
	}
}

void SimAsymmetryAnalyzer::makePlots() {
	drawQuad(qMissedSpectrum,"Energy/");
	for(EventType t=TYPE_0_EVENT; t<=TYPE_III_EVENT; ++t) {
		drawQuadSides(qBCT[EAST][t],qBCT[WEST][t],false,"BetaCosTheta/");
		drawQuadSides(qBeta[EAST][t],qBeta[WEST][t],false,"BetaCosTheta/");
		drawQuadSides(qCosth[EAST][t],qCosth[WEST][t],false,"BetaCosTheta/");
		drawQuadSides(qAsymAcc[EAST][t],qAsymAcc[WEST][t],false,"BetaCosTheta/");
	}
	for(EventType t=TYPE_0_EVENT; t<=TYPE_III_EVENT; ++t)
		drawQuadSides(qWrongSide[EAST][t],qWrongSide[WEST][t],false,"WrongSide/");
}

void combineHists(TH1* h1, TH1* h2) {
	h1->Add(h2);
	h2->Scale(0.);
	h2->Add(h1);
}

void unmix(TH1* hE, TH1* hW, TH1* pEx, TH1* pWx) {
	for(int b=0; b<=hE->GetNbinsX()+1; b++) {
		double xE = pEx->GetBinContent(b);
		double xW = pWx->GetBinContent(b);
		double nE = hE->GetBinContent(b);
		double nW = hW->GetBinContent(b);
		hE->SetBinContent(b,(nE-nW*xW)*(1+xE)/(1-xE*xW));
		hW->SetBinContent(b,(nW-nE*xE)*(1+xW)/(1-xE*xW));
	}
}

std::vector<TH1*> SimAsymmetryAnalyzer::calculateCorrectionsOld(AsymmetryAnalyzer& Asim) {
	// apply backscatter reweightings... ???
	
	// unmodified asymmetry
	Asim.calcSuperCombos();
	std::vector<TH1*> asymStages;
	asymStages.push_back(calculateSR("BaseAsymmetry",Asim.qTotalSpectrum[EAST],Asim.qTotalSpectrum[WEST]));
	
	// mangle Type II's (turn into mis-ID's IIIs) for trigger-side analysis choices
	if(Asim.anChoice==ANCHOICE_A || Asim.anChoice==ANCHOICE_G) {
		for(Side s = EAST; s <= WEST; ++s) {
			*Asim.qEnergySpectra[s][nBetaTubes][TYPE_III_EVENT] += *Asim.qEnergySpectra[s][nBetaTubes][TYPE_II_EVENT];
			*qWrongSide[s][TYPE_II_EVENT] *= -1;
			*Asim.qEnergySpectra[s][nBetaTubes][TYPE_II_EVENT] += *qWrongSide[s][TYPE_II_EVENT];
			*qWrongSide[s][TYPE_III_EVENT] += *Asim.qEnergySpectra[s][nBetaTubes][TYPE_II_EVENT];
			*Asim.qEnergySpectra[s][nBetaTubes][TYPE_II_EVENT] *= 0;
			*qWrongSide[s][TYPE_II_EVENT] *= 0;
			
			*qAsymAcc[s][TYPE_III_EVENT] += *qAsymAcc[otherSide(s)][TYPE_II_EVENT];
			*qAsymAcc[s][TYPE_II_EVENT] *= 0;
			*qBCT[s][TYPE_III_EVENT] += *qBCT[otherSide(s)][TYPE_II_EVENT];
			*qBCT[s][TYPE_II_EVENT] *= 0;
		}
	}
	
	// fix event type assignments based on MC wrong-side assignments
	for(EventType t=TYPE_0_EVENT; t<=TYPE_III_EVENT; ++t) {
		quadHists* qMisCorr[2];
		quadHists* qRightID[2];
		std::vector<TH1*> hToPlot;
		
		for(Side s = EAST; s <= WEST; ++s) {
			qMisCorr[s] = myA->cloneQuadHist(qWrongSide[otherSide(s)][t],sideSubst("SideCorrector_%c_",s)+itos(t),"Mis-ID events correction");
			qRightID[s] = myA->cloneQuadHist(Asim.qEnergySpectra[s][nBetaTubes][t],sideSubst("RightID_%c_",s)+itos(t),"Correctly-ID'd events");
		}
		
		for(Side s = EAST; s <= WEST; ++s) {
			for(AFPState a = AFP_OFF; a<=AFP_ON; ++a) {
				TH1* hMisID = qMisCorr[s]->fgbg[a]->h[GV_OPEN];
				TH1* hRightID = qRightID[s]->fgbg[a]->h[GV_OPEN];
				hRightID->Add(qWrongSide[s][t]->fgbg[a]->h[GV_OPEN],-1.0); // remove counts coming from wrong side
				hMisID->Divide(hRightID);
				
				hMisID->SetLineColor(2+2*s);
				hMisID->SetLineStyle(1+2*a);
				hMisID->SetMinimum(0);
				hMisID->SetMaximum(t<TYPE_II_EVENT?0.2:Asim.anChoice==ANCHOICE_A?2:1);
				hToPlot.push_back(hMisID);
			}
		}
		drawSimulHistos(hToPlot);
		printCanvas("MisId_"+itos(t));
		
		for(AFPState a = AFP_OFF; a<=AFP_ON; ++a) {
			unmix(Asim.qEnergySpectra[EAST][nBetaTubes][t]->fgbg[a]->h[GV_OPEN],
				  Asim.qEnergySpectra[WEST][nBetaTubes][t]->fgbg[a]->h[GV_OPEN],
				  qMisCorr[EAST]->fgbg[a]->h[GV_OPEN],
				  qMisCorr[WEST]->fgbg[a]->h[GV_OPEN]);
		}
		// at this point, all events for this type should be assigned to the correct side; re-calculate intermediate asymmetry
		Asim.calcSuperCombos();
		asymStages.push_back(calculateSR("SR_Delta_2_"+itos(t),Asim.qTotalSpectrum[EAST],Asim.qTotalSpectrum[WEST]));
	}
	
	// side-summed beta cos theta acceptance
	quadHists* qBCTsum[2];
	for(Side s = EAST; s <= WEST; ++s) {
		qBCTsum[s] = myA->cloneQuadHist(qBCT[s][TYPE_0_EVENT], sideSubst("BCTSum_%c",s), "Beta Cos Theta acceptance");
		if(!(ANCHOICE_A <= Asim.anChoice && Asim.anChoice <= ANCHOICE_E))
			*qBCTsum[s] *= 0;	// analysis choices without Type 0 events
		if(ANCHOICE_A <= Asim.anChoice && Asim.anChoice <= ANCHOICE_F && Asim.anChoice != ANCHOICE_D)
			*qBCTsum[s] += *qBCT[s][TYPE_I_EVENT];
		if(Asim.anChoice == ANCHOICE_A || Asim.anChoice == ANCHOICE_G) {
			*qBCT[otherSide(s)][TYPE_II_EVENT] *= -1.;
			*qBCTsum[s] += *qBCT[otherSide(s)][TYPE_II_EVENT];
			*qBCTsum[s] += *qBCT[s][TYPE_III_EVENT];
		}
		if(Asim.anChoice == ANCHOICE_C || Asim.anChoice == ANCHOICE_H) {
			*qBCTsum[s] += *qBCT[s][TYPE_II_EVENT];
			*qBCTsum[s] += *qBCT[s][TYPE_III_EVENT];
		}
		if(Asim.anChoice == ANCHOICE_J)
			*qBCTsum[s] += *qBCT[s][TYPE_II_EVENT];
		if(Asim.anChoice == ANCHOICE_K)
			*qBCTsum[s] += *qBCT[s][TYPE_III_EVENT];
	}
	for(Side s = EAST; s <= WEST; ++s)
		for(AFPState a = AFP_OFF; a <= AFP_ON; ++a)
			if(!(s==EAST && a==AFP_OFF))
				qBCTsum[EAST]->fgbg[AFP_OFF]->h[GV_OPEN]->Add(qBCTsum[s]->fgbg[a]->h[GV_OPEN],s==EAST?1.:-1);
	qBCTsum[EAST]->fgbg[AFP_OFF]->h[GV_OPEN]->Scale(-1.0);
	qBCTsum[EAST]->fgbg[AFP_OFF]->h[GV_OPEN]->GetYaxis()->SetRangeUser(0.,0.5);
	qBCTsum[EAST]->fgbg[AFP_OFF]->h[GV_OPEN]->Draw();
	printCanvas("BetaCosThetaAcceptance");
	
	// Delta 3 correction for average beta cos theta difference
	asymStages.push_back((TH1*)asymStages.back()->Clone("SR_Delta_3"));
	for(int b=1; b<=asymStages.back()->GetNbinsX(); b++) {
		double e0 = asymStages.back()->GetBinCenter(b);
		asymStages.back()->SetBinContent(b, asymStages.back()->GetBinContent(b) *
										 (beta(e0)/2.) / qBCTsum[EAST]->fgbg[AFP_OFF]->h[GV_OPEN]->GetBinContent(b));
	}
	
	return asymStages;
	
}

std::vector<TH1*> SimAsymmetryAnalyzer::calculateCorrections(AsymmetryAnalyzer& Adat, AsymmetryAnalyzer& Asim) {
	// unmodified asymmetry
	Adat.calcSuperCombos();
	std::vector<TH1*> asymStages;
	asymStages.push_back(calculateSR("BaseAsymmetry",Adat.qTotalSpectrum[EAST],Adat.qTotalSpectrum[WEST]));
	
	// mangle Type II's (turn into mis-ID's IIIs) for trigger-side analysis choices
	if(Adat.anChoice==ANCHOICE_A || Adat.anChoice==ANCHOICE_G) {
		for(Side s = EAST; s <= WEST; ++s) {
			*Asim.qEnergySpectra[s][nBetaTubes][TYPE_III_EVENT] += *Asim.qEnergySpectra[s][nBetaTubes][TYPE_II_EVENT];
			*qWrongSide[s][TYPE_II_EVENT] *= -1;
			*Asim.qEnergySpectra[s][nBetaTubes][TYPE_II_EVENT] += *qWrongSide[s][TYPE_II_EVENT];
			*qWrongSide[s][TYPE_III_EVENT] += *Asim.qEnergySpectra[s][nBetaTubes][TYPE_II_EVENT];
			*Asim.qEnergySpectra[s][nBetaTubes][TYPE_II_EVENT] *= 0;
			*qWrongSide[s][TYPE_II_EVENT] *= 0;
			
			*qAsymAcc[s][TYPE_III_EVENT] += *qAsymAcc[otherSide(s)][TYPE_II_EVENT];
			*qAsymAcc[s][TYPE_II_EVENT] *= 0;
			*qBCT[s][TYPE_III_EVENT] += *qBCT[otherSide(s)][TYPE_II_EVENT];
			*qBCT[s][TYPE_II_EVENT] *= 0;
			
			*Adat.qEnergySpectra[s][nBetaTubes][TYPE_III_EVENT] += *Adat.qEnergySpectra[s][nBetaTubes][TYPE_II_EVENT];
			*Adat.qEnergySpectra[s][nBetaTubes][TYPE_II_EVENT] *= 0;
		}
	}
	
	// fix event type assignments based on MC wrong-side assignments
	for(EventType t=TYPE_0_EVENT; t<=TYPE_III_EVENT; ++t) {
		quadHists* qMisCorr[2];
		std::vector<TH1*> hToPlot;
		
		for(Side s = EAST; s <= WEST; ++s)
			qMisCorr[s] = myA->cloneQuadHist(qWrongSide[otherSide(s)][t],sideSubst("SideCorrector_%c_",s)+itos(t),"Mis-ID events correction");
		
		for(Side s = EAST; s <= WEST; ++s) {
			for(AFPState a = AFP_OFF; a<=AFP_ON; ++a) {
				TH1* hMisID = qMisCorr[s]->fgbg[a]->h[GV_OPEN];
				TH1* hRightID = Asim.qEnergySpectra[s][nBetaTubes][t]->fgbg[a]->h[GV_OPEN];
				hRightID->Add(qWrongSide[s][t]->fgbg[a]->h[GV_OPEN],-1.0); // remove counts coming from wrong side
				hMisID->Divide(hRightID);
				
				hMisID->SetLineColor(2+2*s);
				hMisID->SetLineStyle(1+2*a);
				hMisID->SetMinimum(0);
				hMisID->SetMaximum(t<TYPE_II_EVENT?0.2:Adat.anChoice==ANCHOICE_A?2:1);
				hToPlot.push_back(hMisID);
			}
		}
		drawSimulHistos(hToPlot);
		printCanvas("MisId_"+itos(t));
		
		for(AFPState a = AFP_OFF; a<=AFP_ON; ++a) {
			unmix(Adat.qEnergySpectra[EAST][nBetaTubes][t]->fgbg[a]->h[GV_OPEN],
				  Adat.qEnergySpectra[WEST][nBetaTubes][t]->fgbg[a]->h[GV_OPEN],
				  qMisCorr[EAST]->fgbg[a]->h[GV_OPEN],
				  qMisCorr[WEST]->fgbg[a]->h[GV_OPEN]);
		}
		// at this point, all events for this type should be assigned to the correct side; re-calculate intermediate asymmetry
		Adat.calcSuperCombos();
		asymStages.push_back(calculateSR("SR_Delta_2_"+itos(t),Adat.qTotalSpectrum[EAST],Adat.qTotalSpectrum[WEST]));
	}
	
	// Delta 3 correction for average beta cos theta difference
	for(EventType t=TYPE_0_EVENT; t<=TYPE_III_EVENT; ++t) {
		for(Side s = EAST; s <= WEST; ++s) {
			for(AFPState a = AFP_OFF; a<=AFP_ON; ++a) {
				TH1* hDat = Adat.qEnergySpectra[s][nBetaTubes][t]->fgbg[a]->h[GV_OPEN];
				TProfile* pAsAc = (TProfile*)qAsymAcc[s][t]->fgbg[a]->h[GV_OPEN];
				for(int b=1; b<=hDat->GetNbinsX(); b++) {
					double e0 = hDat->GetBinCenter(b);
					hDat->SetBinContent(b, hDat->GetBinContent(b) *
										(1-correctedAsymmetry(e0,(a==AFP_ON?-1:1)*(s==EAST?-0.5:0.5))) /
										(1-(a==AFP_ON?-1:1)*pAsAc->GetBinContent(b)));
				}
			}
		}
	}
	Adat.calcSuperCombos();
	asymStages.push_back(calculateSR("SR_Delta_3",Adat.qTotalSpectrum[EAST],Adat.qTotalSpectrum[WEST]));
	
	return asymStages;
}
