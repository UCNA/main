#include "SimAsymmetryAnalyzer.hh"
#include "BetaSpectrum.hh"
#include "GraphicsUtils.hh"
#include <TProfile.h>
#include <TGraph.h>

SimAsymmetryAnalyzer::SimAsymmetryAnalyzer(OctetAnalyzer* OA): OctetAnalyzerPlugin(OA,"simasymmetry") {
	myA->isSimulated = true;
	int nEnergyBins = 150;
	double energyMax = 1500;
	
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
			
			qWrongSide[s][t] = registerCoreHist("hWrongSide_"+itos(t),"Wrong-side ID events",
												nEnergyBins, 0, energyMax, s);
			qWrongSide[s][t]->setAxisTitle(X_DIRECTION,"Reconstructed Energy [keV]");
			qWrongSide[s][t]->setDrawRange(0,false);
			qWrongSide[s][t]->setRangeUser(0,800);
		}
	}
}

void SimAsymmetryAnalyzer::fillCoreHists(ProcessedDataScanner& PDS, double weight) {
	assert(PDS.isSimulated());
	Sim2PMT& S2P = (Sim2PMT&)PDS;
	Side s = S2P.fSide;
	if(!(s==EAST || s==WEST)) return;
	if(S2P.fPID != PID_BETA) return;
	if(S2P.passesPositionCut(s) && S2P.fType <= TYPE_III_EVENT) {
		((TProfile*)(qBCT[S2P.primSide][S2P.fType]->fillPoint))->Fill(S2P.getEtrue(),beta(S2P.ePrim)*S2P.costheta,weight);
		if( (S2P.fType == TYPE_II_EVENT?otherSide(s):s) != S2P.primSide )
			qWrongSide[s][S2P.fType]->fillPoint->Fill(S2P.getEtrue(),weight);
	}
}

void SimAsymmetryAnalyzer::makePlots() {
	drawQuad(qMissedSpectrum,"Energy/");
	for(EventType t=TYPE_0_EVENT; t<=TYPE_IV_EVENT; ++t)
		drawQuadSides(qBCT[EAST][t],qBCT[WEST][t],false,"BetaCosTheta/");
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

std::vector<TH1*> SimAsymmetryAnalyzer::calculateCorrections(AsymmetryAnalyzer& Adat, AsymmetryAnalyzer& Asim) {
	// unmodified asymmetry
	Adat.calcSuperCombos();
	std::vector<TH1*> asymStages;
	asymStages.push_back(calculateSR("BaseAsymmetry",Adat.qTotalSpectrum[EAST],Adat.qTotalSpectrum[WEST]));
	
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
				hMisID->SetMaximum(t<TYPE_II_EVENT?0.2:1);
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
	
	return asymStages;
}
