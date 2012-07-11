#include "SimAsymmetryAnalyzer.hh"
#include "BetaSpectrum.hh"
#include "GraphicsUtils.hh"
#include <TProfile.h>

SimAsymmetryAnalyzer::SimAsymmetryAnalyzer(OutputManager* pnt, const std::string& nm, const std::string& inflname):
AsymmetryAnalyzer(pnt,nm,inflname) {
	isSimulated = true;
	qMissedSpectrum = registerCoreHist("MissedSpectrum","Missing Events Energy Spectrum",
									   nEnergyBins, 0, energyMax, BOTH);
	for(Side s = EAST; s <= WEST; ++s) {		
		for(EventType t=TYPE_0_EVENT; t<=TYPE_IV_EVENT; ++t) {
			TProfile pBCTTemplate(("pBcT_Type_"+itos(t)).c_str(),
								  ("Type "+itos(t)+" Beta Cos Theta").c_str(),
								  nEnergyBins,0,energyMax);
			qBCT[s][t] = registerCoreHist(pBCTTemplate,s);
			qBCT[s][t].setAxisTitle(X_DIRECTION,"Reconstructed Energy [keV]");
			qBCT[s][t].setAxisTitle(Y_DIRECTION,"<#beta cos #theta>");
			qBCT[s][t].setDrawRange(-1.0,false);
			qBCT[s][t].setDrawRange(1.0,true);
			
			qWrongSide[s][t] = registerCoreHist("hWrongSide_"+itos(t),"Wrong-side ID events",
												nEnergyBins, 0, energyMax, s);

		}
	}
}

void SimAsymmetryAnalyzer::fillCoreHists(ProcessedDataScanner& PDS, double weight) {
	AsymmetryAnalyzer::fillCoreHists(PDS, weight);
	assert(PDS.isSimulated());
	Sim2PMT& S2P = (Sim2PMT&)PDS;
	if(S2P.fType == TYPE_IV_EVENT && S2P.primRadius() < S2P.fiducialRadius) {
		qMissedSpectrum.fillPoint->Fill(S2P.ePrim,weight);
		((TProfile*)qBCT[S2P.costheta<0?EAST:WEST][S2P.fType].fillPoint)->Fill(S2P.ePrim,beta(S2P.ePrim)*S2P.costheta,weight);
	}
	Side s = S2P.fSide;
	if(!(s==EAST || s==WEST)) return;
	if(S2P.fPID != PID_BETA) return;
	if(S2P.passesPositionCut(s) && S2P.fType <= TYPE_III_EVENT)
		((TProfile*)qBCT[s][S2P.fType].fillPoint)->Fill(S2P.getEtrue(),beta(S2P.ePrim)*S2P.costheta,weight);
}

void SimAsymmetryAnalyzer::makePlots() {
	AsymmetryAnalyzer::makePlots();	
	drawQuad(qMissedSpectrum,"Energy/");
	for(EventType t=TYPE_0_EVENT; t<=TYPE_IV_EVENT; ++t)
		drawQuadSides(qBCT[EAST][t],qBCT[WEST][t],false,"BetaCosTheta/");
}

SegmentSaver* SimAsymmetryAnalyzer::makeAnalyzer(const std::string& nm,const std::string& inflname) {
		return isSimulated?new SimAsymmetryAnalyzer(this,nm,inflname):new AsymmetryAnalyzer(this,nm,inflname); 
}

void SimAsymmetryAnalyzer::compareMCtoData(RunAccumulator& OAdata) {
	// re-cast to correct type
	AsymmetryAnalyzer& dat = (AsymmetryAnalyzer&)OAdata;
	
	hAsym->SetLineColor(4);
	dat.hAsym->SetLineColor(2);
	hAsym->Draw("HIST E1");
	dat.hAsym->Draw("SAME HIST E1");
	printCanvas("DataComparison/Asymmetry");
	
	drawHistoPair(dat.hSuperSum,hSuperSum);
	printCanvas("DataComparison/SuperSum");
	for(EventType tp = TYPE_0_EVENT; tp <= TYPE_II_EVENT; ++tp) {
		dat.hEvtSS[tp]->SetMarkerStyle(1);
		hEvtSS[tp]->SetMarkerStyle(1);
		drawHistoPair(dat.hEvtSS[tp],hEvtSS[tp]);
		printCanvas("DataComparison/SuperSum_Type_"+itos(tp));
		
		dat.hTpAsym[tp]->SetMarkerStyle(1);
		hTpAsym[tp]->SetMarkerStyle(1);
		drawHistoPair(dat.hTpAsym[tp],hTpAsym[tp]);
		printCanvas("DataComparison/Asymmetry_Type_"+itos(tp));
	}
	
	for(unsigned int t=TYPE_0_EVENT; t<=TYPE_II_EVENT; t++) {
		std::vector<TH1*> hToPlot;
		for(Side s = EAST; s <= WEST; ++s) {		
			for(AFPState afp = AFP_OFF; afp <= AFP_ON; ++afp) {
				qEnergySpectra[s][nBetaTubes][t].fgbg[afp]->h[GV_OPEN]->SetMarkerColor(2+2*s);
				qEnergySpectra[s][nBetaTubes][t].fgbg[afp]->h[GV_OPEN]->SetMarkerStyle(22+4*afp);
				hToPlot.push_back(qEnergySpectra[s][nBetaTubes][t].fgbg[afp]->h[GV_OPEN]);
				dat.qEnergySpectra[s][nBetaTubes][t].fgbg[afp]->h[GV_OPEN]->SetMarkerColor(2+2*s);
				dat.qEnergySpectra[s][nBetaTubes][t].fgbg[afp]->h[GV_OPEN]->SetMarkerStyle(20+4*afp);
				hToPlot.push_back(dat.qEnergySpectra[s][nBetaTubes][t].fgbg[afp]->h[GV_OPEN]);
			}
		}
		drawSimulHistos(hToPlot,"HIST P");
		printCanvas("DataComparison/Type_"+itos(t));
	}
}

