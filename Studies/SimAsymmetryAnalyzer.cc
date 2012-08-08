#include "SimAsymmetryAnalyzer.hh"
#include "BetaSpectrum.hh"
#include "GraphicsUtils.hh"
#include <TProfile.h>

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
			
			qWrongSide[s][t] = registerCoreHist("hWrongSide_"+itos(t),"Wrong-side ID events",
													 nEnergyBins, 0, energyMax, s);
			qWrongSide[s][t]->setAxisTitle(X_DIRECTION,"Reconstructed Energy [keV]");
		}
	}
}

void SimAsymmetryAnalyzer::fillCoreHists(ProcessedDataScanner& PDS, double weight) {
	assert(PDS.isSimulated());
	Sim2PMT& S2P = (Sim2PMT&)PDS;
	if(S2P.fType == TYPE_IV_EVENT && S2P.primRadius() < S2P.fiducialRadius) {
		qMissedSpectrum->fillPoint->Fill(S2P.ePrim,weight);
		((TProfile*)qBCT[S2P.costheta<0?EAST:WEST][S2P.fType]->fillPoint)->Fill(S2P.ePrim,beta(S2P.ePrim)*S2P.costheta,weight);
	}
	Side s = S2P.fSide;
	if(!(s==EAST || s==WEST)) return;
	if(S2P.fPID != PID_BETA) return;
	if(S2P.passesPositionCut(s) && S2P.fType <= TYPE_III_EVENT) {
		((TProfile*)qBCT[s][S2P.fType]->fillPoint)->Fill(S2P.getEtrue(),beta(S2P.ePrim)*S2P.costheta,weight);
		if( (s==EAST && S2P.costheta > 0) || (s==WEST && S2P.costheta < 0) )
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

