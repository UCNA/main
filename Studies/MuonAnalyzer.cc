#include "MuonAnalyzer.hh"

MuonAnalyzer::MuonAnalyzer(OctetAnalyzer* OA): OctetAnalyzerPlugin(OA,"muon"), nEnergyBins(150), energyMax(1500) {
	for(Side s = EAST; s <= WEST; ++s) {
		qMuonSpectra[s][false] = registerCoreHist("hMuonNoSub", "Tagged Muon Events Energy",
													   nEnergyBins, 0, energyMax, s);
		qMuonSpectra[s][false]->setSubtraction(false);
		qMuonSpectra[s][false]->setAxisTitle(X_DIRECTION,"Energy [keV]");
		qMuonSpectra[s][true] = registerCoreHist("hMuonSpectrum", "Tagged Muon Events Energy",
													  nEnergyBins, 0, energyMax, s);
		qMuonSpectra[s][true]->setAxisTitle(X_DIRECTION,"Energy [keV]");
		
		qBackMuons[s][false] = registerCoreHist("hBackMuNoSub", "Tagged Muon Events Energy",
													 nEnergyBins, 0, energyMax, s);
		qBackMuons[s][false]->setSubtraction(false);
		qBackMuons[s][false]->setAxisTitle(X_DIRECTION,"Energy [keV]");
		qBackMuons[s][true] = registerCoreHist("hBackMu", "Tagged Muon Events Energy",
													nEnergyBins, 0, energyMax, s);
		qBackMuons[s][true]->setAxisTitle(X_DIRECTION,"Energy [keV]");
		
		TH2F hMuPosTemplate("hMuonPos","Muon event positions",100,-65,65,100,-65,65);
		TH2F hBackMuPosTemplate("hBackMuPos","Backing veto muon positions",100,-65,65,100,-65,65);
		pMuonPos[s] = registerFGBGPair(hMuPosTemplate,AFP_OTHER,s);
		pMuonPos[s]->doSubtraction = false;
		pBackMuPos[s] = registerFGBGPair(hBackMuPosTemplate,AFP_OTHER,s);
		pBackMuPos[s]->doSubtraction = false;
		pMuonPos[s]->setAxisTitle(X_DIRECTION,"x Position [mm]");
		pMuonPos[s]->setAxisTitle(Y_DIRECTION,"y Position [mm]");
		pBackMuPos[s]->setAxisTitle(X_DIRECTION,"x Position [mm]");
		pBackMuPos[s]->setAxisTitle(Y_DIRECTION,"y Position [mm]");
	}
}

void MuonAnalyzer::fillCoreHists(ProcessedDataScanner& PDS, double weight) {
	Side s = PDS.fSide;
	if(!(s==EAST || s==WEST)) return;
	if(!(PDS.fPID == PID_MUON && PDS.fType <= TYPE_III_EVENT)) return;
	if(PDS.passesPositionCut(s)) {
		qMuonSpectra[s][false]->fillPoint->Fill(PDS.getEtrue(),weight);
		qMuonSpectra[s][true]->fillPoint->Fill(PDS.getEtrue(),weight);
	}
	((TH2F*)pMuonPos[s]->h[currentGV])->Fill(PDS.wires[s][X_DIRECTION].center,PDS.wires[s][Y_DIRECTION].center,weight);
	
	if(!PDS.fTaggedBack[s]) return;
	if(PDS.passesPositionCut(s)) {
		qBackMuons[s][false]->fillPoint->Fill(PDS.getEtrue(),weight);
		qBackMuons[s][true]->fillPoint->Fill(PDS.getEtrue(),weight);
	}
	((TH2F*)pBackMuPos[s]->h[currentGV])->Fill(PDS.wires[s][X_DIRECTION].center,PDS.wires[s][Y_DIRECTION].center,weight);
}

void MuonAnalyzer::calculateResults() {
	TF1 fLandau("fLandau","landau",500,1500);
	for(Side s = EAST; s <= WEST; ++s) {
		fLandau.SetLineColor(2+2*s);
		for(AFPState afp = AFP_OFF; afp <= AFP_ON; ++afp) {
			qBackMuons[s][false]->fgbg[afp]->h[GV_OPEN]->Fit(&fLandau,"QR+");
			Stringmap m;
			m.insert("side",sideSubst("%c",s));
			m.insert("afp",afpWords(afp));
			m.insert("height",fLandau.GetParameter(0));
			m.insert("d_height",fLandau.GetParError(0));
			m.insert("mpv",fLandau.GetParameter(1));
			m.insert("d_mpv",fLandau.GetParError(1));
			m.insert("sigma",fLandau.GetParameter(2));
			m.insert("d_sigma",fLandau.GetParError(2));
			myA->qOut.insert("muback_fit",m);
		}
	}
}

void MuonAnalyzer::makePlots() {
	drawQuadSides(qMuonSpectra[EAST][true], qMuonSpectra[WEST][true], true, "MuonSpectra");
	drawQuadSides(qMuonSpectra[EAST][false], qMuonSpectra[WEST][false], true, "MuonSpectra");
	drawQuadSides(qBackMuons[EAST][true], qBackMuons[WEST][true], true, "MuonSpectra");
	drawQuadSides(qBackMuons[EAST][false], qBackMuons[WEST][false], true, "MuonSpectra");
	
	// positions
	if(myA->depth <= 0) {
		for(Side s = EAST; s <= WEST; ++s) {
			pMuonPos[s]->h[GV_OPEN]->Draw("COL");
			printCanvas(sideSubst("MuonSpectra/MuPos_%c",s));
			pBackMuPos[s]->h[GV_OPEN]->Draw("COL");
			printCanvas(sideSubst("MuonSpectra/BackMuPos_%c",s));
		}
	}
}
