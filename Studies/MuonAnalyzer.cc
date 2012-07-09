#include "MuonAnalyzer.hh"

std::string MuonAnalyzer::processedLocation = "";	// set this later depending on situtation

MuonAnalyzer::MuonAnalyzer(OutputManager* pnt, const std::string& nm, const std::string& inflname): OctetAnalyzer(pnt,nm,inflname) {
	ignoreMissingHistos = true;
	for(Side s = EAST; s <= WEST; ++s) {
		qMuonSpectra[s][false] = registerCoreHist("hMuonNoSub", "Tagged Muon Events Energy",150, 0, 1500, s, &hMuonSpectra[s][false]);
		qMuonSpectra[s][false].setSubtraction(false);
		qMuonSpectra[s][true] = registerCoreHist("hMuonSpectrum", "Tagged Muon Events Energy",150, 0, 1500, s, &hMuonSpectra[s][true]);
		
		qBackMuons[s][false] = registerCoreHist("hBackMuNoSub", "Tagged Muon Events Energy",150, 0, 1500, s, &hBackMuons[s][false]);
		qBackMuons[s][false].setSubtraction(false);
		qBackMuons[s][true] = registerCoreHist("hBackMu", "Tagged Muon Events Energy",150, 0, 1500, s, &hBackMuons[s][true]);
		
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
	ignoreMissingHistos = false;
}

void MuonAnalyzer::fillCoreHists(ProcessedDataScanner& PDS, double weight) {
	Side s = PDS.fSide;
	if(!(s==EAST || s==WEST)) return;
	if(!(PDS.fPID == PID_MUON && PDS.fType <= TYPE_III_EVENT)) return;
	if(PDS.passesPositionCut(s)) {
		hMuonSpectra[s][false]->Fill(PDS.getEtrue(),weight); 
		hMuonSpectra[s][true]->Fill(PDS.getEtrue(),weight); 
	}
	((TH2F*)pMuonPos[s]->h[currentGV])->Fill(PDS.wires[s][X_DIRECTION].center,PDS.wires[s][Y_DIRECTION].center,weight);
	
	if(!PDS.fTaggedBack[s]) return;
	if(PDS.passesPositionCut(s)) {
		hBackMuons[s][false]->Fill(PDS.getEtrue(),weight); 
		hBackMuons[s][true]->Fill(PDS.getEtrue(),weight); 
	}
	((TH2F*)pBackMuPos[s]->h[currentGV])->Fill(PDS.wires[s][X_DIRECTION].center,PDS.wires[s][Y_DIRECTION].center,weight);
}

void MuonAnalyzer::makePlots() {
	drawQuadSides(qMuonSpectra[EAST][true], qMuonSpectra[WEST][true], true, "MuonSpectra");
	drawQuadSides(qMuonSpectra[EAST][false], qMuonSpectra[WEST][false], true, "MuonSpectra");
	drawQuadSides(qBackMuons[EAST][true], qBackMuons[WEST][true], true, "MuonSpectra");
	drawQuadSides(qBackMuons[EAST][false], qBackMuons[WEST][false], true, "MuonSpectra");
	
	// positions
	if(depth <= 0) {
		for(Side s = EAST; s <= WEST; ++s) {
			pMuonPos[s]->h[GV_OPEN]->Draw("COL");
			printCanvas(sideSubst("MuonSpectra/MuPos_%c",s));
			pBackMuPos[s]->h[GV_OPEN]->Draw("COL");
			printCanvas(sideSubst("MuonSpectra/BackMuPos_%c",s));
		}
	}
}
