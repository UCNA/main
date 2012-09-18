#include "HighEnergyExcess.hh"
#include "GraphUtils.hh"

void fitHighEnergyExcess(QFile& qOut, quadHists* qh, double e0, double e1) {
	for(AFPState afp = AFP_OFF; afp <= AFP_ON; ++afp) {
		TH1* hEn = qh->fgbg[afp]->h[GV_OPEN];
		TH1* hEnBG = qh->fgbg[afp]->h[GV_CLOSED];
		
		int b0 = hEnBG->FindBin(e0);
		int b1 = hEnBG->FindBin(e1)-1;
		if(b1 == hEnBG->GetNbinsX()) ++b1;
		double nBG = hEnBG->Integral(b0,b1);
		double xs = hEn->Integral(b0,b1);
		Double_t d_xs;
		double d_nBG = hEnBG->IntegralAndError(b0,b1,d_xs);
		
		Stringmap m;
		m.insert("side",sideSubst("%c",qh->mySide));
		m.insert("afp",afpWords(afp));
		m.insert("name",qh->name);
		m.insert("nBG",nBG);		// number of BG counts
		m.insert("d_nBG",d_nBG);	// error on BG counts
		m.insert("xs",xs);			// number of excess counts
		m.insert("d_xs",d_xs);		// error on excess counts
		m.insert("eMin",e0);		// lower energy cut
		m.insert("eMax",e1);		// upper energy cut
		m.insert("b0",itos(b0));
		m.insert("b1",itos(b1));
		qOut.insert("bg_subtr_xs",m);
	}
}

HighEnergyExcessAnalyzer::HighEnergyExcessAnalyzer(OctetAnalyzer* OA): OctetAnalyzerPlugin(OA,"highenergy") {
	for(Side s = EAST; s <= WEST; ++s) {
		qExcessr2[s] = registerCoreHist("Excessr2","Excess events radial distribution",10,0,80*80,s);
		qExcessr2[s]->setAxisTitle(X_DIRECTION,"event radius squared [mm^2]");
		qExcessTheta[s] = registerCoreHist("ExcessTheta","Excess events angular distribution",10,-M_PI,M_PI,s);
		qExcessTheta[s]->setAxisTitle(X_DIRECTION,"event position angle [radians]");
		qExcessSpectra[s] = registerCoreHist("ExcessE","Excess high energy event spectrum",72,800,8000,s);
		qExcessSpectra[s]->setAxisTitle(X_DIRECTION,"Energy [keV]");
		qExcessGamma[s] = registerCoreHist("ExcessGamma","Excess high energy gamma event spectrum",72,800,8000,s);
		qExcessGamma[s]->setAxisTitle(X_DIRECTION,"Energy [keV]");
	}
}

void HighEnergyExcessAnalyzer::fillCoreHists(ProcessedDataScanner& PDS, double weight) {
	Side s = PDS.fSide;
	if(!(s==EAST || s==WEST)) return;
	if(PDS.fPID == PID_SINGLE && PDS.fType == TYPE_IV_EVENT)
		qExcessGamma[s]->fillPoint->Fill(PDS.getEtrue(),weight);
	if(PDS.fPID != PID_BETA) return;
	if(PDS.fType <= TYPE_I_EVENT)
		qExcessSpectra[s]->fillPoint->Fill(PDS.getEtrue(),weight);
	if(PDS.fType <= TYPE_I_EVENT && PDS.getEtrue()>1000) {
		qExcessr2[s]->fillPoint->Fill(PDS.radius2(s),weight);
		qExcessTheta[s]->fillPoint->Fill(atan2(PDS.wires[s][Y_DIRECTION].center,PDS.wires[s][X_DIRECTION].center),weight);
	}
}

void HighEnergyExcessAnalyzer::calculateResults() {
	for(Side s = EAST; s <= WEST; ++s) {
		fitHighEnergyExcess(myA->qOut,qExcessSpectra[s],1000,2200);
		fitHighEnergyExcess(myA->qOut,qExcessSpectra[s],2200,7000);
		fitHighEnergyExcess(myA->qOut,qExcessGamma[s],200,1000);
		fitHighEnergyExcess(myA->qOut,qExcessGamma[s],1000,2200);
		fitHighEnergyExcess(myA->qOut,qExcessGamma[s],2200,7000);
	}
}

void HighEnergyExcessAnalyzer::makePlots() {
	if(myA->depth <= 0) {
		drawQuadSides(qExcessGamma[EAST],qExcessGamma[WEST],true,"Energy");
		drawQuadSides(qExcessSpectra[EAST],qExcessSpectra[WEST],true,"Energy");
		drawQuadSides(qExcessr2[EAST],qExcessr2[WEST],true,"Positions");
		drawQuadSides(qExcessTheta[EAST],qExcessTheta[WEST],true,"Positions");
	}
}
