#include "PositionAnalyzer.hh"

PositionAnalyzer::PositionAnalyzer(OctetAnalyzer* OA): OctetAnalyzerPlugin(OA,"position"), offSects(5.0,45.0) {
	myA->ignoreMissingHistos = true;
	for(unsigned int m=0; m<offSects.nSectors(); m++) {
		for(AxisDirection d = X_DIRECTION; d <= Y_DIRECTION; ++d) {
			poff[d].push_back(registerFGBGPair("pOff_"+itos(m), "Type I Position Offsets",50,-25, 25));
			poff[d].back()->setAxisTitle(X_DIRECTION,"E-W offset [mm]");
		}
	}
	myA->ignoreMissingHistos = false;
	
	for(Side s = EAST; s <= WEST; ++s) {
		for(EventType t=TYPE_0_EVENT; t<=TYPE_III_EVENT; ++t) {
			TH2F hPositionsTemplate(("hPos_Type_"+itos(t)).c_str(),
									("Type "+itos(t)+" Positions").c_str(),
									200,-65,65,200,-65,65);
			qPositions[s][t] = registerCoreHist(hPositionsTemplate,s);
			qPositions[s][t]->setAxisTitle(X_DIRECTION,"x Position [mm]");
			qPositions[s][t]->setAxisTitle(Y_DIRECTION,"y Position [mm]");
		}
	}
}

void PositionAnalyzer::fillCoreHists(ProcessedDataScanner& PDS, double weight) {
	Side s = PDS.fSide;
	if(PDS.fPID!=PID_BETA || !(s==EAST||s==WEST)) return;
	if(PDS.fType <= TYPE_III_EVENT)
		((TH2F*)qPositions[s][PDS.fType]->fillPoint)->Fill(PDS.wires[s][X_DIRECTION].center,PDS.wires[s][Y_DIRECTION].center,weight);
	if(PDS.fType != TYPE_I_EVENT) return;
	double x = 0.5*(PDS.wires[EAST][X_DIRECTION].center+PDS.wires[WEST][X_DIRECTION].center);
	double y = 0.5*(PDS.wires[EAST][Y_DIRECTION].center+PDS.wires[WEST][Y_DIRECTION].center);
	unsigned int m = offSects.sector(x,y);
	if(m>=offSects.nSectors()) return;
	for(AxisDirection d = X_DIRECTION; d <= Y_DIRECTION; ++d)
		poff[d][m]->h[currentGV]->Fill(PDS.wires[EAST][d].center-PDS.wires[WEST][d].center,weight);
}

void PositionAnalyzer::calculateResults() {
	
}

void PositionAnalyzer::makePlots() {
	if(myA->depth > 0) return;
	for(Side s = EAST; s <= WEST; ++s) {
		for(unsigned int t=TYPE_0_EVENT; t<=TYPE_II_EVENT; t++) {
			qPositions[s][t]->setDrawRange(0,false);
			drawQuad(qPositions[s][t],"Positions","COL");
		}
	}
}
