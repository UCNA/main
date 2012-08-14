#include "AnodePositionAnalyzer.hh"

AnodePositionAnalyzer::AnodePositionAnalyzer(RunAccumulator* RA, unsigned int nr): PositionBinnedAnalyzer(RA,"AnodePos",nr) {
	for(Side s = EAST; s <= WEST; ++s) {
		TH1F hTemplate("hAnode_","Anode Signal",200,0,4000);
		sectAnode[s] = allocateSegmentHistograms(hTemplate,AFP_OTHER,s);
	}
}

void AnodePositionAnalyzer::fillCoreHists(ProcessedDataScanner& PDS, double weight) {
	Side s = PDS.fSide;
	if(!(PDS.fType == TYPE_0_EVENT && PDS.fPID == PID_BETA)) return;
	unsigned int m = sects.sector(PDS.wires[s][X_DIRECTION].center,PDS.wires[s][Y_DIRECTION].center);
	if(m>=sects.nSectors()) return;
	sectAnode[s].back()->h[currentGV]->Fill(PDS.mwpcs[s].anode,weight);
	sectAnode[s][m]->h[currentGV]->Fill(PDS.mwpcs[s].anode,weight);
}

void AnodePositionAnalyzer::calculateResults() {
	//TODO
}
