#include "BGDecayAnalyzer.hh"

BGDecayAnalyzer::BGDecayAnalyzer(OctetAnalyzer* OA): OctetAnalyzerPlugin(OA,"bgdecay") {
	TH2F hBGDecayTemplate("hBGDecay","energy vs time",60,0,300,80,0,2000);
	for(Side s = EAST; s <= WEST; ++s) {
		qBGDecay[s] = registerCoreHist(hBGDecayTemplate,s);
		qBGDecay[s]->setSubtraction(false);
		qBGDecay[s]->setTimeScaling(false);
		qBGDecay[s]->setAxisTitle(X_DIRECTION,"time [s]");
		qBGDecay[s]->setAxisTitle(Y_DIRECTION,"energy [keV]");
	}
}

void BGDecayAnalyzer::fillCoreHists(ProcessedDataScanner& PDS, double weight) {
	Side s = PDS.fSide;
	if(!(s==EAST || s==WEST)) return;
	if(PDS.fPID == PID_BETA && PDS.fType <= TYPE_III_EVENT && PDS.passesPositionCut(s))
		((TH2F*)qBGDecay[s]->fillPoint)->Fill(PDS.runClock[s],PDS.getEtrue(),weight);
}
