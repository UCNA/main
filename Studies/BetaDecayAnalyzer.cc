#include "BetaDecayAnalyzer.hh"

std::string BetaDecayAnalyzer::processedLocation = "";

BetaDecayAnalyzer::BetaDecayAnalyzer(OutputManager* pnt, const std::string& nm, const std::string& inflName):
OctetAnalyzer(pnt,nm,inflName) {
	addPlugin(myMuons = new MuonAnalyzer(this));
	addPlugin(myPos = new PositionAnalyzer(this));
	addPlugin(myAsym = new AsymmetryAnalyzer(this));
}

//-----------------------------------------

SimBetaDecayAnalyzer::SimBetaDecayAnalyzer(OutputManager* pnt, const std::string& nm, const std::string& inflName):
BetaDecayAnalyzer(pnt,nm,inflName) {
	addPlugin(mySimAsym = new SimAsymmetryAnalyzer(this));
}
