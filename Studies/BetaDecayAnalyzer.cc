#include "BetaDecayAnalyzer.hh"

std::string BetaDecayAnalyzer::processedLocation = "";

BetaDecayAnalyzer::BetaDecayAnalyzer(OutputManager* pnt, const std::string& nm, const std::string& inflName):
OctetAnalyzer(pnt,nm,inflName) {
	ignoreMissingHistos = true;
	addPlugin(myMuons = new MuonAnalyzer(this));
	addPlugin(myPos = new PositionAnalyzer(this));
	addPlugin(myWG = new WirechamberGainAnalyzer(this));
	addPlugin(myAsym = new AsymmetryAnalyzer(this));
	ignoreMissingHistos = false;
}

//-----------------------------------------

SimBetaDecayAnalyzer::SimBetaDecayAnalyzer(OutputManager* pnt, const std::string& nm, const std::string& inflName):
BetaDecayAnalyzer(pnt,nm,inflName) {
	addPlugin(mySimAsym = new SimAsymmetryAnalyzer(this));
}
