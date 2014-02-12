#include "BetaDecayAnalyzer.hh"

BetaDecayAnalyzer::BetaDecayAnalyzer(OutputManager* pnt, const std::string& nm, const std::string& inflName):
OctetAnalyzer(pnt,nm,inflName) {
	addPlugin(myMuons = new MuonPlugin(this));
	addPlugin(myPos = new PositionsPlugin(this));
	addPlugin(myAsym = new AsymmetryPlugin(this));
	ignoreMissingHistos = true;
	addPlugin(myHEE = new HighEnergyExcessPlugin(this));
	ignoreMissingHistos = false;
}

//-----------------------------------------

SimBetaDecayAnalyzer::SimBetaDecayAnalyzer(OutputManager* pnt, const std::string& nm, const std::string& inflName):
OctetAnalyzer(pnt,nm,inflName) {
	addPlugin(myPos = new PositionsPlugin(this));
	addPlugin(myAsym = new AsymmetryPlugin(this));
	addPlugin(mySimAsym = new SimAsymmetryPlugin(this));
}

SegmentSaver* SimBetaDecayAnalyzer::makeAnalyzer(const std::string& nm, const std::string& inflname) {
	if(isSimulated)
		return new SimBetaDecayAnalyzer(this,nm,inflname);
	else
		return new BetaDecayAnalyzer(this,nm,inflname);
}
