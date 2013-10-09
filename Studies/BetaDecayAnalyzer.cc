#include "BetaDecayAnalyzer.hh"

BetaDecayAnalyzer::BetaDecayAnalyzer(OutputManager* pnt, const std::string& nm, const std::string& inflName):
OctetAnalyzer(pnt,nm,inflName) {
	addPlugin(myMuons = new MuonPlugin(this));
	addPlugin(myPos = new PositionsPlugin(this));
	addPlugin(myAG = new AnodeGainPlugin(this));
	addPlugin(myCG = new CathodeGainPlugin(this));
	addPlugin(myAsym = new AsymmetryPlugin(this));
	//addPlugin(myHEE = new HighEnergyExcessPlugin(this));
}

//-----------------------------------------

SimBetaDecayAnalyzer::SimBetaDecayAnalyzer(OutputManager* pnt, const std::string& nm, const std::string& inflName):
OctetAnalyzer(pnt,nm,inflName) {
	addPlugin(myPos = new PositionsPlugin(this));
	addPlugin(myAG = new AnodeGainPlugin(this));
	addPlugin(myAsym = new AsymmetryPlugin(this));

	addPlugin(mySimAsym = new SimAsymmetryPlugin(this));
	addPlugin(mySimID = new WirechamberSimBackscattersPlugin(this));
}

SegmentSaver* SimBetaDecayAnalyzer::makeAnalyzer(const std::string& nm, const std::string& inflname) {
	if(isSimulated)
		return new SimBetaDecayAnalyzer(this,nm,inflname);
	else
		return new BetaDecayAnalyzer(this,nm,inflname);
}
