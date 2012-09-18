#include "BetaDecayAnalyzer.hh"

std::string BetaDecayAnalyzer::processedLocation = "";

BetaDecayAnalyzer::BetaDecayAnalyzer(OutputManager* pnt, const std::string& nm, const std::string& inflName):
OctetAnalyzer(pnt,nm,inflName) {
	addPlugin(myMuons = new MuonAnalyzer(this));
	addPlugin(myPos = new PositionAnalyzer(this));
	addPlugin(myAG = new AnodeGainAnalyzer(this));
	addPlugin(myCG = new CathodeGainAnalyzer(this));
	addPlugin(myAsym = new AsymmetryAnalyzer(this));
	//addPlugin(myHEE = new HighEnergyExcessAnalyzer(this));
}

//-----------------------------------------

SimBetaDecayAnalyzer::SimBetaDecayAnalyzer(OutputManager* pnt, const std::string& nm, const std::string& inflName):
OctetAnalyzer(pnt,nm,inflName) {
	addPlugin(myPos = new PositionAnalyzer(this));
	addPlugin(myAG = new AnodeGainAnalyzer(this));
	addPlugin(myAsym = new AsymmetryAnalyzer(this));

	addPlugin(mySimAsym = new SimAsymmetryAnalyzer(this));
	addPlugin(mySimID = new WirechamberSimTypeID(this));
}

SegmentSaver* SimBetaDecayAnalyzer::makeAnalyzer(const std::string& nm, const std::string& inflname) {
	if(isSimulated)
		return new SimBetaDecayAnalyzer(this,nm,inflname);
	else
		return new BetaDecayAnalyzer(this,nm,inflname);
}
