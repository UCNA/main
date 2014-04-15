#include <cassert>
#include <fstream>

#include "OutputManager.hh"
#include "RunAccumulator.hh"
#include "G4toPMT.hh"
#include "PositionsPlugin.hh"
#include "PathUtils.hh"
#include "BetaDecayAnalyzer.hh"

class PosOffsetAnalyzer: public RunAccumulator {
public:
	PosOffsetAnalyzer(OutputManager* pnt, const std::string& nm, const std::string& inflName = ""): RunAccumulator(pnt,nm,inflName) {
		addPlugin(new PositionOffsetsPlugin(this));
	}
	virtual SegmentSaver* makeAnalyzer(const std::string& nm, const std::string& inflname) { return new PosOffsetAnalyzer(this,nm,inflname); }
};


int main(int argc, char *argv[]) {
	if(argc != 3) {
		assert(argc);
		printf("Use: %s [analyzer name] [analysis directory]\n", argv[0]);
		return 0;
	}
	
	// set up output paths
	auto v = split(argv[2],"/");
	std::string nm = v.back();
	v.pop_back();
	std::string pth = "/"+join(v,"/");
	OutputManager OM("NameUnused",pth);
	std::string slist = std::string(argv[2])+"/simlist.txt";
	
	// set up analyzer
	RunAccumulator* RA = NULL;
	std::string anm(argv[1]);
	if(anm=="beta") RA = new SimBetaDecayAnalyzer(&OM,nm);
	else if(anm=="posoff") RA = new PosOffsetAnalyzer(&OM,nm);
	
	if(!RA) {
		printf("Unknown analyzer selection '%s'!\n", argv[1]);
		return -1;
	}

	if(!fileExists(slist)) {
		// combine prior simulations
		RA->mergeDir();
		RA->makeOutput(true);
	} else {
		// set up input simulations
		G4toPMT g2p;
		ifstream file;
		file.open(slist);
		std::string fname;
		while(!file.fail() && !file.eof()){
			fname = "";
			file >> fname;
			if(fname.size())
				g2p.addFile(fname);
		}
		file.close();
		PMTCalibrator PCal(16000);
		g2p.setCalibrator(PCal);

		// process data and exit
		RA->loadSimData(g2p);
		RA->makeOutput(false);
	}
	delete RA;
	return 0;
}
