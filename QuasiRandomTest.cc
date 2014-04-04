#include "G4toPMT.hh"
#include "PMTGenerator.hh"
#include "OutputManager.hh"
#include "OctetAnalyzer.hh"
#include "AsymmetryPlugin.hh"
#include "PathUtils.hh"
#include "StyleSetup.hh"
#include "AnalysisDB.hh"

/// recursively simulate binary tree of data points for convergence tests
void SimTree(OctetAnalyzer& SimRA, Sim2PMT& SP, unsigned int nlevels, unsigned int plotlevel, unsigned int npts, unsigned int indx = 0) {
	
	AnalysisDB::disableADB = true;
	SimRA.simPerfectAsym = true;
	
	printf("\n------ Generating data in '%s'...\n",SimRA.basePath.c_str());
	
	if(nlevels>0) {
		// generate sub-directories
		for(unsigned int i=0; i<2; i++) {
			OctetAnalyzer* subRA = (OctetAnalyzer*)SimRA.makeAnalyzer(itos(nlevels-1)+"_"+itos(2*indx+i),"");
			SimTree(*subRA,SP,nlevels-1,plotlevel,npts,2*indx+i);
			SimRA.addSegment(*subRA);
			delete(subRA);
		}
	} else {
		// generate simulated data
		SimRA.loadBothAFP(SP, npts);
		for(AFPState afp = AFP_OFF; afp <= AFP_ON; ++afp)
			SimRA.totalTime[afp][GV_OPEN] += 1.0;
	}
	
	// generate output
	SimRA.makeOutput(nlevels>=plotlevel);
	SimRA.write();
}


/// analyzer for beta decay data
class AsymAnalyzer: public OctetAnalyzer {
public:
	/// constructor
	AsymAnalyzer(OutputManager* pnt, const std::string& nm = "AsymAnalyzer", const std::string& inflName = ""): OctetAnalyzer(pnt,nm,inflName) { addPlugin(myAsym = new AsymmetryPlugin(this)); }
	/// create a new instance of this object (cloning self settings) for given directory
	virtual SegmentSaver* makeAnalyzer(const std::string& nm, const std::string& inflname) { return new AsymAnalyzer(this,nm,inflname); }
	AsymmetryPlugin* myAsym;		//< asymmetry plugin
};


//---------------------------------------------------------------------------------------

int main(int, char**) {

	ROOTStyleSetup();
	
	OutputManager OM("nameUnused",getEnvSafe("UCNA_ANA_PLOTS")+"/test/");
	AsymAnalyzer AA(&OM,"SimStats_PR");
	
	G4toPMT g2p;
	g2p.addFile(getEnvSafe("G4WORKDIR")+"/output/2014013_GeomC_n1_f_p/analyzed_*.root");
	PMTCalibrator PCal(16000);
	g2p.setCalibrator(PCal);
	g2p.startScan();
	
	SimTree(AA, g2p, 9, 7, 1000);

	return 0;
}
