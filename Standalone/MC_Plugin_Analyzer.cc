#include <cassert>
#include <fstream>

#include "OutputManager.hh"
#include "RunAccumulator.hh"
#include "G4toPMT.hh"
#include "PositionsPlugin.hh"
#include "PathUtils.hh"
#include "BetaDecayAnalyzer.hh"
#include "GraphicsUtils.hh"

#include <TGraphAsymmErrors.h>

class PosOffsetAnalyzer: public RunAccumulator {
public:
	PosOffsetAnalyzer(OutputManager* pnt, const std::string& nm, const std::string& inflName = ""): RunAccumulator(pnt,nm,inflName) {
		//addPlugin(new PositionOffsetsPlugin(this)); // TODO split this back out (restore lost changes)
	}
	virtual SegmentSaver* makeAnalyzer(const std::string& nm, const std::string& inflname) { return new PosOffsetAnalyzer(this,nm,inflname); }
};

/// plugin with extra plots for simulated runs
class MC_Effic_Plugin: public AnalyzerPlugin {
public:
	/// constructor
	MC_Effic_Plugin(RunAccumulator* RA): AnalyzerPlugin(RA,"MC_Effic"), gEffic(NULL)  {
		for(unsigned int i=0; i<2; i++) {
			energyEffic[i] = registerFGBGPair("hEnergyEffic_"+itos(i), "Event Counts", 80, 0, 800, AFP_OTHER, BOTH);
			energyEffic[i]->setAxisTitle(X_DIRECTION,"energy [keV]");
			energyEffic[i]->setAxisTitle(Y_DIRECTION,"number of events");
		}
	}
	
	/// fill from scan data point
	virtual void fillCoreHists(ProcessedDataScanner& PDS, double weight) {
		Sim2PMT& S2P = dynamic_cast<Sim2PMT&>(PDS);
		if(S2P.primRadius() > 45) return;
		energyEffic[false]->h[currentGV]->Fill(S2P.ePrim, weight);
		if(S2P.fPID == PID_BETA) energyEffic[true]->h[currentGV]->Fill(S2P.ePrim, weight);
	}
	
	/// make combined eloss table
	virtual void calculateResults() {
		gEffic = new TGraphAsymmErrors(energyEffic[true]->h[GV_OPEN]->GetNbinsX());
		gEffic->BayesDivide(energyEffic[true]->h[GV_OPEN], energyEffic[false]->h[GV_OPEN], "w");
	}
	/// plot results
	virtual void makePlots() {
		
		myA->defaultCanvas->SetCanvasSize(300,200);
		myA->defaultCanvas->SetLeftMargin(0.14);
		myA->defaultCanvas->SetRightMargin(0.04);
		
		gEffic->SetMinimum(0.95);
		gEffic->SetMaximum(1.0);
		gEffic->Draw("AP");
		gEffic->SetTitle("Simulated detector efficiency");
		gEffic->GetXaxis()->SetTitle("beta energy [keV]");
		gEffic->GetXaxis()->SetLimits(0,782);
		gEffic->GetYaxis()->SetTitle("detected event fraction");
		gEffic->GetYaxis()->SetTitleOffset(1.8);
		
		gEffic->SetMarkerStyle(33);
		gEffic->SetMarkerSize(0.5);
		gEffic->Draw("AP");
		printCanvas("SimDetEffic");
		
		drawHistoPair(energyEffic[true]->h[GV_OPEN],energyEffic[false]->h[GV_OPEN]);
		printCanvas("SimDetEvts");
	}
	
	fgbgPair* energyEffic[2];		///< detector efficiency as function of primary energy
	TGraphAsymmErrors* gEffic;		///< extracted efficiency curve
};

class MC_Effic_Analyzer: public RunAccumulator {
public:
	MC_Effic_Analyzer(OutputManager* pnt, const std::string& nm, const std::string& inflName = ""): RunAccumulator(pnt,nm,inflName) {
		addPlugin(new MC_Effic_Plugin(this));
	}
	virtual SegmentSaver* makeAnalyzer(const std::string& nm, const std::string& inflname) { return new MC_Effic_Analyzer(this,nm,inflname); }
};



int main(int argc, char *argv[]) {
	if(argc != 4) {
		assert(argc);
		printf("Use: %s [analyzer name] [analysis directory] [run to compare]\n", argv[0]);
		return 0;
	}
	
	// set up output paths
	auto v = split(argv[2],"/");
	std::string nm = v.back();
	v.pop_back();
	std::string pth = "/"+join(v,"/");
	OutputManager OM("NameUnused",pth);
	std::string slist = std::string(argv[2])+"/simlist.txt";
	int run = atoi(argv[3]);
	
	// set up analyzer
	RunAccumulator* RA = NULL;
	std::string anm(argv[1]);
	if(anm=="beta") RA = new SimBetaDecayAnalyzer(&OM,nm);
	else if(anm=="posoff") RA = new PosOffsetAnalyzer(&OM,nm);
	else if(anm=="effic") RA = new  MC_Effic_Analyzer(&OM,nm);
	
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
		PMTCalibrator PCal(run);
		std::cout << "Here I am\n";
		g2p.setCalibrator(PCal);
		std::cout << "Here I am now\n";

		RunInfo RI = CalDBSQL::getCDB()->getRunInfo(run);
		g2p.setAFP(RI.afpState);
		// process data and exit
		RA->loadSimData(g2p);
		
		//RA->makeOutput(true);
	}
	delete RA;
	return 0;
}
