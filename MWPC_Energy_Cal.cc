#include "StyleSetup.hh"
#include "ControlMenu.hh"
#include "PathUtils.hh"
#include "OctetSimuCloneManager.hh"
#include "WirechamberGainMapPlugins.hh"
#include "WirechamberEnergyPlugins.hh"
#include "CathodeTuningAnalyzer.hh"

/// Wirechamber energy calibrations analyzer
class MWPC_Ecal_Analyzer: public RunAccumulator {
public:
	/// constructor
	MWPC_Ecal_Analyzer(unsigned int nr, OutputManager* pnt, const std::string& nm, const std::string& inflname = ""):
	RunAccumulator(pnt,nm,inflname), nRings(nr) {
		ignoreMissingHistos = true;
		addPlugin(anode_plgn = new AnodeGainMapPlugin(this,nRings));
		addPlugin(ccloud_plgn = new CCloudGainMapPlugin(this,nRings));
		addPlugin(nedep_plgn = new WirechamberNullEdepMapPlugin(this,nRings));
		addPlugin(wgain_plgn = new MWPCGainPlugin(this));
		addPlugin(cgain_plgn = new CathodeGainPlugin(this));
		addPlugin(cms_plgn = new WirechamberCathMaxSumPlugin(this));
		ignoreMissingHistos = false;
	}
	
	/// cloning generator: just return another of the same subclass (with any settings you want to preserve)
	virtual SegmentSaver* makeAnalyzer(const std::string& nm,
									   const std::string& inflname) { return new MWPC_Ecal_Analyzer(nRings,this,nm,inflname); }

	const unsigned int nRings;
	AnodeGainMapPlugin* anode_plgn;
	CCloudGainMapPlugin* ccloud_plgn;
	MWPCGainPlugin* wgain_plgn;
	WirechamberNullEdepMapPlugin* nedep_plgn;
	CathodeGainPlugin* cgain_plgn;
	WirechamberCathMaxSumPlugin* cms_plgn;
};

/// Simulated wirechamber data energy calibrations analyzer
class Sim_MWPC_Ecal_Analyzer: public RunAccumulator {
public:
	/// constructor
	Sim_MWPC_Ecal_Analyzer(unsigned int nr, OutputManager* pnt, const std::string& nm, const std::string& inflname = ""):
	RunAccumulator(pnt,nm,inflname), nRings(nr) {
		ignoreMissingHistos = true;
		addPlugin(edep_plgn = new WirechamberEdepMapPlugin(this,nRings));
		addPlugin(wgain_plgn = new MWPCGainPlugin(this));
		addPlugin(sim23_plgn = new WirechamberSimBackscattersPlugin(this));
		addPlugin(cgain_plgn = new CathodeGainPlugin(this));
		addPlugin(cms_plgn = new WirechamberCathMaxSumPlugin(this));
		addPlugin(ste_plgn = new WirechamberSimTrigEfficPlugin(this));
		ignoreMissingHistos = false;
	}
	/// cloning generator: just return another of the same subclass (with any settings you want to preserve)
	virtual SegmentSaver* makeAnalyzer(const std::string& nm,
									   const std::string& inflname) { return new Sim_MWPC_Ecal_Analyzer(nRings,this, nm, inflname); }
	const unsigned int nRings;
	WirechamberEdepMapPlugin* edep_plgn;
	MWPCGainPlugin* wgain_plgn;
	WirechamberSimBackscattersPlugin* sim23_plgn;
	CathodeGainPlugin* cgain_plgn;
	WirechamberCathMaxSumPlugin* cms_plgn;
	WirechamberSimTrigEfficPlugin* ste_plgn;
};

/// run MWPC analyzers over beta decay data or simulation octets
void mi_MWPCCal(std::deque<std::string>&, std::stack<std::string>& stack) {

	unsigned int nRings = streamInteractor::popInt(stack);
	std::string tpSelect = streamInteractor::popString(stack);
	int octn = streamInteractor::popInt(stack);
	
	OctetSimuCloneManager OSCM("MWPC_ECal_"+itos(nRings));
	OutputManager OM("ThisNameIsNotUsedAnywhere",getEnvSafe("UCNA_ANA_PLOTS"));
	
	// simulations input setup
	OSCM.simFile="/data2/mmendenhall/G4Out/2010/20120823_neutronBetaUnpol/analyzed_";
	OSCM.simFactor = 1.0;
	OSCM.nTot = 312;
	OSCM.stride = 73;
	const std::string simOutputDir="MWPC_ECal_"+itos(nRings)+"_Sim0823";
	
	if(tpSelect=="sim") {
		Sim_MWPC_Ecal_Analyzer MEA_Sim(nRings,&OM,simOutputDir);
		if(octn==1000) {
			OSCM.combineSims(MEA_Sim);
			MEA_Sim.edep_plgn->genPosmap("mwpc_sim");
		} else if(octn==1001) {
			MWPC_Ecal_Analyzer BDA(nRings,&OM,OSCM.outputDir,RunAccumulator::processedLocation);
			Sim_MWPC_Ecal_Analyzer BDA_MC(nRings,&OM,simOutputDir,OSCM.baseDir+"/"+simOutputDir+"/"+simOutputDir);
			BDA_MC.compareMCtoData(BDA);
			BDA_MC.write();
		} else {
			OSCM.hoursOld = 24*30;
			OSCM.doCompare = true;
			OSCM.simOct(MEA_Sim,octn);
		}
	} else {
		MWPC_Ecal_Analyzer MEA(nRings,&OM,OSCM.outputDir);
		if(octn==1000) {
			OSCM.combineOcts(MEA);
			MEA.anode_plgn->genPosmap("anode");
			MEA.ccloud_plgn->genPosmap("ccloud");
		} else {
			//OSCM.doPlots = true;
			OSCM.scanOct(MEA, octn);
		}
	}
}



//-------------------------------------------------------------------

int main(int argc, char *argv[]) {

	ROOTStyleSetup(false);
	
	NameSelector selectDatSim("Data/Simulation");
	selectDatSim.addChoice("Beta decay data octet","dat");
	selectDatSim.addChoice("Simulated octet","sim");
	selectDatSim.setDefault("dat");
				
	inputRequester exitMenu("Exit Menu",&menutils_Exit);
	
	inputRequester pm_MWPC_Oct("Process Beta Octets",&mi_MWPCCal);
	pm_MWPC_Oct.addArg("Octet number","","Beta decay octet in list, or 1000 to combine previously processed octets.");
	pm_MWPC_Oct.addArg(&selectDatSim);
	pm_MWPC_Oct.addArg("n rings","8","Number of rings for subdividing fiducial volume");
	
	// main menu
	OptionsMenu OM("Wirechamber Calibrations Menu");
	OM.addChoice(&pm_MWPC_Oct,"oct");
	OM.addChoice(&exitMenu,"x");
	
	// load command line arguments
	std::deque<std::string> args;
	for(int i=1; i<argc; i++)
		args.push_back(argv[i]);
	std::stack<std::string> stack;
	OM.doIt(args,stack);
	
	printf("\n\n\n>>>>> Goodbye. <<<<<\n\n\n");
	return 0;
}
