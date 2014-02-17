#include "StyleSetup.hh"
#include "PathUtils.hh"
#include "CathodeTuningAnalyzer.hh"
#include "OctetSimuCloneManager.hh"

int main(int argc, char *argv[]) {
	
	assert(argc >= 2);
	int octn = atoi(argv[1]);
	ROOTStyleSetup();
	
	OctetSimuCloneManager OSCM("BetaOctetPositions_NewNorm");
	OutputManager OM("ThisNameIsNotUsedAnywhere",getEnvSafe("UCNA_ANA_PLOTS"));
	
	// simulations input setup
	OSCM.simFile="/data2/mmendenhall/G4Out/2010/20120823_neutronBetaUnpol/analyzed_";
	OSCM.simFactor = 1.0;
	OSCM.nTot = 312;
	OSCM.stride = 73;

	//const std::string simOutputDir=OSCM.outputDir+"_Sim0823";
	const std::string simOutputDir="BetaOctetPositions_Sim0823";
	
	WirechamberCalibrator::calibrateCathodes = true;
	ProcessedDataScanner::redoPositions = true;
	
	if(octn < 0) {
	
		SimCathodeTuningAnalyzer CTA_Sim(&OM,simOutputDir);
		if(octn==-1000) OSCM.combineSims(CTA_Sim);
		else OSCM.simOct(CTA_Sim,-octn-1);
	
	} else {
		CathodeTuningAnalyzer CTA(&OM,OSCM.outputDir);
		if(octn>=1000) {
			ROOTStyleSetup(false);
			if(octn==1000) OSCM.combineOcts(CTA);
			else if(octn==1002) {
				// load data and simulation to compare hit distribution shape around each cathode
				CathodeTuningAnalyzer CTAdat(&OM,OSCM.outputDir,RunAccumulator::processedLocation);
				SimCathodeTuningAnalyzer CTsim(&OM,"NameUnused",getEnvSafe("UCNA_ANA_PLOTS")+"/"+simOutputDir+"/"+simOutputDir);
				processCathTweak(*CTAdat.myCT,*CTsim.myCT);
			}
		} else OSCM.scanOct(CTA, octn);
	}

	return 0;
}
