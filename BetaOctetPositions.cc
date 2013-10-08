#include "StyleSetup.hh"
#include "PathUtils.hh"
#include "CathodeTweakAnalyzer.hh"

int main(int argc, char *argv[]) {
	
	assert(argc >= 2);
	int octn = atoi(argv[1]);
	
	ROOTStyleSetup();
	const std::string outputDir="BetaOctetPositions";
	bool doPlots = false;
	RunAccumulator::processedLocation = getEnvSafe("UCNA_ANA_PLOTS")+"/"+outputDir+"/"+outputDir;
	
	if(octn==1000) {
		OutputManager OM("ThisNameIsNotUsedAnywhere",getEnvSafe("UCNA_ANA_PLOTS"));
		MWPCTuningAnalyzer AA(&OM,outputDir);
		processOctets(AA,Octet::loadOctets(QFile(getEnvSafe("UCNA_OCTET_LIST"))),365*24*3600);
	} else {
		Octet oct = Octet::loadOctet(QFile(getEnvSafe("UCNA_OCTET_LIST")),octn);
		if(!oct.getNRuns()) return -1;
		OutputManager OM("ThisNameIsNotUsedAnywhere",getEnvSafe("UCNA_ANA_PLOTS")+"/"+outputDir);
		MWPCTuningAnalyzer AA(&OM,oct.octName());
		processOctets(AA,oct.getSubdivs(oct.divlevel+1,false),0*24*3600, doPlots);
	}
	
	return 0;
}
