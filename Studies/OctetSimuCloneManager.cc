#include "OctetSimuCloneManager.hh"
#include "PathUtils.hh"

OctetSimuCloneManager::OctetSimuCloneManager(const std::string& dname, const std::string& bdir):
outputDir(dname), baseDir(bdir), doPlots(false), simFactor(1.0), nTot(0), stride(0) {
	RunAccumulator::processedLocation = baseDir+"/"+outputDir+"/"+outputDir;
}

void OctetSimuCloneManager::scanOct(RunAccumulator& RA, unsigned int octn) {
		Octet oct = Octet::loadOctet(QFile(getEnvSafe("UCNA_OCTET_LIST")),octn);
		if(!oct.getNRuns()) return;
		RunAccumulator* octRA = (RunAccumulator*)RA.makeAnalyzer(oct.octName(),"");
		processOctets(*octRA,oct.getSubdivs(oct.divlevel+1,false),0*24*3600, doPlots);
		delete octRA;
}

void OctetSimuCloneManager::combineOcts(RunAccumulator& RA) {
	RA.mergeOcts(Octet::loadOctets(QFile(getEnvSafe("UCNA_OCTET_LIST"))));
}

Sim2PMT* OctetSimuCloneManager::getSimdata(unsigned int octn) {
	G4toPMT* simData = new G4toPMT();
	for(unsigned int i=0; i<stride; i++)
		simData->addFile(simFile+itos((stride*octn+i)%nTot)+".root");
	simData->PGen[EAST].xscatter = simData->PGen[WEST].xscatter = 0.005;
	return simData;
}

void OctetSimuCloneManager::simOct(RunAccumulator& SimRA, unsigned int octn) {
		Octet oct = Octet::loadOctet(QFile(getEnvSafe("UCNA_OCTET_LIST")),octn);
		if(!oct.getNRuns()) {
			printf("No runs in octet %i.\n",octn);
			return;
		}
		RunAccumulator* octSim = (RunAccumulator*)SimRA.makeAnalyzer(oct.octName(),"");
		Sim2PMT* simData = getSimdata(octn);
		octSim->simuClone(getEnvSafe("UCNA_ANA_PLOTS")+"/"+outputDir+"/"+oct.octName(), *simData, simFactor, 0.*3600, doPlots);
		delete simData;
		delete octSim;
}

void OctetSimuCloneManager::combineSims(RunAccumulator& SimRA) {
	SimRA.mergeSims(getEnvSafe("UCNA_ANA_PLOTS")+"/"+outputDir);
}
