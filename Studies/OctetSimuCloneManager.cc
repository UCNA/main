#include "OctetSimuCloneManager.hh"
#include "PathUtils.hh"

OctetSimuCloneManager::OctetSimuCloneManager(const std::string& dname, const std::string& bdir):
outputDir(dname), baseDir(bdir), doPlots(false), doCompare(false), hoursOld(0), simFactor(1.0), nTot(0), stride(0), ownSimData(false), simData(NULL) {
	RunAccumulator::processedLocation = baseDir+"/"+outputDir+"/"+outputDir;
}

void OctetSimuCloneManager::scanOct(RunAccumulator& RA, const Octet& oct) {
	if(!oct.getNRuns()) return;
	RunAccumulator* octRA = (RunAccumulator*)RA.makeAnalyzer(oct.octName(),"");
	processOctets(*octRA,oct.getSubdivs(nextDiv(oct.divlevel),false),hoursOld*3600,doPlots);
	delete octRA;
}

void OctetSimuCloneManager::scanOct(RunAccumulator& RA, unsigned int octn) {
		Octet oct = Octet::loadOctet(QFile(getEnvSafe("UCNA_OCTET_LIST")),octn);
		scanOct(RA,oct);
}

void OctetSimuCloneManager::combineOcts(RunAccumulator& RA) {
	RA.mergeOcts(Octet::loadOctets(QFile(getEnvSafe("UCNA_OCTET_LIST"))));
}

unsigned int OctetSimuCloneManager::recalcAllOctets(RunAccumulator& RA, bool doPlots) {
	return recalcOctets(RA, Octet::loadOctets(QFile(getEnvSafe("UCNA_OCTET_LIST"))), doPlots);
}

void OctetSimuCloneManager::setSimData(Sim2PMT* s2p) {
	if(simData && ownSimData) delete simData;
	ownSimData=false;
	simData = s2p;
}

void OctetSimuCloneManager::setOctetSimdata(unsigned int octn) {
	if(simData && ownSimData) delete simData;
	ownSimData = true;
	G4toPMT* G2P = new G4toPMT();
	for(unsigned int i=0; i<stride; i++)
		G2P->addFile(simFile+itos((stride*octn+i)%nTot)+".root");
	G2P->PGen[EAST].xscatter = G2P->PGen[WEST].xscatter = 0.005;
	G2P->runCathodeSim();
	simData = G2P;
}

void OctetSimuCloneManager::simOct(RunAccumulator& SimRA, const Octet& oct) {
	if(!oct.getNRuns()) return;
	smassert(simData);
	RunAccumulator* octSim = (RunAccumulator*)SimRA.makeAnalyzer(oct.octName(),"");
	octSim->simuClone(getEnvSafe("UCNA_ANA_PLOTS")+"/"+outputDir+"/"+oct.octName(), *simData, simFactor, hoursOld*3600, doPlots, doCompare);
	delete octSim;
}

void OctetSimuCloneManager::simOct(RunAccumulator& SimRA, unsigned int octn) {
		Octet oct = Octet::loadOctet(QFile(getEnvSafe("UCNA_OCTET_LIST")),octn);
		setOctetSimdata(octn);
		simOct(SimRA,oct);
}

void OctetSimuCloneManager::combineSims(RunAccumulator& SimRA, RunAccumulator* OrigRA) {
	SimRA.mergeSims(getEnvSafe("UCNA_ANA_PLOTS")+"/"+outputDir, OrigRA);
}
