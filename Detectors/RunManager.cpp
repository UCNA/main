#include "RunManager.hh"
#include "PathUtils.hh"
#include "Subsystem.hh"
#include "EnergyCalibrator.hh"
#include "CalDBSQL.hh"
#include "CalDBFake.hh"
#include "strutils.hh"

RunManager::RunManager(RunNum r, RunMode RM, std::string bp, bool noload): OutputManager("RunManager",""), runMode(RM), manualAnalyze(false),
RI(CalDBSQL::getCDB()->getRunInfo(r)), SI(SourcesInfo(r)) {
	
	plotPath = basePath = bp+"/Plots/"+replace(RI.groupName,' ','_')+"/"+itos(r)+"_"+RI.roleName+"/";
	dataPath = bp+"/RunData/";
	rootPath = bp+"/OutputTree/";
	makePath(dataPath);
	qOut.setOutfile(dataPath+"Run_"+itos(RI.runNum)+".txt");
	
	Stringmap M;
	M.insert("basePath",basePath);
	M.insert("plotPath",plotPath);
	M.insert("dataPath",dataPath);
	M.insert("rootPath",rootPath);
	qOut.insert("RunPaths",M);
	
	if(noload) {
		treeFile = NULL;
		myTree = NULL;
		nEvents = 0;
	} else {
		std::string infile = RunManager::getDataPath(r);
		assert(infile.size());
		printf("Loading run %s...\n",infile.c_str());
		treeFile = new TFile(infile.c_str(),"READ");
		myTree = (TTree*)(treeFile->Get("h1"));
		nEvents = myTree->GetEntries();
	}
		
	CDB = CalDBSQL::getCDB();
	if(!CDB->isValid(r)) {
		printf("\n**** WARNING: Fake Calibrations in use for run %i! ****\n\n",r);
		assert(IGNORE_DEAD_DB);
		CDB = new CalDBFake();
		Stringmap M;
		M.insert("level","severe");
		M.insert("description","Cal_DB_Missing");
		qOut.insert("Warning",M);
	}
	PCal = new PMTCalibrator(r,CDB);
	PCal->printSummary();
}

bool RunManager::checkBranch(std::string s) const {
	bool hasBranch =  myTree->GetBranch(s.c_str());
	if(!hasBranch)
		printf("**** Data tree does not contain branch '%s'!\n",s.c_str());
	return hasBranch;
}

float* RunManager::readBranch(std::string bname) {
	float* bdat = new float[nEvents];
	TBranch* b = myTree->GetBranch(bname.c_str());
	assert(b);
	float x;
	b->SetAddress(&x);
	b->LoadBaskets();
	for(unsigned int e=0; e<nEvents; e++) {
		b->GetEntry(e);
		bdat[e] = x;
	}
	return bdat;
}

std::string RunManager::getDataPath(RunNum rn) {
	std::vector<std::string> paths;
	paths.push_back("../rootfiles");						// universal, make a link to correct location from here
	paths.push_back("/data/ucnadata/2010/rootfiles");		// 2010 data on UCN
	paths.push_back("/data/ucnadata/2011/rootfiles");		// new data on UCN
	paths.push_back("../SampleRuns");						// on Praetorius
	paths.push_back("/home/data_analyzed/2011/rootfiles");	// on PCUCN*
	paths.push_back("/home/data_analyzed/2010/rootfiles");	// on PCUCN*
	paths.push_back("/home/data_analyzed/2009/rootfiles");	// on PCUCN*
	for(unsigned int i=0; i<paths.size(); i++) {
		std::string testpath = paths[i]+"/full"+itos(rn)+".root";
		if(fileExists(testpath)) return testpath;
	}
	return "";
}

void RunManager::write(std::string) {
	// record run and sources info
	qOut.insert("RunInfo",RI.toStringmap());
	for(Side m=EAST; m<=WEST; m = nextSide(m))
		for(std::vector<Source>::const_iterator it = SI.sourcesBegin(m); it != SI.sourcesEnd(m); it++)
			qOut.insert("Source", it->toStringmap());
	
	OutputManager::write(std::string("Run_")+itos(RI.runNum)+".txt");
}

void RunManager::makeRootOut(Subsystem* trigger) {
	
	if(!(runMode & RUNMODE_PLOTSOUT || runMode & RUNMODE_TREEOUT))
		return;
	
	printf("\n--------- Building output .root file... ----------\n");
	makePath(rootPath);
	char tmp[512];
	sprintf(tmp,"%s/Run_%i.root",rootPath.c_str(),RI.runNum);
	TFile* outputFile = new TFile(tmp,"RECREATE");
	outputFile->cd();
	
	if(runMode & RUNMODE_PLOTSOUT)
		writeItems();
	
	
	if(trigger && (runMode & RUNMODE_TREEOUT)) {
		// set up output tree
		printf("Generating processed events tree...\n");
		TTree* T = new TTree("OutTree","Summary variables tree");
		for(std::vector<Subsystem*>::iterator it = subsystems.begin(); it != subsystems.end(); it++) {
			(*it)->addOutBranches(T);
			printf("%s\n",(*it)->name.c_str());
		}
		
		// fill events
		for(UInt_t e = 0; e<nEvents; e++) {
			if(!trigger->triggered(e))
				continue;
			for(std::vector<Subsystem*>::iterator it = subsystems.begin(); it != subsystems.end(); it++)
				(*it)->fillEvent(e);
			T->Fill();
		}
		T->Write();
		delete(T);
		printf("Done.\n");
	}

	outputFile->Close();
	

	printf("---------          Done.          ----------\n");
}

