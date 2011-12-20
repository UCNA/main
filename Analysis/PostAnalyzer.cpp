#include "PostAnalyzer.hh"
#include "PathUtils.hh"
#include <stdio.h>
#include <utility>

unsigned int PostAnalyzer::addRun(RunNum r) {
	
	std::string fileName = std::string("../OutputTree/Run_")+itos(r)+".root";
	if(!addFile(fileName))
		return 0;
	ProcessedDataScanner::addRun(r);

	if(withCals)
		PCals.insert(std::make_pair(r,new PMTCalibrator(r,CDB)));
	
	assert(false);
	//RunFile rf(r,true);
	//for(Side s=EAST; s<NONE; s=nextSide(s)) {
	//	if(!rf.getRuntime(s))
	//		printf("\n* Warning: zero length run for %i.%i\n",r,(int)s);
	//	totalTime.t[s] += rf.getRuntime(s);
	//}
	//runTimes.add(r,rf.getRuntime(BOTH));
	return 1;
}


void PostAnalyzer::setReadpoints() {
	char pl[] = {'X','Y'};
	char tmp[512];
	for(Side s=EAST; s<=WEST; s=nextSide(s)) {
		
		for(int p=0; p<2; p++) {
			sprintf(tmp,"Wires_%c%c",sideNames(s),pl[p]);
			Tch->SetBranchAddress(tmp,&wires[s][p]);
		}
		
		sprintf(tmp,"MWPC_%c",sideNames(s));
		Tch->SetBranchAddress(tmp,&mwpcs[s]);
		
		sprintf(tmp,"BetaSc%c",sideNames(s));
		Tch->SetBranchAddress(tmp,&scints[s]);
		sprintf(tmp,"BetaSc%c_led_pd",sideNames(s));
		Tch->SetBranchAddress(tmp,&led_pd[s]);
	}
	Tch->SetBranchAddress("Trigger",&trig);	
}
