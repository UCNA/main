#include "SiDet_Analyzer.hh"
#include <TMath.h>
#include <cfloat>

void SiDet_Analyzer::setupOutputTree() {
	printf("Adding branches for SiDet_Analyzer...\n");
	anaTree->Branch("Edep",&Edep,"Edep/D");
	anaTree->Branch("Pos",Pos,"Pos[3]/D");
}

void SiDet_Analyzer::resetAnaEvt() {
	Edep = 0;
	for(unsigned int i=0; i<3; i++) Pos[i] = 0;
}

void SiDet_Analyzer::processTrack() {
	
	// detector ID numbers
	const int ID_Si = 0;

	// scintillator deposited energy, position, hit time
	if(detectorID==ID_Si) {
		Edep += trackinfo->Edep;
		for(unsigned int j=0; j<3; j++) Pos[j] += trackinfo->edepPos[j];
	}
	
}

void SiDet_Analyzer::processEvent() {
	// normalize position variables
	for(int i = 0; i < 3; i++)
		if(Edep>0)
			Pos[i] /= Edep;
}

int main(int argc, char** argv) {
	
	if(argc<3) {
		cout<<"Syntax: "<<argv[0]<<" <filename containing list of raw root files> <output root file name>"<<endl;
		exit(1);
	}
	
	SiDet_Analyzer SDA(argv[2]);
	SDA.analyzeFileList(argv[1]);
	
	return 0;
}
