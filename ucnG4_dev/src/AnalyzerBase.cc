#include "AnalyzerBase.hh"

ucnG4_analyzer::ucnG4_analyzer(const std::string& outfname): anaTree(NULL),
outf(new TFile(outfname.c_str(),"RECREATE")), myevt(new bmMCEvent()) { }

void ucnG4_analyzer::analyzeFileList(const string& flist) {
	ifstream file;
	file.open(flist.c_str());
	string fname;
	while(!file.fail() && !file.eof()){
		fname = "";
		file >> fname;
		if(fname.size())
			analyzeFile(fname);
	}  
}

ucnG4_analyzer::~ucnG4_analyzer() {
	if(outf) {
		outf->cd();
		anaTree->Write();
		outf->Close();
		//delete outf; //TODO
	}
	//delete anaTree; //TODO
	//delete myevt; //TODO
}

void ucnG4_analyzer::analyzeFile(const string& fname) {
	
	if(!anaTree) {
		outf->cd();
		anaTree = new TTree("anaTree", "tree for analysis");
		anaTree->Branch("primKE",&primKE,"primKE/D");
		anaTree->Branch("primTheta",&primTheta,"primTheta/D");
		anaTree->Branch("primPos",primPos,"primPos[4]/D");
		anaTree->Branch("primWeight",&primWeight,"primWeight/D");
		anaTree->Branch("trapped",&fTrapped,"trapped/I");
		anaTree->Branch("compTime",&fCompTime,"compTime/D");
		anaTree->Branch("seed",&seed,"seed/I");	
		setupOutputTree();		
	}
	
	cout << "opening " << fname << endl;
	TFile f(fname.c_str(),"read");
	if(f.IsZombie()) {
		cout << "*** Error opening" << fname << endl;
		return;
	}

	TTree* tree = (TTree*) f.Get("EventTree");
	tree->SetBranchAddress("MC_event_output",&myevt);
	
	for(int ii=0;ii<tree->GetEntriesFast();ii++) {		
				
		// reset analyzer variables
		seed=fTrapped=0;
		primKE=primTheta=0;
		for(unsigned int i=0; i<4; i++) primPos[i]=0;
		resetAnaEvt(); 
		
		// get event data
		tree->GetEntry(ii);
		fTrapped = myevt->trapped;
		fCompTime = myevt->compTime;
		TClonesArray& priminfo_array = *(myevt->primaryInfo);    
		TClonesArray& trackinfo_array = *(myevt->trackInfo);
		
		// scan each track segment
		Int_t nTracks = trackinfo_array.GetEntries();			
		for(int n=0; n<nTracks; n++) {
			trackinfo = (bmTrackInfo*)trackinfo_array[n];
			priminfo = (bmPrimaryInfo*)priminfo_array[n];
			pID = trackinfo->pID;
			detectorID = trackinfo->hcID;
			trackID = trackinfo->trackID;
			processTrack();
		}
		
		// primary event info
		priminfo = (bmPrimaryInfo*)priminfo_array[0];
		primKE = priminfo->KE;
		primWeight = priminfo->weight;
		for(unsigned int i=0; i<3; i++)
			primPos[i] = priminfo->vertex[i];
		primPos[3] = sqrt(primPos[0]*primPos[0]+primPos[1]*primPos[1]);
		TVector3 mom(priminfo->p[0],priminfo->p[1],priminfo->p[2]);
		primTheta = mom.Theta();
		seed = priminfo->seed;

		processEvent();
		if(saveEvent())
			anaTree->Fill();
	}
		
	f.Close();
}
