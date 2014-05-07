#include "AnalyzerBase.hh"

ucnG4_analyzer::ucnG4_analyzer(const std::string& outfname): anaTree(NULL),
outf(new TFile(outfname.c_str(),"RECREATE")), myevt(new MCEvent()) { }

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
		anaTree->Branch("primKE",&vtx.primKE,"primKE/D");
		anaTree->Branch("primTheta",&vtx.primTheta,"primTheta/D");
		anaTree->Branch("primPos",vtx.primPos,"primPos[4]/D");
		anaTree->Branch("primWeight",&vtx.primWeight,"primWeight/D");
		anaTree->Branch("trapped",&vtx.fTrapped,"trapped/I");
		anaTree->Branch("compTime",&vtx.fCompTime,"compTime/D");
		anaTree->Branch("seed",&vtx.seed,"seed/I");
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
	
	for(int ii=0; ii<tree->GetEntriesFast(); ii++) {
				
		// reset analyzer variables
		vtx.seed = vtx.fTrapped = 0;
		vtx.primKE = vtx.primTheta = 0;
		for(unsigned int i=0; i<4; i++) vtx.primPos[i]=0;
		pass_number = 1;
		resetAnaEvt(); 
		
		// get event data
		tree->GetEntry(ii);
		vtx.fTrapped = myevt->trapped;
		vtx.fCompTime = myevt->compTime;
		TClonesArray& priminfo_array = *(myevt->primaryInfo);    
		TClonesArray& trackinfo_array = *(myevt->trackInfo);
		
		// primary event info, summed over primaries where applicable
		vtx.primKE = 0;
		for(unsigned int n=0; n<priminfo_array.GetEntries(); n++) {
			priminfo = (PrimaryInfo*)priminfo_array[n];
			vtx.primKE += priminfo->KE;
			
			if(n>0) continue;
			
			vtx.primWeight = priminfo->weight;
			for(unsigned int i=0; i<3; i++)
				vtx.primPos[i] = priminfo->vertex[i];
			vtx.primPos[3] = sqrt(vtx.primPos[0]*vtx.primPos[0]+vtx.primPos[1]*vtx.primPos[1]);
			TVector3 mom(priminfo->p[0],priminfo->p[1],priminfo->p[2]);
			vtx.primTheta = mom.Theta();
			vtx.seed = priminfo->seed;
		}
		
		// scan each track segment
		Int_t nTracks = trackinfo_array.GetEntries();
		while(pass_number) {
			for(int n=0; n<nTracks; n++) {
				trackinfo = (TrackInfo*)trackinfo_array[n];
				priminfo = (PrimaryInfo*)priminfo_array[n];
				pID = trackinfo->pID;
				detectorID = trackinfo->hcID;
				trackID = trackinfo->trackID;
				processTrack();
			}
			processEvent();
			--pass_number;
		}
		
		writeEvent();			
	}
		
	f.Close();
}
