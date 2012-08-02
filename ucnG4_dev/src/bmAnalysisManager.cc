///This is a simple class that will handle the MC-ROOT interface.
///The main point is to create a global interface object which 
///renders the flexibility for user "rooting" at various level 
///in the MC, e.g. detector construction, tracking, etc.
///Created based on IcaG4 from http://root.cern.ch/root/HowtoMC.html
///Jianglai 05/26/2006
///Updated MPM 10/2011

#include "G4ios.hh"
#include "G4Run.hh"
#include "G4Event.hh"
#include "G4Track.hh"
#include "G4VVisManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4IonTable.hh"
#include "G4SDManager.hh"
#include "G4PrimaryVertex.hh"
#include "G4PrimaryParticle.hh"
#include "G4SDManager.hh"
#include "G4EventManager.hh"
#include "Randomize.hh"
#include "SMExcept.hh"

#include "bmAnalysisManager.hh"
#include "bmPrimaryGeneratorAction.hh"

#include <vector>
#include <string>
#include <cassert>

using namespace CLHEP;
using namespace std;

//a global analysis manager
bmAnalysisManager *gbmAnalysisManager = (bmAnalysisManager *)0;

////////////////////////////////////////////////////////////////////////////////////
bmAnalysisManager::bmAnalysisManager(): pMcEvent(&mcEvent) {
	fROOTOutputFile=0;
	fEventTree=0;
	if (gbmAnalysisManager)
		delete gbmAnalysisManager;
	gbmAnalysisManager = this;
	fRunNumber = 0; //initialized to zero
	G4cout<<"Analysis manager constructed "<<G4endl;
}

void bmAnalysisManager::OpenFile(const G4String filename) {
	if(fROOTOutputFile != 0)
		CloseFile();  
	
	G4cout<<"Opening root file "<<filename<<G4endl;
	
	fROOTOutputFile= new TFile(filename.c_str(), "RECREATE", "Geant4 benchmark simulation output file");  
	
	if (fROOTOutputFile == 0) {
		SMExcept e("CannotOpenOutputFile");
		e.insert("filename",filename);
		throw(e);
    }
	CreateTrees();
}

void bmAnalysisManager::CloseFile() { 
	if(fROOTOutputFile != 0) {
		G4cout << "Closing " << fROOTOutputFile->GetName() << G4endl;
		fROOTOutputFile->cd();
		fEventTree->Write();
		delete fROOTOutputFile;
		fROOTOutputFile=0;
		fEventTree=0;
    }
}

void bmAnalysisManager::CreateTrees() {
	fEventTree = new TTree("EventTree","tree of MC event");
	fEventTree->SetMaxVirtualSize(10000000);
	fEventTree->Branch("MC_event_output","bmMCEvent",&pMcEvent,64000,99); 
}

//do this only once for one run
void bmAnalysisManager::StoreHitCollectionIDs() {
	G4SDManager * SDman = G4SDManager::GetSDMpointer();
	detectorIDs.clear();
	for(size_t ii=0;ii<fSDNames.size();ii++){
		int id = SDman->GetCollectionID(fSDNames[ii]+"/trackerCollection");
		if(id<0) {
			SMExcept e("BadCollectionID");
			e.insert("id",id);
			e.insert("name",fSDNames[ii]+"/trackerCollection");
			throw(e);
		}
		detectorIDs.push_back(id);
		G4cout<<fSDNames[ii]<<" : Id = "<<id<<G4endl;
	} 
}

void bmAnalysisManager::FillPrimaryData(const G4Event* evt_in, const long seed) {
	mcEvent.eventID = evt_in->GetEventID();
	double w = 1.0;
	if(evt_in->GetUserInformation())
		w = ((PrimEvtWeighting*)evt_in->GetUserInformation())->w;
	
	// convert each primary to ROOT-friendly form and store
	for (int ivert=0; ivert<evt_in->GetNumberOfPrimaryVertex(); ivert++){
		bmPrimaryInfo prim_info;
		G4PrimaryVertex* pv= evt_in->GetPrimaryVertex(ivert);
		G4PrimaryParticle* pp = pv->GetPrimary();
		for(unsigned int i=0; i<3; i++) {
			prim_info.vertex[i] = pv->GetPosition()(i)/m;
			prim_info.p[i] = pp->GetMomentum()(i)/keV;
		}
		G4double momentum = pp->GetMomentum().mag();
		G4double mass = pp->GetMass();
		G4double energy = sqrt(momentum*momentum+mass*mass);
		prim_info.KE = (energy - mass)/keV;
		prim_info.seed = seed;
		prim_info.weight = w;
		mcEvent.AddPrimaryInfo(prim_info);
	}
}

void bmAnalysisManager::FillTrackerData(const G4Event *evt) {
	
	G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
	bmSteppingAction* USA = (bmSteppingAction*)G4EventManager::GetEventManager()->GetUserSteppingAction();
	
	for(size_t nn=0; nn<detectorIDs.size();nn++){
		if (HCE) {
			bmTrackerHitsCollection* HC_detector = (bmTrackerHitsCollection*)(HCE->GetHC(detectorIDs[nn]));
			G4int n_hit = 0;
			if(HC_detector!=NULL){
				n_hit = HC_detector->entries();
				for(int ii=0;ii<n_hit;ii++){
					bmTrackInfo track_info;
					track_info.hcID = detectorIDs[nn];
					((*HC_detector)[ii])->fillTrackInfo(track_info);
					mcEvent.AddTrackInfo(track_info);
				}
			}
		}
	}
	
	mcEvent.trapped = USA->GetTrappedFlag();
	mcEvent.compTime = USA->GetTimeSpent();
}

void bmAnalysisManager::FillEventTree() {
	fEventTree->Fill();
	mcEvent.ClearEvent();
}
