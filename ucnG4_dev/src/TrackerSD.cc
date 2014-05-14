#include "SMExcept.hh"
#include "strutils.hh"
#include <cmath>
#include <cassert>

#include "TrackerSD.hh"

#include <G4SystemOfUnits.hh>
#include <G4HCofThisEvent.hh>
#include <G4Step.hh>
#include <G4ThreeVector.hh>
#include <G4SDManager.hh>
#include <G4ios.hh>
#include <G4VProcess.hh>
#include <G4LossTableManager.hh>
#include <G4ParticleDefinition.hh>
#include <G4Gamma.hh>

TrackerSDMessenger::TrackerSDMessenger(TrackerSD* T): mySD(T) {
	sdDir = new G4UIdirectory(("/SD/"+mySD->GetName()+"/").c_str());
	sdDir->SetGuidance("Sensitive detector response settings");
	
	kbCmd = new G4UIcmdWithADouble((sdDir->GetCommandPath()+"kb").c_str(), this);
	kbCmd->SetGuidance("Birk's Law quenching constant in cm/MeV");
	kbCmd->SetDefaultValue(0.01907);
	kbCmd->AvailableForStates(G4State_Idle);
}

TrackerSDMessenger::~TrackerSDMessenger() {
	delete kbCmd;
	delete sdDir;
}

void TrackerSDMessenger::SetNewValue(G4UIcommand* command, G4String newValue) {
	if( command == kbCmd ) {
		G4double k = kbCmd->GetNewDoubleValue(newValue);
		G4cout << "Setting Birk's Law kb = " << k << " cm/MeV for " << mySD->GetName() << G4endl;
		mySD->SetKb(k * cm/MeV);
	}
}

//----------------------------------------------------------------

TrackerSD::TrackerSD(G4String name): G4VSensitiveDetector(name), kb(0.01907*cm/MeV), rho(1.032*g/cm3) {
	new TrackerSDMessenger(this);
	collectionName.insert("trackerCollection");
}

void TrackerSD::Initialize(G4HCofThisEvent* HCE) {
	// make a new hits collection and register it for this event
	assert(collectionName.size());
	trackerCollection = new TrackerHitsCollection(SensitiveDetectorName,collectionName[0]); 
	G4int HCID = G4SDManager::GetSDMpointer()->GetCollectionID(trackerCollection); 
	HCE->AddHitsCollection(HCID, trackerCollection); 
	tracks.clear();
	originEnergy.clear();
}

// quenching calculation... see Junhua's thesis
double TrackerSD::quenchFactor(double E) const {
	const G4double a = 116.7*MeV*cm*cm/g;		// dEdx fit parameter a*e^(b*E)
	const G4double b = -0.7287;					// dEdx fit parameter a*e^(b*E)
	const G4double dEdx = a*rho*pow(E/keV,b);	// estimated dE/dx
	return 1.0/(1+kb*dEdx);
}

//If the track is already stored, simply update dedx
//otherwise add a new entry into the hit collection
G4bool TrackerSD::ProcessHits(G4Step* aStep,G4TouchableHistory*) {
	assert(aStep);
	G4Track* aTrack = aStep->GetTrack();
	assert(aTrack);

	G4String creator_proc = "";
	// Check if the track has the creator process (not the case for primaries)
	const G4VProcess* creatorProcess = aTrack->GetCreatorProcess();
	if(!creatorProcess) creator_proc="original";
	else creator_proc = creatorProcess->GetProcessName();  
	
	G4StepPoint* preStep = aStep->GetPreStepPoint();
	assert(preStep);
	G4StepPoint* postStep = aStep->GetPostStepPoint();
	assert(postStep);
	
	G4ThreeVector prePos = preStep->GetPosition();
	G4ThreeVector postPos = postStep->GetPosition();
	G4double E0 = preStep->GetKineticEnergy();
	G4double E1 = postStep->GetKineticEnergy();
	G4double Ec = 0.5*(E0+E1);

	// get prior track, or initialize a new one
	G4int thisTrackID = aTrack->GetTrackID();
	std::map<G4int,TrackerHit*>::iterator myTrack = tracks.find(thisTrackID);
	if(myTrack==tracks.end()) {
		TrackerHit* newHit = new TrackerHit();
		newHit->SetTrackID(thisTrackID);
		newHit->SetPID(aTrack->GetDefinition()->GetPDGEncoding());
		newHit->SetProcessName(creator_proc);
		newHit->SetIncidentEnergy(preStep->GetKineticEnergy());
		newHit->SetPos(postPos);
		newHit->SetHitTime(preStep->GetGlobalTime());
		newHit->SetIncidentMomentum(preStep->GetMomentum());
		G4VPhysicalVolume* preVolume = preStep->GetPhysicalVolume();
		newHit->SetVolumeName(preVolume?preVolume->GetName():"Unknown");
		newHit->SetVertex(aTrack->GetVertexPosition());
		newHit->SetCreatorVolumeName(aTrack->GetLogicalVolumeAtVertex()->GetName());
		newHit->nSecondaries = 0;
		std::map<const G4Track*,double>::iterator itorig = originEnergy.find(aTrack);
		if(itorig == originEnergy.end()) {
			// originEnergy = 0 for primary tracks (not a sub-track of another track in this volume)
			newHit->originEnergy = 0;
		} else {
			// get previously stored origin energy for secondary tracks; remove listing entry
			newHit->originEnergy = itorig->second;
			originEnergy.erase(itorig);
		}
		unsigned int hitn = trackerCollection->insert(newHit);
		tracks.insert(std::pair<G4int,TrackerHit*>(thisTrackID,(TrackerHit*)trackerCollection->GetHit(hitn-1)));
		myTrack = tracks.find(thisTrackID);
	}
	
	// accumulate edep, edepq, local position for this step
	G4double edep = aStep->GetTotalEnergyDeposit();
	G4double edepQ = edep*quenchFactor(myTrack->second->originEnergy==0?Ec:myTrack->second->originEnergy);
	G4ThreeVector localPosition = preStep->GetTouchableHandle()->GetHistory()->GetTopTransform().TransformPoint(prePos);
	myTrack->second->AddEdep(edep,localPosition);
	myTrack->second->AddEdepQuenched(edepQ);
	myTrack->second->SetExitMomentum(postStep->GetMomentum());
	
	// record origin energy for secondaries in same volume
	const G4TrackVector* secondaries = aStep->GetSecondary();
	assert(secondaries);
	while(myTrack->second->nSecondaries < secondaries->size()) {
		const G4Track* sTrack = (*secondaries)[myTrack->second->nSecondaries++];
		if(sTrack->GetVolume() != aTrack->GetVolume())
			continue;
		const G4double eOrig = myTrack->second->originEnergy>0?myTrack->second->originEnergy:Ec;
		if(originEnergy.find(sTrack) != originEnergy.end()) {
			SMExcept e("duplicateSecondary");
			e.insert("eOrig",eOrig);
			e.insert("pID",myTrack->second->GetPID());
			e.insert("nSec",myTrack->second->nSecondaries++);
			e.insert("eOrig_old",originEnergy.find(sTrack)->second);
			throw(e);
		}
		originEnergy.insert(std::pair<const G4Track*,double>(sTrack,eOrig));
	}
	
	return true;
}

void TrackerSD::EndOfEvent(G4HCofThisEvent*) {
	if (verboseLevel>0) { 
		G4int NbHits = trackerCollection->entries();
		G4cout << "\n-------->Hits Collection: in this event they are " << NbHits 
		<< " hits in the tracker chambers: " << GetName() << G4endl;
		for (G4int i=0;i<NbHits;i++) (*trackerCollection)[i]->Print();
	} 
}
