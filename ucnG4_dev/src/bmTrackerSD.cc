//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: bmTrackerSD.cc,v 1.11 2011-10-12 19:31:16 mmendenhall Exp $
// GEANT4 tag $Name:  $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "bmTrackerSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4VProcess.hh"
#include "G4LossTableManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4Gamma.hh"
#include "SMExcept.hh"
#include "strutils.hh"
#include <cmath>
#include <cassert>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

bmTrackerSD::bmTrackerSD(G4String name): G4VSensitiveDetector(name) {
	collectionName.insert("trackerCollection");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void bmTrackerSD::Initialize(G4HCofThisEvent* HCE) {
	// make a new hits collection and register it for this event
	assert(collectionName.size());
	trackerCollection = new bmTrackerHitsCollection(SensitiveDetectorName,collectionName[0]); 
	G4int HCID = G4SDManager::GetSDMpointer()->GetCollectionID(trackerCollection); 
	HCE->AddHitsCollection(HCID, trackerCollection); 
	tracks.clear();
	originEnergy.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// quenching calculation... see Junhua's thesis
double quenchFactor(double E) {
	const G4double kb = 0.01907*cm/MeV;			// Birk's law quenching constant
	const G4double a = 116.7*MeV*cm*cm/g;		// dEdx fit parameter a*e^(b*E)
	const G4double b = -0.7287;					// dEdx fit parameter a*e^(b*E)
	const G4double rho = 1.032*g/cm3;			// scintillator density
	const G4double dEdx = a*rho*pow(E/keV,b);	// estimated dE/dx
	return 1.0/(1+kb*dEdx);
}

//If the track is already stored, simply update dedx
//otherwise add a new entry into the hit collection
G4bool bmTrackerSD::ProcessHits(G4Step* aStep,G4TouchableHistory*) {
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
	std::map<G4int,bmTrackerHit*>::iterator myTrack = tracks.find(thisTrackID);
	if(myTrack==tracks.end()) {
		bmTrackerHit* newHit = new bmTrackerHit();
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
		tracks.insert(std::pair<G4int,bmTrackerHit*>(thisTrackID,(bmTrackerHit*)trackerCollection->GetHit(hitn-1)));
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






//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void bmTrackerSD::EndOfEvent(G4HCofThisEvent*) {
	if (verboseLevel>0) { 
		G4int NbHits = trackerCollection->entries();
		G4cout << "\n-------->Hits Collection: in this event they are " << NbHits 
		<< " hits in the tracker chambers: " << GetName() << G4endl;
		for (G4int i=0;i<NbHits;i++) (*trackerCollection)[i]->Print();
	} 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

