#include "BetaSpectrum.hh"
#include "Enums.hh"
#include "PathUtils.hh"
#include "SMExcept.hh"

#include <TRandom.h>
#include <TF2.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <TString.h>
#include <unistd.h>
#include <vector>
#include <algorithm>
#include <cassert>
#include <numeric>
#include <bitset>  // djaffe 2aug06

#include "AnalysisManager.hh"
#include "SurfaceGenerator.hh"
#include "PrimaryGeneratorAction.hh"

#include <globals.hh>
#include <Randomize.hh>
#include <G4RunManager.hh>
#include <G4ParticleTable.hh>
#include <G4ParticleDefinition.hh>

/// generate a random position in a disk
void diskRandom(G4double radius, G4double& x, G4double& y) {
	while(true) {
		x = (2.0*G4UniformRand()-1.)*radius;
		y = (2.0*G4UniformRand()-1.)*radius;
		if(x*x+y*y<=radius*radius) break;
	}
}

PrimaryGeneratorAction::PrimaryGeneratorAction(DetectorConstruction* myDC):
myDetector(myDC), ETS(NULL), sourceRadius(0), relToSourceHolder(false) {
	particleGun = new G4ParticleGun();
	myMessenger = new PrimaryGeneratorMessenger(this);
	
	// default to electron
	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	G4ParticleDefinition* particle = particleTable->FindParticle("e-");
	particleGun->SetParticleDefinition(particle);
	particleGun->SetParticleTime(0.0*ns);
}

PrimaryGeneratorAction::~PrimaryGeneratorAction() {
	delete particleGun;
	if(ETS) delete ETS;
}

void PrimaryGeneratorAction::throwEvents(const std::vector<NucDecayEvent>& evts, G4Event* anEvent) {
	if(!evts.size()) return;
	
	G4ThreeVector direction;
	G4ThreeVector vtx;
	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	double wavg = 1.0;
	for(std::vector<NucDecayEvent>::const_iterator it = evts.begin(); it != evts.end(); it++) {
		if(it->d == D_ELECTRON) particleGun->SetParticleDefinition(particleTable->FindParticle("e-"));
		else if(it->d == D_GAMMA) particleGun->SetParticleDefinition(particleTable->FindParticle("gamma"));
		else continue;
		for(AxisDirection d = X_DIRECTION; d <= Z_DIRECTION; ++d) {
			direction[d] = it->p[d];
			vtx[d] = it->x[d]*m;
		}
		wavg *= it->w;
		particleGun->SetParticleEnergy(it->E*keV);
		particleGun->SetParticleMomentumDirection(direction);
		particleGun->SetParticlePosition(vtx);
		particleGun->SetParticleTime(it->t*s);
		displayGunStatus();
		particleGun->GeneratePrimaryVertex(anEvent);
	}
	// record event weight
	if(wavg != 1) {
		PrimEvtWeighting* w = new PrimEvtWeighting(pow(wavg,1./evts.size()));
		anEvent->SetUserInformation(w);
	}
}

void PrimaryGeneratorAction::SetEventFile(G4String val) {
	printf("Setting event generator input from '%s'\n",val.data());
	if(ETS) delete ETS;
	ETS = NULL;
	if(val=="") return;
	ETS = new EventTreeScanner();
	ETS->addFile(val.data());
}

void PrimaryGeneratorAction::displayGunStatus() {
	G4cout << particleGun->GetParticleDefinition()->GetParticleName() << " gun from " << particleGun->GetParticlePosition()/m
	<< "m towards " << particleGun->GetParticleMomentumDirection() << " at " << particleGun->GetParticleTime()/ns
	<< "ns : " << particleGun->GetParticleEnergy()/keV << "keV" << G4endl;
}

void PrimaryGeneratorAction::setVertices(std::vector<NucDecayEvent>& v) {
	G4ThreeVector v0 = posOffset;
	if(relToSourceHolder) {
		v0 += myDetector->getHolderPos();
	}
	if(sourceRadius) {
		G4double x,y;
		diskRandom(sourceRadius,x,y);
		v0[X_DIRECTION] += x;
		v0[Y_DIRECTION] += y;
	}
	for(unsigned int i=0; i<v.size(); i++)
		for(AxisDirection d = X_DIRECTION; d <= Z_DIRECTION; ++d)
			v[i].x[d] += v0[d];
}

void PrimaryGeneratorAction::initEventRandomSeed(G4Event* anEvent) {
	std::bitset<32> a ( (long) 0);
	//use the UI fRunNumber in AnalysisManager which got set
	//in RunAction, instead of the GEANT4 default RunID
	std::bitset<32> run ( gAnalysisManager->GetRunNumber() );
	std::bitset<32> evt ( anEvent->GetEventID() );
	
	//this is the "host id", unique for each host
	std::bitset<32> site ( (long) gethostid() ) ;
	
	for ( int i = 32-1-1; i >= 0 ; i-- ) {
		if ( run.test(i) ) a.set(31-1-i) ;
	}
	
	// create seed = (a XOR site) XOR evt
	std::bitset<32> seed = (a^site)^evt ;
	
	// set highest bit to zero to avoid negative seed
	if ( seed.test(31) ) seed.reset(31) ;
	
	myseed = seed.to_ulong();
	CLHEP::HepRandom::setTheSeed(myseed);	// random seed for Geant
	gRandom->SetSeed(myseed);		// random seed for ROOT
	G4cout<<"run "<<gAnalysisManager->GetRunNumber()<<" evt "<<anEvent->GetEventID()<<" seed "<<myseed<<G4endl;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {

	initEventRandomSeed(anEvent);
		
	if(ETS) {
		static std::vector<NucDecayEvent> v;
		v.clear();
		ETS->loadEvt(v);
		if(!ETS->firstpass) {
			// quit processing once new events have been used up
			G4cout << "Event input file exhausted; stopping run." << G4endl;
			G4RunManager::GetRunManager()->AbortRun();
		}
		setVertices(v);
		throwEvents(v,anEvent);
	} else {
		displayGunStatus();
		particleGun->GeneratePrimaryVertex(anEvent);
	}
	
	gAnalysisManager->FillPrimaryData(anEvent,myseed);
}
