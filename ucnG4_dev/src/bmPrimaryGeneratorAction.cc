//// ********************************************************************
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
// $Id: bmPrimaryGeneratorAction.cc,v 1.34 2011-12-12 01:29:37 mmendenhall Exp $
// GEANT4 tag $Name:  $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "bmPrimaryGeneratorAction.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "BetaSpectrum.hh"
#include "globals.hh"
#include "Randomize.hh"
#include "bmAnalysisManager.hh"
#include "PathUtils.hh"
#include "SMExcept.hh"
#include <iostream>
#include <fstream>
#include <vector>
#include <TString.h>
#include <TRandom.h>
#include <TF2.h>

#include <unistd.h>
#include <vector>
#include <algorithm>
#include <cassert>
#include <numeric>
#include <bitset>  // djaffe 2aug06

namespace CLHEP {}
using namespace CLHEP; 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//for beta spectrum
TF1* funcBetaSpectrum = NULL;
TF1* funcHeavyBeta = NULL;

double rand_outof_list(const double* line, const double* weight, const int npoints, int& selected) {
	// build cumulative list
	std::vector<double> division(npoints+1,0.0);
	for(int i=0; i<npoints; i++)
		division[i+1] = division[i]+weight[i];
	// pick selected line
	std::vector<double>::iterator itsel = std::upper_bound(division.begin(),division.end(),division.back()*G4UniformRand());
	selected = int(itsel-division.begin())-1;	
	assert(0<=selected && selected<npoints);
	return line[selected];
}

void RandomizeMomentum(G4ThreeVector& mom) {
	// isotropic emission
	G4double phi = 2.0*M_PI*G4UniformRand();
	G4double costheta = 2.0*G4UniformRand()-1.0;
	G4double sintheta = sqrt(1.0-costheta*costheta);
	mom.set(cos(phi)*sintheta,sin(phi)*sintheta,costheta);
}  

/// generate a random position in a tube
void randomTubePosition(const G4ThreeVector centerpos, const G4double radius, const G4double halfz, G4ThreeVector& pos) {
	G4double x0,y0;
	while(true) {
		x0 = (2.0*G4UniformRand()-1.)*radius;
		y0 = (2.0*G4UniformRand()-1.)*radius;
		if(x0*x0+y0*y0<=radius*radius) break;
	}
	G4double z0=(2.0*G4UniformRand()-1.0)*halfz;
	pos = centerpos+G4ThreeVector(x0,y0,z0);
}

/// generate a random position with uniform number in radial bins
void randomUniformRadialBins(const G4ThreeVector centerpos, const G4double radius, const G4double halfz, G4ThreeVector& pos) {
	G4double r = G4UniformRand()*radius;
	G4double theta = 2*PI*G4UniformRand();
	G4double z0=(2.0*G4UniformRand()-1.0)*halfz;
	pos = centerpos+G4ThreeVector(r*cos(theta),r*sin(theta),z0);
}

double genericBetaSpectrum(double* x, double *par) {
	double KE = x[0]; // particle kinetic energy
	double Q = par[0]; // spectrum endpoint
	return plainPhaseSpace((KE+m_e)/m_e,(Q+m_e)/m_e);
}

double heavyBetaSpectrum(double* x, double* par) {
	double KE = x[0];	// beta kinetic energy
	double Q = par[0];	// spectrum endpoint
	double A = par[1];	// nucleus A
	double Z = par[2];	// nucleus Z
	double W = (KE+m_e)/m_e;
	double W0 = (Q+m_e)/m_e;
	double R = pow(A,1./3.)*neutron_R0;
	
	// TODO: recoil/weak magnetism terms???
	if(0<KE && KE<Q)
		return plainPhaseSpace(W,W0)*WilkinsonF0(Z,W,R)*WilkinsonL0(Z,W,R)*(1.+Wilkinson_g(W,W0));
	return 0;
}

void bmPrimaryGeneratorAction::
throwEvents(const std::vector<NucDecayEvent>& evts, G4Event* anEvent) {
	G4ThreeVector direction;
	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	for(std::vector<NucDecayEvent>::const_iterator it = evts.begin(); it != evts.end(); it++) {
		if(it->d == D_ELECTRON) particleGun->SetParticleDefinition(particleTable->FindParticle("e-"));
		else if(it->d == D_GAMMA) particleGun->SetParticleDefinition(particleTable->FindParticle("gamma"));
		else continue;
		G4cout << "\t" << it->d << " = " << it->E << "keV" << G4endl;
		RandomizeMomentum(direction);
		particleGun->SetParticleEnergy(it->E*keV);
		particleGun->SetParticleMomentumDirection(direction);
		particleGun->GeneratePrimaryVertex(anEvent);
	}
}

void bmPrimaryGeneratorAction::
throwElectronsAndGammas(const std::vector<G4double>& electrons,
						const std::vector<G4double>& gammas,
						G4Event* anEvent) {
	
	
	G4cout << "Throwing " << electrons.size() << " electrons and " << gammas.size() << " gammas." << G4endl;
	
	G4ThreeVector direction;
	
	// generate electrons
	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	particleGun->SetParticleDefinition(particleTable->FindParticle("e-"));
	for(std::vector<G4double>::const_iterator it = electrons.begin(); it != electrons.end(); it++) {
		G4cout << "\te- = " << *it/keV << "keV" << G4endl;
		RandomizeMomentum(direction);
		particleGun->SetParticleEnergy(*it);
		particleGun->SetParticleMomentumDirection(direction);
		particleGun->GeneratePrimaryVertex(anEvent);
	}
	
	// generate photons
	particleGun->SetParticleDefinition(particleTable->FindParticle("gamma"));
	for(std::vector<G4double>::const_iterator it = gammas.begin(); it != gammas.end(); it++) {	
		G4cout << "\tgamma = " << *it/keV << "keV" << G4endl;
		RandomizeMomentum(direction);
		particleGun->SetParticleEnergy(*it);
		particleGun->SetParticleMomentumDirection(direction);
		particleGun->GeneratePrimaryVertex(anEvent);
	}
}

// lines based on NuDat 2.6
void bmPrimaryGeneratorAction::Cd113mSourceGenerator(G4Event* anEvent) {
	
	std::vector<G4double> electrons;
	std::vector<G4double> gammas;
	
	funcHeavyBeta->FixParameter(0,585.7);
	funcHeavyBeta->FixParameter(1,113.);
	funcHeavyBeta->FixParameter(2,49.);
	electrons.push_back(funcHeavyBeta->GetRandom()*keV);
	
	throwElectronsAndGammas(electrons,gammas,anEvent);
}

void bmPrimaryGeneratorAction::In114SourceGenerator(G4Event* anEvent) {
	
	G4double gunEnergy;
	
	std::vector<G4double> electrons;
	std::vector<G4double> gammas;
	int selected(-1);
	
	if(G4UniformRand() < 0.5) {
		// 5+ decay to:				114Cd ------------	or 1+ state ------------------
		const double line190[] =	{725.24,	558.43,	190.27,	189.44,	186.03,	162.33};
		const double branch190[] =	{3.2,		3.2,	15.56,	6.71,	31.9,	40.1};
		gunEnergy = rand_outof_list(line190, branch190, 6, selected)*keV;
		if(selected < 3)
			gammas.push_back(gunEnergy);
		else
			electrons.push_back(gunEnergy);
		// auger K
		if(selected >= 2 && G4UniformRand() < 0.0598/0.94)
			electrons.push_back(20.1*keV);
	} else {
		// 1+ state 1988.7keV beta decay to 114Sn
		if(G4UniformRand() < 0.9936*0.94) {
			funcHeavyBeta->FixParameter(0,1988.7);
			funcHeavyBeta->FixParameter(1,114.);
			funcHeavyBeta->FixParameter(2,50.);
			electrons.push_back(funcHeavyBeta->GetRandom()*keV);
		}
	}
	
	throwElectronsAndGammas(electrons,gammas,anEvent);
}


bmPrimaryGeneratorAction::bmPrimaryGeneratorAction(bmDetectorConstruction* myDC): myDetector(myDC) {
	particleGun = new G4ParticleGun();
	gunMessenger = new bmPrimaryGeneratorMessenger(this);
	
	// default to electron
	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	G4ParticleDefinition* particle = particleTable->FindParticle("e-");
	particleGun->SetParticleDefinition(particle);
	particleGun->SetParticleTime(0.0*ns);
	
	// beta spectra TF1s
	funcBetaSpectrum = new TF1("funcBetaSpectrum",genericBetaSpectrum,0,3000,1);
	funcBetaSpectrum->SetNpx(3000);
	funcHeavyBeta = new TF1("funcHeavyBeta",heavyBetaSpectrum,0,5000,3);
	funcHeavyBeta->SetNpx(500);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

bmPrimaryGeneratorAction::~bmPrimaryGeneratorAction() {
	delete particleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void bmPrimaryGeneratorAction::displayGunStatus() {
	G4cout << "Electron gun from " << particleGun->GetParticlePosition()/m
	<< " towards " << particleGun->GetParticleMomentumDirection()
	<< ": " << particleGun->GetParticleEnergy()/keV << "keV" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void bmPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {
	
	std::bitset<32> a ( (long) 0);
	//use the UI fRunNumber in bmAnalysisManager which got set 
	//in bmRunAction, instead of the GEANT4 default RunID
	std::bitset<32> run ( gbmAnalysisManager->GetRunNumber() );
	std::bitset<32> evt ( anEvent->GetEventID() );
	
	//   if one wants to use the hostname somehow, this is the way to do it
	//   char hostname[256];
	//   gethostname(hostname,255);
	//   cout<<hostname<<endl;
	
	//this is the "host id", unique for each host
	std::bitset<32> site ( (long) gethostid() ) ; 
	
	for ( int i = 32-1-1; i >= 0 ; i-- ) {
		if ( run.test(i) ) a.set(31-1-i) ;
	}
	
	// create seed = (a XOR site) XOR evt
	std::bitset<32> seed = (a^site)^evt ;
	
	// set highest bit to zero to avoid negative seed
	if ( seed.test(31) ) seed.reset(31) ; 
	
	long myseed = seed.to_ulong();
	
	HepRandom::setTheSeed(myseed);	// random seed for Geant
	gRandom->SetSeed(myseed);		// random seed for ROOT
	G4cout<<"run "<<gbmAnalysisManager->GetRunNumber()<<" evt "<<anEvent->GetEventID()<<" seed "<<myseed<<G4endl;
	
	if(particleType.size()) {
		G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
		G4ParticleDefinition* particle = particleTable->FindParticle(particleType);
		if(!particle) {
			SMExcept e("UnknownParticle");
			e.insert("name",particleType);
			throw(e);
		}
		particleGun->SetParticleDefinition(particle);
	}
	particleGun->SetParticleTime(0.0*ns);
	
	static NucDecayLibrary NDL(getEnvSafe("UCNA_AUX")+"/NuclearDecays/",1e-6);
	std::vector<NucDecayEvent> evts;
	
	// set vertex position for event according to position generator
	// assumed position has been set in /gun/position!
	static G4ThreeVector position_saved = particleGun->GetParticlePosition();
	G4ThreeVector vertex_position = position_saved;
	if(positioner=="Fixed") {
		// fixed position set from macro file
		vertex_position = position_saved;
	} else if(positioner=="SourceDrop") {
		// nominally 3mm disc around pre-set position for calibration sources
		randomTubePosition(position_saved, sourceRadius, 0*um, vertex_position); 
	} else if(positioner=="DecayTrapUniform") {
		// uniform fill of decay trap volume
		randomTubePosition(G4ThreeVector(0,0,0), 62.2*mm, 1.5*m, vertex_position);
	} else if(positioner=="DecayTrapFiducial") {
		// uniform fill of 5.5cm radius volume from which events reach detectors
		randomTubePosition(G4ThreeVector(0,0,0), 5.5*cm, 1.5*m, vertex_position);
	} else if(positioner=="SpectrometerVolumeUniform") {
		// uniform fill of (visible) spectrometer volume, for Xe
		randomTubePosition(G4ThreeVector(0,0,0), 7.5*cm, 2.17*m, vertex_position);
	} else if(positioner=="UniformRadialGasFill") {
		randomUniformRadialBins(G4ThreeVector(0,0,0),7.5*cm, 2.15*m, vertex_position);
	} else {
		G4cout << "********* WARNING: Undefined positioner type! Defaulting to 'Fixed'! **********" << G4endl;
	}
	particleGun->SetParticlePosition(vertex_position);
	
	if(gunType=="Cd113m") {
		Cd113mSourceGenerator(anEvent);
	} else if (gunType=="In114" || gunType=="In114E" || gunType=="In114W") {
		In114SourceGenerator(anEvent);
	} else if (gunType=="endpoint" || gunType=="neutronBetaUnpol") {
		double eOrig = particleGun->GetParticleEnergy();
		if(gunType=="neutronBetaUnpol") {
			funcBetaSpectrum->FixParameter(0,neutronBetaEp);
		} else {
			// assume endpoint is set by particle gun energy if not for neutrons
			funcBetaSpectrum->FixParameter(0,eOrig/keV);
		}
		particleGun->SetParticleEnergy(funcBetaSpectrum->GetRandom()*keV);
		G4ThreeVector direction;
		RandomizeMomentum(direction);
		particleGun->SetParticleMomentumDirection(direction);
		displayGunStatus();
		particleGun->GeneratePrimaryVertex(anEvent);
		particleGun->SetParticleEnergy(eOrig);
	} else if(gunType=="uniformRandMomentum") {
		double eOrig = particleGun->GetParticleEnergy();
		particleGun->SetParticleEnergy(G4UniformRand()*eOrig);
		G4ThreeVector direction;
		RandomizeMomentum(direction);
		particleGun->SetParticleMomentumDirection(direction);
		displayGunStatus();
		particleGun->GeneratePrimaryVertex(anEvent);
		particleGun->SetParticleEnergy(eOrig);
	} else if (gunType=="eGunRandMomentum") {
		//for both eGun generator, assume gun energy is set by the
		//standard gun energy command!!!
		G4ThreeVector direction;
		RandomizeMomentum(direction);
		particleGun->SetParticleMomentumDirection(direction);
		displayGunStatus();
		particleGun->GeneratePrimaryVertex(anEvent);
	} else if (gunType=="eGun") {
		//for both eGun generator, assume gun energy is set by the   
		//standard gun energy command!!!
		particleGun->SetParticleMomentumDirection(G4ThreeVector(0,0,1));
		displayGunStatus();
		particleGun->GeneratePrimaryVertex(anEvent);
	} else if(NDL.hasGenerator(gunType)) {
		NucDecaySystem& NDS = NDL.getGenerator(gunType);
		while(!evts.size())
			NDS.genDecayChain(evts);
	} else {
		G4cout << "********* WARNING: Undefined generator type! No events generated!! **********" << G4endl;
	}
	
	throwEvents(evts,anEvent);
	
	//now fill the primary tree here
	gbmAnalysisManager->FillPrimaryData(anEvent,myseed);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

