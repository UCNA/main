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
// $Id: bmPrimaryGeneratorAction.hh,v 1.14 2011-12-12 01:29:36 mmendenhall Exp $
// GEANT4 tag $Name:  $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef bmPrimaryGeneratorAction_h
#define bmPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4String.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleGun.hh"
#include "G4Event.hh"
#include "G4VUserEventInformation.hh"
#include "bmDetectorConstruction.hh"
#include "bmPrimaryGeneratorMessenger.hh"
#include "SurfaceGenerator.hh"
#include "ElectronBindingEnergy.hh"
#include "NuclEvtGen.hh"
#include <vector>
#include <Rtypes.h>
#include <TF1.h>

/// User event information for recording primary event weighting
class PrimEvtWeighting: public G4VUserEventInformation {
public:
	/// constructor
	PrimEvtWeighting(double W): w(W) {}
	/// print info
	void Print() const { G4cout << "Primary weighting " << w <<  G4endl; }
	
	double w;	//< event primary weight
};

using namespace std;
class bmPrimaryGeneratorMessenger;

/// (uncorrected) beta spectrum probability with arbitrary endpoint (for use with TF1)
double genericBetaSpectrum(double* x, double *par);
/// corrected beta decay spectrum from heavy nucleus (for use with TF1)
double heavyBetaSpectrum(double* x, double* par);
	
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class bmPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
public:
	/// constructor
    bmPrimaryGeneratorAction(bmDetectorConstruction*);  
	/// destructor
	~bmPrimaryGeneratorAction();
	
	void GeneratePrimaries(G4Event*);
	void SetGunType(G4String val) { gunType = val; }
	void SetParticleType(G4String val) { particleType = val; }
	void SetPositioner(G4String val) { positioner = val; }
	void SetSourceRadius(double r) { sourceRadius = r; }
	
private:
	G4ParticleGun* particleGun;
	bmDetectorConstruction* myDetector;
	bmPrimaryGeneratorMessenger* gunMessenger;	//< messenger of this class
	G4String gunType;							//< event generator gun to use
	G4String particleType;						//< particle type to throw
	G4String positioner;						//< how to position initial events
	double sourceRadius;						//< radius for sealed source generator
	G4ThreeVector vertex_position;				//< event vertex position
	
	/// choose event vertex position
	void selectVertex();
	
	/// throw multiple electrons and gammas in one event
	void throwElectronsAndGammas(const std::vector<G4double>& electrons,
								 const std::vector<G4double>& gammas,
								 G4Event* anEvent);
	/// throw a gamma towards some point on specified surface, recording weight relative to 4pi uniform
	void throwGammaAt(SurfaceSeg* S, double eGamma, G4Event* anEvent);

	/// throw a cluster of events
	void throwEvents(const std::vector<NucDecayEvent>& evts, G4Event* anEvent);
	
	/// print what the particle gun is set up to do
	void displayGunStatus();
	
	/// Cd113 Metastable 11/2-, 14.1year HL
	void Cd113mSourceGenerator(G4Event* anEvent);
	/// In114, based on NuDat 2.6
	void In114SourceGenerator(G4Event* anEvent);
	/// approximation for neutron capture on Cu gammas, based on probabilities in Robby's eLog 134
	void nCaptureCuGammas(G4Event* anEvent, SurfaceAssembly* S);
	/// approximation for neutron capture on Fe gammas
	void nCaptureFeGammas(G4Event* anEvent, SurfaceAssembly* S);
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


