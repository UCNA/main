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
#include "bmDetectorConstruction.hh"
#include "bmPrimaryGeneratorMessenger.hh"
#include <vector>
#include <Rtypes.h>
#include <TF1.h>

using namespace std;
class bmPrimaryGeneratorMessenger;

/// (uncorrected) beta spectrum probability with arbitrary endpoint (for use with TF1)
double genericBetaSpectrum(double* x, double *par);
/// (uncorrected) 2D polarized beta spectrum as a function of kinetic energy and cos theta (for use with TF2)
double asymNeutronBetaSpectrum(double* x, double*);


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class bmPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
public:
	/// constructor
    bmPrimaryGeneratorAction(bmDetectorConstruction*);  
	/// destructor
	~bmPrimaryGeneratorAction();
	
	void GeneratePrimaries(G4Event*);
	void SetGunType(G4String val) { gunType = val; }
	void SetPositioner(G4String val) { positioner = val; }
	void SetSourceRadius(double r) { sourceRadius = r; }
	
private:
	G4ParticleGun* particleGun;
	bmDetectorConstruction* myDetector;
	bmPrimaryGeneratorMessenger* gunMessenger;	// messenger of this class
	G4String gunType;		// event generator gun to use
	G4String positioner;	// how to position initial events
	double sourceRadius;	// radius for sealed source generator
	
	/// polarized neutron beta decay with given flipper state
	void polarizedNeutronBetaDecayGenerator(G4Event* anEvent, bool flipper);
	
	/// throw multiple electrons and gammas in one event
	void throwElectronsAndGammas(const std::vector<G4double>& electrons,
								 const std::vector<G4double>& gammas,
								 G4Event* anEvent);

	/// print what the particle gun is set up to do
	void displayGunStatus();
	
	void Sn113SourceGenerator(G4Event* anEvent);
	void Sn113GSourceGenerator(G4Event* anEvent);
	void Bi207SourceGenerator(G4Event* anEvent);
	void Bi207GSourceGenerator(G4Event* anEvent);
	void Ce139SourceGenerator(G4Event* anEvent);
	void Ce139GSourceGenerator(G4Event* anEvent);
	void Sr85SourceGenerator(G4Event* anEvent);
	void Cd109SourceGenerator(G4Event* anEvent);
	void Cd109GSourceGenerator(G4Event* anEvent);
	void Cd113mGSourceGenerator(G4Event* anEvent);
	void Mn54SourceGenerator(G4Event* anEvent);
	
	// Indium source and backgrounds
	void In114SourceGenerator(G4Event* anEvent);
	void In114GSourceGenerator(G4Event* anEvent);
	void Sc46GSourceGenerator(G4Event* anEvent);
	void Co60GSourceGenerator(G4Event* anEvent);
	void Ag110GSourceGenerator(G4Event* anEvent);
	
	// Xenon isotopes
	void Xe125_1_2p_SourceGenerator(G4Event* anEvent);
	void Xe133_11_2m_SourceGenerator(G4Event* anEvent);
	void Xe133_3_2p_SourceGenerator(G4Event* anEvent);
	void Xe135_11_2m_SourceGenerator(G4Event* anEvent);
	void Xe135_3_2p_SourceGenerator(G4Event* anEvent);
	void Xe137_7_2m_SourceGenerator(G4Event* anEvent);
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


