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
// $Id: bmPrimaryGeneratorMessenger.cc,v 1.4 2011-09-14 00:23:06 mmendenhall Exp $
// GEANT4 tag $Name:  $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "bmPrimaryGeneratorMessenger.hh"

#include "bmPrimaryGeneratorAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

bmPrimaryGeneratorMessenger::bmPrimaryGeneratorMessenger(bmPrimaryGeneratorAction* bmGun)
:bmAction(bmGun) {
	gunDir = new G4UIdirectory("/benchmark/gun/");
	gunDir->SetGuidance("PrimaryGenerator control");
	
	gunTypeCmd = new G4UIcmdWithAString("/benchmark/gun/type",this);
	gunTypeCmd->SetGuidance("Set the generator gun type.");
	gunTypeCmd->SetGuidance(" Choices: eGun eGunRandMomentum");
	gunTypeCmd->SetDefaultValue("eGun");
	gunTypeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
	
	gunPtclCmd = new G4UIcmdWithAString("/benchmark/gun/particle",this);
	gunPtclCmd->SetGuidance("Set the gun particle thrown.");
	gunPtclCmd->SetDefaultValue("e-");
	gunPtclCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
	
	positionerCmd = new G4UIcmdWithAString("/benchmark/gun/positioner",this);
	positionerCmd->SetGuidance("Set the generator gun positioner.");
	positionerCmd->SetGuidance(" Choice : Fixed, SourceDrop, DecayTrapUniform, DecayTrapFiducial");
	positionerCmd->SetDefaultValue("Fixed");
	positionerCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
	
	sourceRadiusCmd = new G4UIcmdWithADoubleAndUnit("/benchmark/gun/sourceRadius",this);
	sourceRadiusCmd->SetGuidance("Radius for SourceDrop generator");
	sourceRadiusCmd->SetDefaultValue(1.5*mm);
	sourceRadiusCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

bmPrimaryGeneratorMessenger::~bmPrimaryGeneratorMessenger()
{
	delete gunTypeCmd;
	delete positionerCmd;
	delete sourceRadiusCmd;
	delete gunDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void bmPrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command, G4String newValue) { 
	if( command == gunTypeCmd )
		bmAction->SetGunType(newValue);
	if( command == gunPtclCmd )
		bmAction->SetParticleType(newValue);
	if( command == positionerCmd )
		bmAction->SetPositioner(newValue);
	if( command == sourceRadiusCmd )
		bmAction->SetSourceRadius(sourceRadiusCmd->GetNewDoubleValue(newValue));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

