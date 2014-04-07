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
// $Id: PrimaryGeneratorMessenger.hh,v 1.4 2011-09-28 01:29:04 mmendenhall Exp $
// GEANT4 tag $Name:  $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef PrimaryGeneratorMessenger_h
#define PrimaryGeneratorMessenger_h

#include "G4UImessenger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIdirectory.hh"
#include "PrimaryGeneratorAction.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PrimaryGeneratorAction;

/// messenger UI for primary event generator
class PrimaryGeneratorMessenger: public G4UImessenger {
public:
	/// constructor
    PrimaryGeneratorMessenger(PrimaryGeneratorAction*);
	/// destructor
	~PrimaryGeneratorMessenger();
    
	/// receive command
    void SetNewValue(G4UIcommand*, G4String);
    
private:
    PrimaryGeneratorAction* Action;				//< primary generator action being controlled
    G4UIdirectory*			gunDir;					//< '/gun/' commands directory
    G4UIcmdWithAString*		gunTypeCmd;				//< control momentum generator (beam, random, file-specified)
	G4UIcmdWithAString*		gunPtclCmd;				//< control default particle to throw
	G4UIcmdWithAString*		positionerCmd;			//< control event positioner to use
	G4UIcmdWithABool*		srcrelCmd;				//< control positions relative to source holder
	G4UIcmdWithAString*		eventFileCmd;			//< control file to read events from
	G4UIcmdWithADoubleAndUnit* sourceRadiusCmd;		//< control radius of "source drop" positioner
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

