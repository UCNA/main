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
// $Id: bmRunAction.hh,v 1.2 2011-10-01 18:48:37 mmendenhall Exp $
// GEANT4 tag $Name:  $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef bmRunAction_h
#define bmRunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"
#include "G4UImessenger.hh"
#include "G4UIdirectory.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4Run;
class G4UIcommand;

/// user actions for whole run
class bmRunAction : public G4UserRunAction, public G4UImessenger {
public:
	/// constructor
    bmRunAction();
	/// destructor
	~bmRunAction();
	
	/// at begin of run
    void BeginOfRunAction(const G4Run*);
	/// at end of run (write/close files)
    void EndOfRunAction(const G4Run*);
	/// UI controls
	void SetNewValue(G4UIcommand * command,G4String newValue);

protected:
	G4UIdirectory  *fFileDir;		//< UI directory for opening output file
	G4UIcommand    *fFileCmd;		//< UI command for opening output file
	G4UIdirectory  *fRunDir;		//< UI directory for setting run number
	G4UIcommand    *fRunNumberCmd;	//< UI command for setting run number
	G4int          fRunNumber;		//< run number
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif





