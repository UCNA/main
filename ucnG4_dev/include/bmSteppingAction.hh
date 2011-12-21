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
// $Id: bmSteppingAction.hh,v 1.5 2011-10-01 18:48:38 mmendenhall Exp $
// GEANT4 tag $Name:  $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef bmSteppingAction_h
#define bmSteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/// user stepping action to check for and abort "trapped" events
class bmSteppingAction : public G4UserSteppingAction {
public:
	/// constructor
    bmSteppingAction();
	
	/// custom per-step action: checks computation time not exceeded
    void UserSteppingAction(const G4Step*);
	
    int GetTrappedFlag() const { return fTrappedFlag; }
    void SetTrappedFlag(int flag) { fTrappedFlag = flag; }
	double GetTimeSpent() const { return timeSpentSoFar; }
	
	/// reset trapping flags
    void Reset(){ fTrappedFlag = timeSpentSoFar = 0.; }
	
private:
    int fTrappedFlag;		//< whether current event is "trapped"
	double timeSpentSoFar;	//< CPU time spent on current event
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
