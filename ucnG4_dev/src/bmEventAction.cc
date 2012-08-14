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
// $Id: bmEventAction.cc,v 1.5 2011-10-02 15:30:16 mmendenhall Exp $
// GEANT4 tag $Name:  $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "bmEventAction.hh"
#include "bmAnalysisManager.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4ios.hh"
#include "bmSteppingAction.hh"
#include <cassert>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void bmEventAction::BeginOfEventAction(const G4Event* evt) {
	timer.Start();	
	G4cout<<"Beginning of event "<<evt->GetEventID()<<G4endl;
	assert(!system("date"));
	((bmSteppingAction*)fpEventManager->GetUserSteppingAction())->Reset();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void bmEventAction::EndOfEventAction(const G4Event* evt) {
	timer.Stop();  
	if(evt->IsAborted())
		G4cout << "** Event aborted. **" << G4endl;
	G4cout<<"End of event "<<evt->GetEventID()<<G4endl;
	gbmAnalysisManager->FillTrackerData(evt);
	gbmAnalysisManager->FillEventTree();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double bmEventAction::getCPUTime() {
	timer.Stop();
	double t = timer.CpuTime();
	timer.Continue();
	return t;
}

