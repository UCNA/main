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
// $Id: bmTrackerSD.hh,v 1.4 2011-10-03 17:52:31 mmendenhall Exp $
// GEANT4 tag $Name:  $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef bmTrackerSD_h
#define bmTrackerSD_h 1

#include "G4VSensitiveDetector.hh"
#include "bmTrackerHit.hh"
#include <map>

class G4Step;
class G4HCofThisEvent;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/// stores/evaluates tracks passing through each SD
class bmTrackerSD : public G4VSensitiveDetector {
public:
	/// constructor
	bmTrackerSD(G4String);
	
	/// run at start of each event
	void Initialize(G4HCofThisEvent*);
	/// run for each step during event
	G4bool ProcessHits(G4Step* aStep, G4TouchableHistory*);
	/// run at end of event
	void EndOfEvent(G4HCofThisEvent*);
	
private:
	std::map<const G4Track*,double> originEnergy;		//< energy at track origin, for Equenched calculaiton
	std::map<G4int,bmTrackerHit*> tracks;				//< event tracks listed by ID			
	bmTrackerHitsCollection* trackerCollection;			//< hits objects
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

