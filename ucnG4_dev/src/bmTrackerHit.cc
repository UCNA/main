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
// $Id: bmTrackerHit.cc,v 1.8 2011-10-01 16:09:56 mmendenhall Exp $
// GEANT4 tag $Name:  $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "bmTrackerHit.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

G4Allocator<bmTrackerHit> bmTrackerHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

bmTrackerHit::bmTrackerHit(): eDepSoFar(0), eDepQuenchedSoFar(0), hitPosition(), edepWeightedPosition(),
edepWeightedPosition2(), incidentMomentum(), exitMomentum(), vertex() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void bmTrackerHit::Print()
{
  G4cout << "  trackID: " << trackID 
	 << "  vertex: "<< G4BestUnit(vertex,"Length")
	 << "  created in "<<creatorVolumeName
	 << "  in "<<volumeName
         << "  incident energy " << G4BestUnit(incidentEnergy,"Energy")
	 << "  position: " << G4BestUnit(hitPosition,"Length") 
	 << "  time: "<<G4BestUnit(hitTime,"Time")
	 << "  edep: "<<G4BestUnit(eDepSoFar,"Energy")
	 << "  edep quenched: "<<G4BestUnit(eDepQuenchedSoFar,"Energy")
	 << "  incident momentm: "<<G4BestUnit(incidentMomentum,"Energy")
	 << "  exit momentum "<<G4BestUnit(exitMomentum,"Energy")
	 << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void bmTrackerHit::fillTrackInfo(bmTrackInfo& h) const {
	h.trackID = GetTrackID();
	h.hitTime = GetHitTime()/ns;
	h.KE = GetIncidentEnergy()/keV;
	h.Edep = GetEdep()/keV; 
	h.EdepQuenched = GetEdepQuenched()/keV; 
	
	for(unsigned int i=0; i<3; i++) {
		h.edepPos[i] = GetEdepPos()[i]/(cm*keV);
		h.edepPos2[i] = GetEdepPos2()[i]/(cm*cm*keV);
		h.vertexPos[i] = GetVertex()[i]/cm;
		h.inPos[i] = GetPos()[i]/cm;
		h.pIn[i] = GetIncidentMomentum()[i]/keV;
		h.pOut[i] = GetExitMomentum()[i]/keV;
	}
	
	h.pID = GetPID();
}
