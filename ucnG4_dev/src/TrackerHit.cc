#include "Enums.hh"
#include "TrackerHit.hh"

#include <G4SystemOfUnits.hh>
#include <G4UnitsTable.hh>
#include <G4VVisManager.hh>
#include <G4Circle.hh>
#include <G4Colour.hh>
#include <G4VisAttributes.hh>

G4Allocator<TrackerHit> TrackerHitAllocator;

TrackerHit::TrackerHit(): eDepSoFar(0), eDepQuenchedSoFar(0), hitPosition(), edepWeightedPosition(),
edepWeightedPosition2(), incidentMomentum(), exitMomentum(), vertex() {}

void TrackerHit::Print() {
  G4cout	<< "  trackID: " << trackID
			<< "  vertex: " << G4BestUnit(vertex,"Length")
			<< "  created in " << creatorVolumeName
			<< "  in " << volumeName
			<< "  incident energy " << G4BestUnit(incidentEnergy,"Energy")
			<< "  position: " << G4BestUnit(hitPosition,"Length")
			<< "  time: " << G4BestUnit(hitTime,"Time")
			<< "  edep: " << G4BestUnit(eDepSoFar,"Energy")
			<< "  edep quenched: " << G4BestUnit(eDepQuenchedSoFar,"Energy")
			<< "  incident momentm: " << G4BestUnit(incidentMomentum,"Energy")
			<< "  exit momentum " << G4BestUnit(exitMomentum,"Energy")
			<< G4endl;
}

void TrackerHit::fillTrackInfo(TrackInfo& h) const {
	h.trackID = GetTrackID();
	h.hitTime = GetHitTime()/ns;
	h.KE = GetIncidentEnergy()/keV;
	h.Edep = GetEdep()/keV; 
	h.EdepQuenched = GetEdepQuenched()/keV; 
	h.isEntering = (originEnergy==0);
	
	for(AxisDirection d=X_DIRECTION; d<=Z_DIRECTION; ++d) {
		h.edepPos[d] = GetEdepPos()[d]/(cm*keV);
		h.edepPos2[d] = GetEdepPos2()[d]/(cm*cm*keV);
		h.vertexPos[d] = GetVertex()[d]/cm;
		h.inPos[d] = GetPos()[d]/cm;
		h.pIn[d] = GetIncidentMomentum()[d]/keV;
		h.pOut[d] = GetExitMomentum()[d]/keV;
	}
	
	h.pID = GetPID();
}
