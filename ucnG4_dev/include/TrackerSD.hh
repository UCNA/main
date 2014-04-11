#ifndef TrackerSD_h
#define TrackerSD_h

#include "TrackerHit.hh"
#include <map>

#include <G4VSensitiveDetector.hh>
#include <G4HCofThisEvent.hh>
#include <G4Step.hh>

/// stores/evaluates tracks passing through each SD
class TrackerSD : public G4VSensitiveDetector {
public:
	/// constructor
	TrackerSD(G4String);
	
	/// run at start of each event
	void Initialize(G4HCofThisEvent*);
	/// run for each step during event
	G4bool ProcessHits(G4Step* aStep, G4TouchableHistory*);
	/// run at end of event
	void EndOfEvent(G4HCofThisEvent*);
	
private:
	std::map<const G4Track*,double> originEnergy;	///< energy at track origin, for Equenched calculaiton
	std::map<G4int,TrackerHit*> tracks;				///< event tracks listed by ID			
	TrackerHitsCollection* trackerCollection;		///< hits objects
};

#endif
