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
	
	/// set kb
	void SetKb(double c) { kb = c; }
	/// sed density
	void SetRho(double c) { rho = c; }
	
	/// calculate quenching factor for electron at given energy
	double quenchFactor(double E) const;
	
private:
	G4double kb;									///< Birk's law quenching constant
	G4double rho;									///< material density
	std::map<const G4Track*,double> originEnergy;	///< energy at track origin, for Equenched calculaiton
	std::map<G4int,TrackerHit*> tracks;				///< event tracks listed by ID			
	TrackerHitsCollection* trackerCollection;		///< hits objects
};


#include <G4UImessenger.hh>
#include <G4UIcmdWithADouble.hh>

/// UI for TrackerSD
class TrackerSDMessenger: public G4UImessenger {
public:
	/// constructor
    TrackerSDMessenger(TrackerSD*);
	/// destructor
	~TrackerSDMessenger();
    
	/// receive command
    void SetNewValue(G4UIcommand*, G4String);
    
private:
    TrackerSD* mySD;			///< TrackerSD being controlled
    G4UIdirectory*		sdDir;	///< '/SD/<name>/' commands directory
	G4UIcmdWithADouble*	kbCmd;	///< Birk's Law coefficient command
};


#endif
