#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h

#include <Rtypes.h>
#include <TF1.h>
#include <vector>

#include "PrimaryGeneratorMessenger.hh"
#include "SurfaceGenerator.hh"
#include "ElectronBindingEnergy.hh"
#include "NuclEvtGen.hh"
#include "DetectorConstruction.hh"

#include <G4VUserPrimaryGeneratorAction.hh>
#include <G4String.hh>
#include <G4ThreeVector.hh>
#include <G4ParticleGun.hh>
#include <G4Event.hh>
#include <G4VUserEventInformation.hh>

/// User event information for recording primary event weighting
class PrimEvtWeighting: public G4VUserEventInformation {
public:
	/// constructor
	PrimEvtWeighting(double W): w(W) {}
	/// print info
	void Print() const { G4cout << "Primary weighting " << w <<  G4endl; }
	
	double w;	///< event primary weight
};

using namespace std;
class PrimaryGeneratorMessenger;

/// (uncorrected) beta spectrum probability with arbitrary endpoint (for use with TF1)
double genericBetaSpectrum(double* x, double *par);
	
class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
public:
	/// constructor
    PrimaryGeneratorAction(DetectorConstruction*);  
	/// destructor
	~PrimaryGeneratorAction();
	
	/// generate primaries for next event
	void GeneratePrimaries(G4Event*);
	
	/// set input events tree
	void SetEventFile(G4String val);
	/// set radius for "source drop" positioner
	void SetSourceRadius(double r) { sourceRadius = r; }
	/// set positioning relative to holder
	void SetPosRelHolder(bool b) { relToSourceHolder = b; }
	/// set positioning offset
	void SetPosOffset(const G4ThreeVector& v) { posOffset = v; }
	
private:
	G4ParticleGun* particleGun;				///< particle gun primary event thrower
	DetectorConstruction* myDetector;		///< detector geometry
	PrimaryGeneratorMessenger* myMessenger;	///< UI messenger of this class
	EventTreeScanner* ETS;					///< Reader for input primary events
	
	
	/// set vertex positions for each primary
	void setVertices(std::vector<NucDecayEvent>& evts);
	G4ThreeVector posOffset;				///< base positioning offset
	double sourceRadius;					///< spread radius for source droplets
	bool relToSourceHolder;					///< make positions relative to source holder, instead of geometry origin
	
	/// throw a cluster of events
	void throwEvents(const std::vector<NucDecayEvent>& evts, G4Event* anEvent);

	/// print what the particle gun is set up to do
	void displayGunStatus();
	
	/// initialize random seed for event
	void initEventRandomSeed(G4Event* anEvent);
	long myseed;							///< random seed for event
};

#endif


