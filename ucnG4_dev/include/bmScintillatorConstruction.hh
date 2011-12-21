#ifndef bmScintillatorConstruction_HH
#define bmScintillatorConstruction_HH 1

#include "bmDetectorConstructionUtils.hh"

/// nitrogen volume with scintillators
class bmScintillatorConstruction: public MaterialUser {
public:
	/// constructor
	bmScintillatorConstruction() { }
	
	/// get position of scintillator face
	G4double getScintFacePos() const { return scintFacePos; }
	/// get width of container
	G4double GetWidth() const { return N2_volume_Z; }
	
	G4double mwpc_entrance_R;
	G4double mwpc_exit_R;
	G4LogicalVolume* container_log;
	G4LogicalVolume* Dscint_log;
	G4LogicalVolume* scint_log;
	G4LogicalVolume* backing_log;
	
	/// construct logical container volume
	void Construct(Side s);
	
protected:
	G4double scintFacePos;
	G4double N2_volume_Z;
	G4VPhysicalVolume* Dscint_phys;
	G4VPhysicalVolume* scint_phys;
	G4VPhysicalVolume* backing_phys;
};

#endif
