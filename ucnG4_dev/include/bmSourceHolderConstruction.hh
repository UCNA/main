#ifndef bmSourceHolderConstruction_HH
#define bmSourceHolderConstruction_HH 1

#include "bmDetectorConstructionUtils.hh"

/// gas-filled region containing anode, cathode planes
class bmSourceHolderConstruction: public MaterialUser {
public:
	/// constructor
	bmSourceHolderConstruction(): fWindowThick(3.6*um), fCoatingThick(0.1*um) { }
	
	G4double fWindowThick;
	G4double fCoatingThick;
	
	G4LogicalVolume* container_log;
	G4LogicalVolume* window_log;
	G4LogicalVolume* coating_log[2];
	
	/// construct holder logical volume
	void Construct();
	
protected:
	G4VPhysicalVolume* window_phys;
	G4VPhysicalVolume* coating_phys[2];
};

#endif