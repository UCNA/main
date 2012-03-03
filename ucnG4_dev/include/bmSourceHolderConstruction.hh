#ifndef bmSourceHolderConstruction_HH
#define bmSourceHolderConstruction_HH 1

#include "bmDetectorConstructionUtils.hh"

/// gas-filled region containing anode, cathode planes
class bmSourceHolderConstruction: public MaterialUser {
public:
	/// constructor
	bmSourceHolderConstruction(): fWindowThick(3.6*um), fCoatingThick(0.1*um),
	fWindowMat(Mylar), fCoatingMat(Al), fSourceHolderThickness(3./16.*inch) { }
	
	/// get thickness
	G4double getHolderThick() const { return fSourceHolderThickness; }
	
	G4double fWindowThick;		//< source foil window single-side thickness
	G4double fCoatingThick;		//< source foil coating thickness
	G4Material* fWindowMat;		//< source foil window material
	G4Material* fCoatingMat;	//< source foil coating material
	
	G4LogicalVolume* container_log;
	G4LogicalVolume* window_log;
	G4LogicalVolume* coating_log[2];
	
	/// construct holder logical volume
	void Construct();
	
protected:
	G4VPhysicalVolume* window_phys;
	G4VPhysicalVolume* coating_phys[2];
	G4double fSourceHolderThickness;
};

#endif