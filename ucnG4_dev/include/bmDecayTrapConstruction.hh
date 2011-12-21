#ifndef bmDecayTrapConstruction_HH
#define bmDecayTrapConstruction_HH 1

#include "bmDetectorConstructionUtils.hh"

/// gas-filled region containing anode, cathode planes
class bmDecayTrapConstruction: public MaterialUser {
public:
	/// constructor
	bmDecayTrapConstruction(): fWindowThick(0.7*um), fCoatingThick(0.3*um),
	fTubeMat(Cu), fWindowMat(Mylar), fCoatingMat(Be) { }
	
	G4double fWindowThick;		//< window thickness
	G4double fCoatingThick;		//< window Be coating thickness
	G4Material* fTubeMat;		//< decay tube material
	G4Material* fWindowMat;		//< decay tube window material
	G4Material* fCoatingMat;	//< decay tube coating material

	G4LogicalVolume* decayTube_log;
	G4VPhysicalVolume* decayTube_phys;
	G4LogicalVolume* trap_win_log[2];
	G4VPhysicalVolume* trap_win_phys[2];
	G4LogicalVolume* mylar_win_log[2];
	G4VPhysicalVolume* mylar_win_phys[2];	
	G4LogicalVolume* be_win_log[2];
	G4VPhysicalVolume* be_win_phys[2];	
	G4LogicalVolume* trap_monitor_log[2];
	G4VPhysicalVolume* trap_monitor_phys[2];
	
	/// construct in given world volume
	void Construct(G4LogicalVolume* world);
};

#endif
