#ifndef bmDecayTrapConstruction_HH
#define bmDecayTrapConstruction_HH 1

#include "bmDetectorConstructionUtils.hh"
#include "G4LogicalVolume.hh"
#include "WiggleSheet.hh"

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

	G4LogicalVolume* decayTube_log;			//< decay trap tube
	G4LogicalVolume* trap_win_log[2];		//< trap window volume
	G4LogicalVolume* mylar_win_log[2];		//< mylar layer of window
	G4LogicalVolume* be_win_log[2];			//< berillium layer of window
	G4LogicalVolume* trap_monitor_log[2];	//< extra event monitoring region
	
	WiggleSheet wigglefoils[2];				//< optional replacement crinkly endcap foils
	
	/// construct in given world volume
	void Construct(G4LogicalVolume* world, double crinkleAngle = 0);
};

#endif
