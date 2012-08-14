#ifndef bmDecayTrapConstruction_HH
#define bmDecayTrapConstruction_HH 1

#include "bmDetectorConstructionUtils.hh"
#include "G4LogicalVolume.hh"
#include "WiggleSheet.hh"

/// gas-filled region containing anode, cathode planes
class bmDecayTrapConstruction: public MaterialUser {
public:
	/// constructor
	bmDecayTrapConstruction(): fWindowThick(0.7*um), fCoatingThick(0.3*um), fIRtrap(2.45*inch), decayTube_Wall(2*mm),
	fIRcollimator(2.3*inch), fTubeMat(Cu), fCollimatorMat(Polyethylene), fWindowMat(Mylar), fCoatingMat(Be) { }
	
	G4double fWindowThick;		//< window thickness
	G4double fCoatingThick;		//< window Be coating thickness
	G4double fIRtrap;			//< decay trap IR
	G4double decayTube_Wall;	//< decay trap wall thickness
	G4double fIRcollimator;		//< collimator IR
	G4Material* fTubeMat;		//< decay tube material
	G4Material* fCollimatorMat;	//< collimator material
	G4Material* fWindowMat;		//< decay tube window material
	G4Material* fCoatingMat;	//< decay tube coating material

	G4LogicalVolume* decayTube_log;			//< decay trap tube
	G4LogicalVolume* trap_win_log[2];		//< trap window volume
	G4LogicalVolume* mylar_win_log[2];		//< mylar layer of window
	G4LogicalVolume* be_win_log[2];			//< berillium layer of window
	G4LogicalVolume* trap_monitor_log[2];	//< extra event monitoring region
	G4LogicalVolume* collimator_log[2];		//< collimator
	G4LogicalVolume* collimatorBack_log[2];	//< bracket behind collimator
	G4LogicalVolume* plug_log;				//< decay trap plug edge discontinuity
	
	WiggleSheet wigglefoils[2];				//< optional replacement crinkly endcap foils
	
	/// construct in given world volume
	void Construct(G4LogicalVolume* world, double crinkleAngle = 0);
};

#endif
