#ifndef DecayTrapConstruction_HH
#define DecayTrapConstruction_HH

#include "DetectorConstructionUtils.hh"
#include "WiggleSheet.hh"

#include <G4LogicalVolume.hh>

/// gas-filled region containing anode, cathode planes
class DecayTrapConstruction: public MaterialUser {
public:
	/// constructor
	DecayTrapConstruction(): fIRtrap(2.45*inch), decayTube_Wall(2*mm),
				 fIRcollimator(2.3*inch), fColThickness(0.8*inch), fColLength(0.8*inch), fTubeMat(Cu), fCollimatorMat(Delrin), fWindowMat(Mylar), fCoatingMat(Be) { }
	
	G4double fWindowThick[2];		///< window thickness
	G4double fCoatingThick[2];		///< window Be coating thickness
        G4double thicknessOfTrapWindow[2];
        G4double beWinPosZ[2];
        G4double mylarWinPosZ[2];
        G4double trap_winPosZ[2];
	G4double fIRtrap;			///< decay trap IR
	G4double decayTube_Wall;	///< decay trap wall thickness
	G4double fIRcollimator;	        ///< collimator IR
        G4double fColThickness;         ///< from IR to OR of collimator
        G4double fColLength;            //< length in z direction
	G4Material* fTubeMat;		///< decay tube material
	G4Material* fCollimatorMat;	///< collimator material
	G4Material* fWindowMat;		///< decay tube window material
	G4Material* fCoatingMat;	///< decay tube coating material
        G4Tubs* trap_win_tube[2];
        G4Tubs* mylarTube[2];        
        G4Tubs* beTube[2];

	G4LogicalVolume* decayTube_log;			///< decay trap tube
	G4LogicalVolume* trap_win_log[2];		///< trap window volume
	G4LogicalVolume* mylar_win_log[2];		///< mylar layer of window
	G4LogicalVolume* be_win_log[2];			///< berillium layer of window
	G4LogicalVolume* trap_monitor_log[2];	///< extra event monitoring region
	G4LogicalVolume* collimator_log[2];		///< collimator
	G4LogicalVolume* collimatorBack_log[2];	///< bracket behind collimator
	G4LogicalVolume* plug_log;				///< decay trap plug edge discontinuity
	
	WiggleSheet wigglefoils[2];				///< optional replacement crinkly endcap foils
	
	/// construct in given world volume
	void Construct(G4LogicalVolume* world, double crinkleAngle = 0);
};

#endif
