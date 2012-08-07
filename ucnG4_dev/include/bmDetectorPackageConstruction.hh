#ifndef bmDetectorPackageConstruction_HH
#define bmDetectorPackageConstruction_HH 1

#include "bmDetectorConstructionUtils.hh"
#include "bmWirechamberConstruction.hh"
#include "bmScintillatorConstruction.hh"

class bmDetectorPackageConstruction: public MaterialUser {
public:
	/// constructor
	bmDetectorPackageConstruction(): detPackageRadius(6.0*inch),
	mwpc_entrance_thickness(0.375*inch), mwpc_entrance_r(3.0*inch),
	mwpc_entrance_depth(5.0*inch), frontwin_frame_thick(1.0*inch), backwin_frame_thick(0.5*inch) {}
	
	/// get scintillator face position within unit
	G4double getScintFacePos() const { return 0; }
	
	G4double detPackageRadius;
	G4double mwpc_entrance_thickness;	//< MWPC entrance tube wall thickness
	G4double mwpc_entrance_r;			//< MWPC entrance tube radius
	G4double mwpc_entrance_depth;		//< MWPC entrance tube depth
	G4double frontwin_frame_thick;		//< MWPC front window frame thickness
	G4double backwin_frame_thick;		//< MWPC exit window frame thickness
	
	bmScintillatorConstruction scint;	//< scintillator assembly
	bmWirechamberConstruction mwpc;		//< wirechamber assembly
	G4LogicalVolume* container_log;		//< overall positioning container
	G4LogicalVolume* mwpc_entrance_log;	//< entrance port container
	G4LogicalVolume* entrance_front_log;//< entrance port front plate
	G4LogicalVolume* entrance_mid_log;	//< entrance port tube
	G4LogicalVolume* entrance_back_log;	//< entrance port back plate (MWPC box cover)
	G4LogicalVolume* mwpc_exit_log;		//< aluminum exit window from wirechamber
	G4LogicalVolume* mwpc_exit_N2_log;	//< N2 between exit window and scintillator
	G4LogicalVolume* backstuff_log;		//< miscellaneous mass behind detectors
	
	G4double entrance_face_pos;			//< entrance window port entrance relative to scint face
	G4double entrance_win_pos;			//< MWPC entrance window position relative to scint face
	G4double exit_frame_pos;			//< exit window frame pos
	
	/// construct logical container volume
	void Construct(Side s);
	
protected:
	G4VPhysicalVolume* scint_phys;
	G4VPhysicalVolume* mwpc_phys;
	G4VPhysicalVolume* mwpc_entrance_phys;
	G4VPhysicalVolume* entrance_front_phys;
	G4VPhysicalVolume* entrance_mid_phys;
	G4VPhysicalVolume* entrance_back_phys;
	G4VPhysicalVolume* mwpc_exit_phys;
	G4VPhysicalVolume* mwpc_exit_N2_phys;
	G4VPhysicalVolume* backstuff_phys;
};


#endif
