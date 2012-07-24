#ifndef bmDetectorPackageConstruction_HH
#define bmDetectorPackageConstruction_HH 1

#include "bmDetectorConstructionUtils.hh"
#include "bmWirechamberConstruction.hh"
#include "bmScintillatorConstruction.hh"

class bmDetectorPackageConstruction: public MaterialUser {
public:
	/// constructor
	bmDetectorPackageConstruction() {}
	/// get scintillator face position within unit
	G4double getScintFacePos() const { return 0; }
	
	bmScintillatorConstruction scint;	//< scintillator assembly
	bmWirechamberConstruction mwpc;		//< wirechamber assembly
	G4LogicalVolume* container_log;		//< overall positioning container
	G4LogicalVolume* mwpc_entrance_log;	//< entrance port container
	G4LogicalVolume* entrance_front_log;//< entrance port front plate
	G4LogicalVolume* entrance_mid_log;	//< entrance port tube
	G4LogicalVolume* entrance_back_log;	//< entrance port back plate (MWPC box cover)
	G4LogicalVolume* mwpc_exit_log;		//< aluminum exit window from wirechamber
	G4LogicalVolume* backstuff_log;		//< miscellaneous mass behind detectors
	G4double entrance_face_pos;
	
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
	G4VPhysicalVolume* backstuff_phys;
};


#endif
