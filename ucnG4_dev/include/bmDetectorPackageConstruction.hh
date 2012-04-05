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
	
	bmScintillatorConstruction scint;
	bmWirechamberConstruction mwpc;
	G4LogicalVolume* container_log;
	G4LogicalVolume* mwpc_entrance_log;
	G4LogicalVolume* mwpc_exit_log;
	
	/// construct logical container volume
	void Construct(Side s);
	
protected:
	G4VPhysicalVolume* scint_phys;
	G4VPhysicalVolume* mwpc_phys;
	G4VPhysicalVolume* mwpc_entrance_phys;
	G4VPhysicalVolume* mwpc_exit_phys;
};


#endif
