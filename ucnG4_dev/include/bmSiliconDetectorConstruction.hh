#ifndef bmSiliconDetectorConstruction_HH
#define bmSiliconDetectorConstruction_HH 1

#include "bmDetectorConstructionUtils.hh"

/// model for circular silicon detector surrounded by holder, possibly with backing
class bmSiliconDetectorConstruction: public MaterialUser {
public:
	/// constructor
	bmSiliconDetectorConstruction(): fDetThick(1.*mm), fDetRadius(5.*mm), pDetMat(Si),
	fHolderThick(10.*mm), fHolderRadius(10.*mm), pHolderMat(Al),
	fBackingThick(4.5*mm), pBackingMat(Vacuum) {}
	
	G4double fDetThick;				//< active detector thickness
	G4double fDetRadius;			//< active detector radius
	G4Material* pDetMat;			//< active detector material
	
	G4double fHolderThick;			//< thickness of holder ring
	G4double fHolderRadius;			//< radius of holder ring
	G4Material* pHolderMat;			//< holder material
	
	G4double fBackingThick;			//< backing thickness
	G4Material* pBackingMat;		//< backing material
		
	G4LogicalVolume* container_log;	//< overall container
	G4LogicalVolume* holder_log;	//< outer holding ring
	G4LogicalVolume* det_log;		//< active detector volume
	G4LogicalVolume* backing_log;	//< backing
	
	/// construct logical container volume
	void Construct();
	
protected:
	G4VPhysicalVolume* det_phys;		//< placed active detector
	G4VPhysicalVolume* holder_phys;		//< placed holder
	G4VPhysicalVolume* backing_phys;	//< placed backing
};

#endif
