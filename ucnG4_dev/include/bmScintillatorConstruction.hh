#ifndef bmScintillatorConstruction_HH
#define bmScintillatorConstruction_HH 1

#include "bmDetectorConstructionUtils.hh"

/// nitrogen volume with scintillators
class bmScintillatorConstruction: public MaterialUser {
public:
	/// constructor
	bmScintillatorConstruction(): scint_Radius(7.5*cm), backing_Radius(10*cm),
	scint_thick(3.5*mm), dead_thick(3.0*um), backing_thick(1.*inch), lightguide_thick(1.0*cm) { }
	
	/// get position of scintillator face
	G4double getScintFacePos() const { return scintFacePos; }
	/// get width of container
	G4double GetWidth() const { return N2_volume_Z; }
	
	G4double scint_Radius;				//< scintillator disc radius
	G4double backing_Radius;			//< backing veto (and overall volume) radius
	G4double scint_thick;				//< scintillator disc thickness
	G4double dead_thick;				//< dead scintillator thickness, 3 um according to Junhua's thesis
	G4double backing_thick;				//< backing vecto thickness (guess)
	G4double lightguide_thick;			//< light guide thickness at scintillator edge (guess), sets scintillator to backing distance

	G4LogicalVolume* container_log;		//< overall container (nitrogen volume)
	G4LogicalVolume* Dscint_log;		//< scintillator dead layer logical volume
	G4LogicalVolume* scint_log;			//< scintillator logical volume
	G4LogicalVolume* backing_log;		//< backing veto logical volume
	G4LogicalVolume* lightguide_log;	//< lightguide material logical volume
	
	/// construct logical container volume
	void Construct(Side s);
	
protected:
	G4double scintFacePos;				//< position of scintillator face in container
	G4double N2_volume_Z;				//< z width of N2 volume
	G4VPhysicalVolume* Dscint_phys;		//< dead layer physical volume
	G4VPhysicalVolume* scint_phys;		//< scintillator physical volume
	G4VPhysicalVolume* backing_phys;	//< backing veto physical volume
	G4VPhysicalVolume* lightguide_phys;	//< lightguide material physical volume
};

#endif
