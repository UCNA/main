#ifndef bmDetectorConstructionUtils_HH
#define bmDetectorConstructionUtils_HH 1

#include "globals.hh"
#include "Enums.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4RotationMatrix.hh"
#include "G4VisAttributes.hh"

const G4double inch = 2.54*cm;
const G4double torr = atmosphere/760.;

/// class for building with common materials
class MaterialUser {
public:
	/// constructor
	MaterialUser();
	
	/// set vacuum pressure
	static void setVacuumPressure(G4double pressure);
	
	//static G4UnitDefinition torr_def;
	
	static G4Material* Be;				//< Beryllium for trap windows
	static G4Material* Al;				//< Aluminum
	static G4Material* Si;				//< Silicon
	static G4Material* Cu;				//< Copper for decay trap
	static G4Material* Wu;				//< Tungsten for anode wires
	static G4Material* Au;				//< Gold for cathode wires coating
	
	static G4Material* Vacuum;			//< our slightly crappy vacuum
	static G4Material* Brass;			//< brass for source holder
	static G4Material* SS304;			//< 304 Stainless Steel
	static G4Material* Kevlar;			//< kevlar for wirechamber window support strings
	static G4Material* Mylar;			//< mylar for windows
	static G4Material* Polyethylene;	//< poly for collimator
	static G4Material* WCPentane;		//< Wirechamber fill: (neo)pentane @ 100torr
	static G4Material* WCNitrogen;		//< Wirechamber fill: Nitrogen @ 100torr
	static G4Material* Sci;				//< scintillator material
};


#endif
