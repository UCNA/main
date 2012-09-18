#ifndef bmWireVolumeConstruction_HH
#define bmWireVolumeConstruction_HH 1

#include "bmDetectorConstructionUtils.hh"

/// gas-filled region containing anode, cathode planes
class bmWireVolumeConstruction: public MaterialUser {
public:
	/// constructor
	bmWireVolumeConstruction(): anode_R(5*um), cathode_R(25*um), plating_thick(0.2*um), spacing(2.54*mm), NbOfWires(64), planeSpacing(1*cm) { }
	
	/// get width
	G4double GetWidth() const { return 2*cm; }
	
	G4Material* fMWPCGas;				//< MWPC fill gas
	G4LogicalVolume* gas_log;			//< constructed logical volume containing wireplanes
	G4LogicalVolume* cathSeg_log;		//< cathode "segment" containing one wire in gas
	G4LogicalVolume* anodeSeg_log;		//< anode "segment" containing one wire in gas
	G4LogicalVolume* cathode_wire_log;	//< cathode wires logical volume
	G4LogicalVolume* cath_plate_log;	//< gold plating on cathode segment
	G4LogicalVolume* anode_wire_log;	//< anode wires logical volume
	
	/// construct logical container volume
	void Construct(Side s);
	
	const G4double anode_R;				//< anode wire radius
	const G4double cathode_R;			//< cathode wire radius
	const G4double plating_thick;		//< thickness of gold plating on wires
	const G4double spacing;				//< wire spacing
	const G4int NbOfWires;				//< number of wires
	const G4double planeSpacing;		//< spacing between wireplanes
};
	
#endif
