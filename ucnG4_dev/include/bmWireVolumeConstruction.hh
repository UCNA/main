#ifndef bmWireVolumeConstruction_HH
#define bmWireVolumeConstruction_HH 1

#include "bmDetectorConstructionUtils.hh"

/// gas-filled region containing anode, cathode planes
class bmWireVolumeConstruction: public MaterialUser {
public:
	/// constructor
	bmWireVolumeConstruction() { }
	
	/// get width
	G4double GetWidth() const { return 2*cm; }
	
	G4Material* fMWPCGas;				//< MWPC fill gas
	G4LogicalVolume* gas_log;			//< constructed logical volume containing wireplanes
	G4LogicalVolume* cathSeg_log;		//< cathode "segment" containing one wire in gas
	G4LogicalVolume* anodeSeg_log;		//< anode "segment" containing one wire in gas
	G4LogicalVolume* cathode_wire_log;	//< cathode wires logical volume
	G4LogicalVolume* anode_wire_log;	//< anode wires logical volume
	
	/// construct logical container volume
	void Construct(Side s);
	
protected:
	
};
	
#endif
