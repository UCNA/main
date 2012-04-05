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

	/// construct logical container volume
	void Construct(Side s);
	
protected:
	
};
	
#endif
