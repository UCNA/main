#ifndef bmSourceHolderConstruction_HH
#define bmSourceHolderConstruction_HH 1

#include "bmDetectorConstructionUtils.hh"
#include "G4UImessenger.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

/// gas-filled region containing anode, cathode planes
class bmSourceHolderConstruction: public MaterialUser, G4UImessenger {
public:
	/// constructor
	bmSourceHolderConstruction();
	
	/// get thickness
	G4double getHolderThick() const { return fSourceHolderThickness; }
	
	G4double fWindowThick;		//< source foil window single-side thickness
	G4double fCoatingThick;		//< source foil coating thickness
	G4Material* fWindowMat;		//< source foil window material
	G4Material* fCoatingMat;	//< source foil coating material
	
	G4LogicalVolume* container_log;
	G4LogicalVolume* window_log;
	G4LogicalVolume* coating_log[2];
	
	/// construct holder logical volume
	void Construct();
	
	/// UI interface
	virtual void SetNewValue(G4UIcommand * command,G4String newValue);
	
protected:
	G4VPhysicalVolume* window_phys;
	G4VPhysicalVolume* coating_phys[2];
	G4double fSourceHolderThickness;

private:
	G4UIdirectory* pUIdir;						//< UI Directory for source holder construction
	G4UIcmdWithADoubleAndUnit* pWindowThickCmd;	//< source holder window thickness command
};

#endif