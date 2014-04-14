#ifndef DETECTORCONSTRUCTION_HH
#define DETECTORCONSTRUCTION_HH

#include <Rtypes.h>
#include <TString.h>

#include "TrackerSD.hh"
#include "Field.hh"
#include "DecayTrapConstruction.hh"
#include "SourceHolderConstruction.hh"
#include "DetectorPackageConstruction.hh"
#include "SiliconDetectorConstruction.hh"
#include "DetectorConstructionUtils.hh"

#include <G4VUserDetectorConstruction.hh>
#include <G4UImessenger.hh>
#include <G4UIdirectory.hh>
#include <G4UIcommand.hh>
#include <G4UIcmdWithADoubleAndUnit.hh>
#include <G4UIcmdWithADouble.hh>
#include <G4UIcmdWith3VectorAndUnit.hh>
#include <G4UIcmdWithABool.hh>
#include <G4UIcmdWithAString.hh>

class DetectorConstruction : public G4VUserDetectorConstruction, G4UImessenger, MaterialUser {
public:
	/// constructor
	DetectorConstruction();
	
	/// construct detector geometry
	G4VPhysicalVolume* Construct();
	/// UI interface
	virtual void SetNewValue(G4UIcommand * command,G4String newValue);
	
	// world volume
	G4LogicalVolume* experimentalHall_log;	
	G4VPhysicalVolume* experimentalHall_phys;
	// components
	DecayTrapConstruction trap;
	SourceHolderConstruction source;
	G4VPhysicalVolume* source_phys;
	DetectorPackageConstruction dets[2];
	G4VPhysicalVolume* detPackage_phys[2];
	SiliconDetectorConstruction siDet;
	G4VPhysicalVolume* siDet_phys;
	
	/// get source holder position
	G4ThreeVector getHolderPos() const { return fSourceHolderPos; }
	
private:
	/// construct detector (Electro-)Magnetic Field
	void ConstructField();  
	
	Field* fpMagField;								///< magnetic field
	
	// sensitive volumes
	TrackerSD* scint_SD[2];
	TrackerSD* Dscint_SD[2];
	TrackerSD* backing_SD[2];
	
	TrackerSD* winIn_SD[2];
	TrackerSD* winOut_SD[2];
	TrackerSD* trap_win_SD[2];
	TrackerSD* kevlar_SD[2];
	
	TrackerSD* mwpc_SD[2];
	TrackerSD* mwpc_planes_SD[2];
	TrackerSD* mwpcDead_SD[2];
	
	TrackerSD* source_SD;
	TrackerSD* trap_monitor_SD[2];
	
	TrackerSD* hall_SD;
	
	TrackerSD* siDet_SD;
	
	// UI commands
	G4UIdirectory* fDetectorDir;					///< UI Directory for detector-related commands
	
	G4UIcmdWithAString* fDetectorGeometry;			///< which detector geometry to construct
	G4String sGeometry;
	
	G4UIcmdWith3VectorAndUnit* fSourceHolderPosCmd;	///< source holder position
	G4ThreeVector fSourceHolderPos;
	
	G4UIcmdWith3VectorAndUnit* fDetOffsetCmd;		///< Symmetrical detector offset from center axis
	G4ThreeVector fDetOffset;
	
	G4UIcmdWithADouble* fDetRotCmd;					///< Symmetrical detector rotation angle around Z axis (radians)
	Float_t fDetRot;
	
	G4UIcmdWithABool* fInFoilCmd;					///< construction of Indium 10um Al source foil
	bool makeInFoil;
	
	G4UIcmdWithADoubleAndUnit* fVacuumLevelCmd;		///< SCS bore vacuum
	Float_t fVacuumPressure;
	
	G4UIcmdWithADoubleAndUnit* fScintStepLimitCmd;	///< step size limiter in scintillator
	Float_t fScintStepLimit;
	
	G4UIcmdWithADoubleAndUnit* fMWPCBowingCmd;		///< additional width of MWPC due to window bowing
	Float_t fMWPCBowing;
	
	G4UIcmdWithADoubleAndUnit* fSourceFoilThickCmd;	///< source foil full thickness
	Float_t fSourceFoilThick;

	G4UIcmdWithADouble* fCrinkleAngleCmd;			///< decay trap foil crinkle angle
	Float_t fCrinkleAngle;
};

#endif

