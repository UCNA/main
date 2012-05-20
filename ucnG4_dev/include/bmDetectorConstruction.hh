//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: bmDetectorConstruction.hh,v 1.11 2011-10-03 04:12:16 mmendenhall Exp $
// GEANT4 tag $Name:  $
//

#ifndef bmDetectorConstruction_H
#define bmDetectorConstruction_H 1

#include "bmDetectorConstructionUtils.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4UImessenger.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithABool.hh"
#include "Rtypes.h"
#include <TString.h>
#include "bmTrackerSD.hh"
#include "bmField.hh"
#include "bmDecayTrapConstruction.hh"
#include "bmSourceHolderConstruction.hh"
#include "bmDetectorPackageConstruction.hh"
#include "bmSiliconDetectorConstruction.hh"

class bmDetectorConstruction : public G4VUserDetectorConstruction, G4UImessenger, MaterialUser {
public:
	/// constructor
	bmDetectorConstruction();
	
	/// construct detector geometry
	G4VPhysicalVolume* Construct();
	/// UI interface
	virtual void SetNewValue(G4UIcommand * command,G4String newValue);
	
private:
	/// construct detector (Electro-)Magnetic Field
	void ConstructField(const TString filename);  
	
	// world volume
	G4LogicalVolume* experimentalHall_log;	
	G4VPhysicalVolume* experimentalHall_phys;
	// components
	bmDecayTrapConstruction trap;
	bmSourceHolderConstruction source;
	G4VPhysicalVolume* source_phys;
	bmDetectorPackageConstruction dets[2];
	G4VPhysicalVolume* detPackage_phys[2];
	bmSiliconDetectorConstruction siDet;
	G4VPhysicalVolume* siDet_phys;
	
	// sensitive volumes
	bmTrackerSD* scint_SD[2];
	bmTrackerSD* Dscint_SD[2];
	bmTrackerSD* backing_SD[2];
	
	bmTrackerSD* winIn_SD[2];
	bmTrackerSD* winOut_SD[2];
	bmTrackerSD* trap_win_SD[2];
	
	bmTrackerSD* mwpc_SD[2];
	bmTrackerSD* mwpcDead_SD[2];
	
	bmTrackerSD* source_SD;
	bmTrackerSD* trap_monitor_SD[2];
	
	bmTrackerSD* hall_SD;
	
	bmTrackerSD* siDet_SD;
	
	// UI commands
	G4UIdirectory* fDetectorDir;					//< UI Directory for detector-related commands
	
	G4UIcommand* fDetectorGeometry;					//< which detector geometry to construct
	G4String sGeometry;
	
	G4UIcommand* fFieldCmd;							//< whether to turn on/off the magnetic field
	G4String fieldSwitch;
	
	G4UIcommand* fFieldMapFileCmd;					//< which field map to use
	TString sFieldMapFile;	
	
	G4UIcmdWith3VectorAndUnit* fSourceHolderPosCmd;	//< source holder position
	G4ThreeVector fSourceHolderPos;
	
	G4UIcmdWithABool* fInFoilCmd;					//< construction of Indium 10um Al source foil
	bool makeInFoil;
	
	G4UIcmdWithADoubleAndUnit* fVacuumLevelCmd;		//< SCS bore vacuum
	Float_t fVacuumPressure;
	
	G4UIcmdWithADoubleAndUnit* fScintStepLimitCmd;	//< step size limiter in scintillator
	Float_t fScintStepLimit;
	
	G4UIcmdWithADouble* fMatterScaleCmd[2];			//< matter interaction scaling factor
	G4double fMatterScale[2];
	
	/// turn field on/off
	void SetFieldOnOff(G4String);
	
	bmField* fpMagField;
};

#endif

