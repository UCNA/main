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
// $Id: bmDetectorConstruction.cc,v 1.41 2011-12-18 10:15:42 mmendenhall Exp $
// GEANT4 tag $Name:  $
//

#include "bmDetectorConstruction.hh"
#include "bmAnalysisManager.hh"

#include "G4SubtractionSolid.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"

#include "G4MagneticField.hh"
#include "G4FieldManager.hh"
#include "G4ChordFinder.hh"
#include "G4PropagatorInField.hh"
#include "G4TransportationManager.hh"

#include "G4UserLimits.hh"
#include "G4PVParameterised.hh"

#include "bmTrackerSD.hh"
#include "bmField.hh"

#include "SMExcept.hh"

#include <cassert>

bmDetectorConstruction::bmDetectorConstruction(): fpMagField(NULL) {
	
	fDetectorDir = new G4UIdirectory("/detector/");
	fDetectorDir->SetGuidance("/detector control");
	
	fDetectorGeometry = new G4UIcmdWithAString("/detector/geometry",this);
	fDetectorGeometry->SetGuidance("Set the geometry of the detector");
	fDetectorGeometry->AvailableForStates(G4State_PreInit);
	sGeometry = "C";
	
	fFieldCmd = new G4UIcmdWithAString("/detector/field",this);
	fFieldCmd->SetGuidance("Set B field switch");
	fieldSwitch = "on";
	
	fAFPFieldCmd = new G4UIcmdWithABool("/detector/afpfield",this);
	fAFPFieldCmd->SetGuidance("Set true to add AFP fringe field to magnetic field model");
	fAFPFieldCmd->SetDefaultValue(false);
	fAddAFPField = false;
	
	fDetOffsetCmd = new G4UIcmdWith3VectorAndUnit("/detector/offset",this);
	fDetOffsetCmd->SetGuidance("antisymmetric offset of detector packages from central axis");
	fDetOffsetCmd->SetDefaultValue(G4ThreeVector());
	fDetOffsetCmd->AvailableForStates(G4State_PreInit);
	fDetOffset = G4ThreeVector();
	
	fDetRotCmd = new G4UIcmdWithADouble("/detector/rotation",this);
	fDetRotCmd->SetGuidance("Antisymmetric rotation of detector packages around z axis");
	fDetRotCmd->SetDefaultValue(0.);
	fDetRotCmd->AvailableForStates(G4State_PreInit);
	fDetRot = 0.;

	fFieldMapFileCmd = new G4UIcmdWithAString("/detector/fieldmapfile",this);
	fFieldMapFileCmd->SetGuidance("Set B field map file");
	sFieldMapFile = "";
	
	fVacuumLevelCmd = new G4UIcmdWithADoubleAndUnit("/detector/vacuum",this);
	fVacuumLevelCmd->SetGuidance("Set SCS vacuum pressure");
	fVacuumPressure = 0;
	
	fMWPCBowingCmd = new G4UIcmdWithADoubleAndUnit("/detector/MWPCBowing",this);
	fMWPCBowingCmd->SetGuidance("Set extra wirechamber width from bowing");
	fMWPCBowingCmd->SetDefaultValue(0.);
	fMWPCBowingCmd->AvailableForStates(G4State_PreInit);
	
	fSourceHolderPosCmd = new G4UIcmdWith3VectorAndUnit("/detector/sourceholderpos",this);
	fSourceHolderPosCmd->SetGuidance("position of the source holder");
	fSourceHolderPosCmd->SetDefaultValue(G4ThreeVector());
	fSourceHolderPosCmd->AvailableForStates(G4State_PreInit);    
	
	fInFoilCmd = new G4UIcmdWithABool("/detector/infoil",this);
	fInFoilCmd->SetGuidance("Set true to build In source foil instead of usual sealed sources");
	fInFoilCmd->SetDefaultValue(false);
	
	fSourceFoilThickCmd = new G4UIcmdWithADoubleAndUnit("/detector/sourcefoilthick",this);
	fSourceFoilThickCmd->SetGuidance("Set source foil full thickness");
	fSourceFoilThick = 7.2*um;
	fSourceFoilThickCmd->SetDefaultValue(fSourceFoilThick);
	fSourceFoilThickCmd->AvailableForStates(G4State_PreInit);
	
	fCrinkleAngleCmd = new G4UIcmdWithADouble("/detector/foilcrinkle",this);
	fCrinkleAngleCmd->SetGuidance("Decay trap foil crinkle angle");
	fCrinkleAngleCmd->SetDefaultValue(0.);
	fCrinkleAngleCmd->AvailableForStates(G4State_PreInit);
	fCrinkleAngle = 0.;

	for(Side s = EAST; s <= WEST; ++s) {
		fMatterScaleCmd[s] = new G4UIcmdWithADouble(sideSubst("/detector/matterscale%c",s).c_str(),this);
		fMatterScaleCmd[s]->SetGuidance("Matter interaction scaling factor");
		fMatterScaleCmd[s]->SetDefaultValue(1.0);
		fMatterScale[s] = 1.0;
	}
	
	fScintStepLimitCmd = new G4UIcmdWithADoubleAndUnit("/detector/scintstepsize",this);
	fScintStepLimitCmd->SetGuidance("step size limit in scintillator, windows");
	fScintStepLimitCmd->SetDefaultValue(1.0*mm);
	
	experimentalHall_log = NULL;
	experimentalHall_phys = NULL;
}

void bmDetectorConstruction::SetNewValue(G4UIcommand * command, G4String newValue) {
	if (command == fDetectorGeometry) {
		sGeometry = G4String(newValue);
	} else if (command == fFieldCmd) {
		fieldSwitch = G4String(newValue);
	} else if (command == fFieldMapFileCmd) {
		sFieldMapFile = TString(newValue);
	} else if(command == fAFPFieldCmd) {
		fAddAFPField = fAFPFieldCmd->GetNewBoolValue(newValue);
		if(fpMagField) fpMagField->addAFP = fAddAFPField;
		G4cout << "Setting AFP field inclusion to " << fAddAFPField << G4endl;
	} else if (command == fDetOffsetCmd) {
		fDetOffset = fDetOffsetCmd->GetNew3VectorValue(newValue);
		G4cout << "Setting detector offsets to " << fDetOffset/mm << " mm" << G4endl;
	} else if (command == fDetRotCmd) {
		fDetRot = fDetRotCmd->GetNewDoubleValue(newValue);
		G4cout << "Setting detector rotation to " << fDetRot << " radians" << G4endl;
	} else if (command == fSourceHolderPosCmd) {
		fSourceHolderPos = fSourceHolderPosCmd->GetNew3VectorValue(newValue);
		G4cout<<"setting the source at "<<fSourceHolderPos/m<<G4endl;
	} else if (command == fVacuumLevelCmd) {
		fVacuumPressure = fVacuumLevelCmd->GetNewDoubleValue(newValue);
	} else if (command == fMatterScaleCmd[EAST]) {
		fMatterScale[EAST] = fMatterScaleCmd[EAST]->GetNewDoubleValue(newValue);
	} else if (command == fMatterScaleCmd[WEST]) {
		fMatterScale[WEST] = fMatterScaleCmd[WEST]->GetNewDoubleValue(newValue);
	} else if (command == fInFoilCmd) {
		makeInFoil = fInFoilCmd->GetNewBoolValue(newValue);
		G4cout << "Setting In source foil construction to " << makeInFoil << G4endl;
	} else if(command == fSourceFoilThickCmd) {
		fSourceFoilThick = fSourceFoilThickCmd->GetNewDoubleValue(newValue);
	} else if(command == fMWPCBowingCmd) {
		fMWPCBowing = fMWPCBowingCmd->GetNewDoubleValue(newValue);
		G4cout << "Adding " << fMWPCBowing/mm << "mm bowing to MWPC volume" << G4endl;
	} else if (command == fScintStepLimitCmd) {
		fScintStepLimit = fScintStepLimitCmd->GetNewDoubleValue(newValue);
		G4cout << "Setting step limit in solids to " << fScintStepLimit/mm << "mm" << G4endl;
	} else if (command == fCrinkleAngleCmd) {
		fCrinkleAngle = fCrinkleAngleCmd->GetNewDoubleValue(newValue);
		G4cout << "Setting decay trap foil crinkle angle to " << fCrinkleAngle << G4endl;
	}else {
		G4cerr << "Unknown command:" << command->GetCommandName() << " passed to bmDetectorConstruction::SetNewValue\n";
    }
}

bmTrackerSD* registerSD(G4String sdName) {
	bmTrackerSD* sd = new bmTrackerSD(sdName);
	G4SDManager::GetSDMpointer()->AddNewDetector(sd);
	gbmAnalysisManager->SaveSDName(sdName);
	return sd;
}

G4VPhysicalVolume* bmDetectorConstruction::Construct() {
	
	//------------------------------------------------------ materials
	setVacuumPressure(fVacuumPressure);
	
	////////////////////////////////////////
	// user step limits
	////////////////////////////////////////	
	G4UserLimits* bmUserCoarseLimits = new G4UserLimits();
	bmUserCoarseLimits->SetMaxAllowedStep(10*m);
	G4UserLimits* bmUserGasLimits = new G4UserLimits();
	bmUserGasLimits->SetMaxAllowedStep(1*cm);
	G4UserLimits* bmUserSolidLimits = new G4UserLimits();
	bmUserSolidLimits->SetMaxAllowedStep(fScintStepLimit);
		
	///////////////////////////////////////
	//experimental Hall
	///////////////////////////////////////
	const G4double expHall_halfx=1.0*m;
	const G4double expHall_halfy=1.0*m;
	const G4double expHall_halfz=4.0*m;  	
	G4Box* experimentalHall_box = new G4Box("expHall_box",expHall_halfx,expHall_halfy,expHall_halfz);
	experimentalHall_log = new G4LogicalVolume(experimentalHall_box,Vacuum,"World_Log");  
	experimentalHall_log->SetVisAttributes(G4VisAttributes::Invisible);
	experimentalHall_log->SetUserLimits(bmUserCoarseLimits);
	experimentalHall_phys = new G4PVPlacement(NULL,G4ThreeVector(),"World_Phys", experimentalHall_log,0,false,0);  
	
	
	///////////////////////////////////////
	// source holder
	///////////////////////////////////////
	source.fWindowThick = fSourceFoilThick/2.;
	if(makeInFoil) {
		G4cout << "Constructing In source foil" << G4endl;
		source.fWindowMat = source.Al;
		source.fCoatingMat = Vacuum;
		source.fWindowThick=5.0*um;
		source.fCoatingThick=0.1*um;
	}
	source.Construct();
	source_phys = new G4PVPlacement(NULL,fSourceHolderPos,source.container_log,"source_container_phys",experimentalHall_log,false,0);
	
	///////////////////////////////////////
	// geometry-dependent settings
	///////////////////////////////////////
	G4cout<<"Using geometry '"<<sGeometry<<"' ..."<<G4endl;
	if(sGeometry=="A"){
		// thick MWPC windows
		dets[EAST].mwpc.fWindowThick=dets[WEST].mwpc.fWindowThick=25*um;
	} else if(sGeometry=="B") {
		// thick MWPC, decay trap windows
		dets[EAST].mwpc.fWindowThick=dets[WEST].mwpc.fWindowThick=25*um;
		trap.fWindowThick=(0.7+12.5)*um;
	} else if(sGeometry=="C") {
		// "default" thin-windows configuration
	} else if(sGeometry=="D") {
		// open-ended decay trap
		trap.fWindowMat=Vacuum;
		trap.fCoatingMat=Vacuum;
	} else if(sGeometry=="2007") {
		// thick windows
		dets[EAST].mwpc.fWindowThick=dets[WEST].mwpc.fWindowThick=25*um;
		trap.fWindowThick=2.5*um;
	} else if(sGeometry=="siDet") {
		
	} else {
		SMExcept e("UnknownGeometry");
		e.insert("name",sGeometry);
		throw(e);
	}
	
	if(sGeometry=="siDet") {

		////////////////////////////////////////
		// silicon detector setup components
		////////////////////////////////////////
		siDet.pBackingMat = Vacuum;
		siDet.Construct();
		siDet_phys = new G4PVPlacement(NULL,G4ThreeVector(0,0,siDet.fHolderThick*0.5+source.getHolderThick()*0.5),
									   siDet.container_log,"silicon_detector_phys",experimentalHall_log,false,0);	
		
		siDet_SD = registerSD("siDet_SD");
		siDet.det_log->SetSensitiveDetector(siDet_SD);
		
	} else {
		
		////////////////////////////////////////
		// beta decay setup components
		////////////////////////////////////////
		trap.Construct(experimentalHall_log, fCrinkleAngle);
		for(Side s = EAST; s <= WEST; ++s) {
			dets[s].mwpc.entranceToCathodes += fMWPCBowing;
			//dets[s].mwpc.exitToCathodes += fMWPCBowing/2.;
			
			dets[s].Construct(s);
			G4RotationMatrix* sideFlip = new G4RotationMatrix();
			sideFlip->rotateZ(fDetRot*ssign(s)*rad);
			if(s==EAST)
				sideFlip->rotateY(M_PI*rad);
			G4ThreeVector sideTrans = G4ThreeVector(0.,0.,ssign(s)*(2.2*m-dets[s].getScintFacePos()))+fDetOffset*ssign(s);
			
			detPackage_phys[s] = new G4PVPlacement(sideFlip,sideTrans,
												   dets[s].container_log,sideSubst("detPackage_phys%c",s),experimentalHall_log,false,0);
			
			dets[s].mwpc.myRotation = sideFlip;
			dets[s].mwpc.myTranslation = (*sideFlip)(dets[s].mwpc.myTranslation);
			dets[s].mwpc.myTranslation += sideTrans;
			dets[s].mwpc.setPotential(2700*volt);
			
			trap.trap_win_log[s]->SetUserLimits(bmUserSolidLimits);
			dets[s].mwpc.container_log->SetUserLimits(bmUserGasLimits);
			dets[s].mwpc.winIn_log->SetUserLimits(bmUserSolidLimits);
			dets[s].mwpc.winOut_log->SetUserLimits(bmUserSolidLimits);
			dets[s].mwpc.kevStrip_log->SetUserLimits(bmUserSolidLimits);
			dets[s].scint.container_log->SetUserLimits(bmUserSolidLimits);
		}
		
		for(Side s = EAST; s <= WEST; ++s ) {
			
			scint_SD[s] = registerSD(sideSubst("scint_SD%c",s));
			dets[s].scint.scint_log->SetSensitiveDetector(scint_SD[s]);
			
			Dscint_SD[s] = registerSD(sideSubst("Dscint_SD%c",s));
			dets[s].scint.Dscint_log->SetSensitiveDetector(Dscint_SD[s]);
			dets[s].scint.container_log->SetSensitiveDetector(Dscint_SD[s]);
			dets[s].mwpc_exit_N2_log->SetSensitiveDetector(Dscint_SD[s]);		// include N2 volume here
			dets[s].scint.lightguide_log->SetSensitiveDetector(Dscint_SD[s]);	// and also light guides

			backing_SD[s] = registerSD(sideSubst("backing_SD%c",s));
			dets[s].scint.backing_log->SetSensitiveDetector(backing_SD[s]);
			
			winOut_SD[s] = registerSD(sideSubst("winOut_SD%c",s));
			dets[s].mwpc.winOut_log->SetSensitiveDetector(winOut_SD[s]);
			
			winIn_SD[s] = registerSD(sideSubst("winIn_SD%c",s));
			dets[s].mwpc.winIn_log->SetSensitiveDetector(winIn_SD[s]);
			
			trap_win_SD[s] = registerSD(sideSubst("trap_win_SD%c",s));
			trap.mylar_win_log[s]->SetSensitiveDetector(trap_win_SD[s]);
			trap.be_win_log[s]->SetSensitiveDetector(trap_win_SD[s]);
			trap.wigglefoils[s].SetSensitiveDetector(trap_win_SD[s]);
			
			mwpc_SD[s] = registerSD(sideSubst("mwpc_SD%c",s));
			dets[s].mwpc.activeRegion.gas_log->SetSensitiveDetector(mwpc_SD[s]);
			dets[s].mwpc.activeRegion.anodeSeg_log->SetSensitiveDetector(mwpc_SD[s]);
			dets[s].mwpc.activeRegion.cathSeg_log->SetSensitiveDetector(mwpc_SD[s]);
			
			mwpc_planes_SD[s] = registerSD(sideSubst("mwpc_planes_SD%c",s));
			dets[s].mwpc.activeRegion.cathode_wire_log->SetSensitiveDetector(mwpc_planes_SD[s]);
			dets[s].mwpc.activeRegion.cath_plate_log->SetSensitiveDetector(mwpc_planes_SD[s]);
			dets[s].mwpc.activeRegion.anode_wire_log->SetSensitiveDetector(mwpc_planes_SD[s]);
			
			mwpcDead_SD[s] = registerSD(sideSubst("mwpcDead_SD%c",s));
			dets[s].mwpc.container_log->SetSensitiveDetector(mwpcDead_SD[s]);
			
			kevlar_SD[s] = registerSD(sideSubst("kevlar_SD%c",s));
			dets[s].mwpc.kevStrip_log->SetSensitiveDetector(kevlar_SD[s]);
			
		}
		
		// source holder
		source_SD = registerSD("source_SD");
		source.window_log->SetSensitiveDetector(source_SD);
		for(Side s = EAST; s <= WEST; ++s)
			source.coating_log[s]->SetSensitiveDetector(source_SD);
		
		// decay trap monitor volumes
		for(Side s = EAST; s <= WEST; ++s ) {
			trap_monitor_SD[s] = registerSD(sideSubst("trap_monitor_SD%c",s));
			trap.trap_monitor_log[s]->SetSensitiveDetector(trap_monitor_SD[s]);
		}
		
		// experimental hall vacuum, decay tube, other inert parts
		hall_SD = registerSD("hall_SD");
		experimentalHall_log->SetSensitiveDetector(hall_SD);
		trap.decayTube_log->SetSensitiveDetector(hall_SD);
		for(Side s = EAST; s <= WEST; ++s ) {
			dets[s].mwpc_entrance_log->SetSensitiveDetector(hall_SD);
			dets[s].mwpc_exit_log->SetSensitiveDetector(hall_SD);
			dets[s].container_log->SetSensitiveDetector(hall_SD);
			trap.collimator_log[s]->SetSensitiveDetector(hall_SD);
			trap.collimatorBack_log[s]->SetSensitiveDetector(hall_SD);
		}
		
		// construct magnetic field
		cout<<"##### "<<sFieldMapFile<<" #####"<<endl;
		ConstructField(sFieldMapFile);
		//then switch the field on or off based on the UI
		SetFieldOnOff(fieldSwitch);
		fpMagField->addAFP = fAddAFPField;
	}
	
	return experimentalHall_phys;
}

#include "G4MagIntegratorStepper.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4ClassicalRK4.hh"
#include "G4SimpleHeum.hh"
#include "G4HelixHeum.hh"
#include "G4HelixImplicitEuler.hh"
#include "G4HelixExplicitEuler.hh"
#include "G4HelixSimpleRunge.hh"
#include "G4HelixMixedStepper.hh"

void bmDetectorConstruction::ConstructField(const TString filename) {
	
	static G4bool fieldIsInitialized = false;
	
	if(!fieldIsInitialized) {
		cout<<"##### Constructing Field #####"<<endl;
		
		// get magnetic field profile
		fpMagField = new bmField(filename);
		// set up field manager for this profile
		G4FieldManager* fieldMgr = G4TransportationManager::GetTransportationManager()->GetFieldManager();
		fieldMgr->SetDetectorField(fpMagField);
		// set up default chord finder
		fieldMgr->CreateChordFinder(fpMagField);
		
		// Select stepper
		G4MagIntegratorStepper* pStepper;
		G4Mag_UsualEqRhs* fEquation = new G4Mag_UsualEqRhs(fpMagField); // equation of motion in magnetic field
		//pStepper = new G4ClassicalRK4 (fEquation);		// general case for "smooth" EM fields
		//pStepper = new G4SimpleHeum( fEquation );			// for slightly less smooth EM fields
		//pStepper = new G4HelixHeum( fEquation );			// for "smooth" pure-B fields
		//pStepper = new G4HelixImplicitEuler( fEquation );	// for less smooth pure-B fields; appears ~50% faster than above
		//pStepper = new G4HelixSimpleRunge( fEquation );	// similar speed to above
		//pStepper = new G4HelixExplicitEuler( fEquation );	// about twice as fast as above
		pStepper = new G4HelixMixedStepper(fEquation,6);	// avoids "Stepsize underflow in Stepper" errors
		fieldMgr->GetChordFinder()->GetIntegrationDriver()->RenewStepperAndAdjust(pStepper);
		
		// set required accuracy for finding intersections
		fieldMgr->GetChordFinder()->SetDeltaChord(100.0*um);
		// set integration relative error limits for small and large steps
		fieldMgr->SetMinimumEpsilonStep(1e-6);
		fieldMgr->SetMaximumEpsilonStep(1e-5);
		// set integration absolute error limit
		fieldMgr->SetDeltaOneStep(0.1*um);
		// allow lots of looping
		G4TransportationManager::GetTransportationManager()->GetPropagatorInField()->SetMaxLoopCount(INT_MAX);
		
		fieldIsInitialized = true;
		
		for(Side s = EAST; s <= WEST; ++s) {
			dets[s].mwpc.myBField = fpMagField;
			dets[s].mwpc.ConstructField();
		}
	}
}

void bmDetectorConstruction::SetFieldOnOff(G4String aSwitch) {
	if(aSwitch=="on") {
		G4cout<<"##### Setting the magnetic field scale to full strength ..."<<G4endl;
		fpMagField->SetFieldScale(1.0);
	} else if (aSwitch=="off") {
		G4cout<<"##### Switching off the magnetic field ..."<<G4endl;
		fpMagField->SetFieldToZero();
	} else {
		G4cout<<"##### Setting the magnetic field scale to "<<aSwitch<<" ..."<<G4endl;
		Double_t scale = atof(aSwitch.data());
		fpMagField->SetFieldScale(scale);
	}
	cout<<"Done setting field scale"<<endl;
}
