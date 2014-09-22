#include "SMExcept.hh"
#include <cassert>

#include "AnalysisManager.hh"
#include "TrackerSD.hh"
#include "Field.hh"
#include "DetectorConstruction.hh"

#include <G4SubtractionSolid.hh>
#include <G4SDManager.hh>
#include <G4RunManager.hh>

#include <G4MagneticField.hh>
#include <G4FieldManager.hh>
#include <G4ChordFinder.hh>
#include <G4PropagatorInField.hh>
#include <G4TransportationManager.hh>
#include <G4UserLimits.hh>
#include <G4PVParameterised.hh>

DetectorConstruction::DetectorConstruction(): fpMagField(NULL) {
	
	fDetectorDir = new G4UIdirectory("/detector/");
	fDetectorDir->SetGuidance("/detector control");
	
	fDetectorGeometry = new G4UIcmdWithAString("/detector/geometry",this);
	fDetectorGeometry->SetGuidance("Set the geometry of the detector");
	fDetectorGeometry->AvailableForStates(G4State_PreInit);
	sGeometry = "C";
	
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
	
	fScintStepLimitCmd = new G4UIcmdWithADoubleAndUnit("/detector/scintstepsize",this);
	fScintStepLimitCmd->SetGuidance("step size limit in scintillator, windows");
	fScintStepLimitCmd->SetDefaultValue(1.0*mm);
	
	experimentalHall_log = NULL;
	experimentalHall_phys = NULL;
}

void DetectorConstruction::SetNewValue(G4UIcommand * command, G4String newValue) {
	if (command == fDetectorGeometry) {
		sGeometry = G4String(newValue);
	} else if (command == fDetOffsetCmd) {
		fDetOffset = fDetOffsetCmd->GetNew3VectorValue(newValue);
		G4cout << "Setting detector offsets to " << fDetOffset/mm << " mm" << G4endl;
	} else if (command == fDetRotCmd) {
		fDetRot = fDetRotCmd->GetNewDoubleValue(newValue);
		G4cout << "Setting detector rotation to " << fDetRot << " radians" << G4endl;
	} else if (command == fSourceHolderPosCmd) {
		fSourceHolderPos = fSourceHolderPosCmd->GetNew3VectorValue(newValue);
		G4cout<<"setting the source at "<<fSourceHolderPos/mm << " mm" << G4endl;
	} else if (command == fVacuumLevelCmd) {
		fVacuumPressure = fVacuumLevelCmd->GetNewDoubleValue(newValue);
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
		G4cerr << "Unknown command:" << command->GetCommandName() << " passed to DetectorConstruction::SetNewValue\n";
    }
}

TrackerSD* registerSD(G4String sdName) {
	TrackerSD* sd = new TrackerSD(sdName);
	G4SDManager::GetSDMpointer()->AddNewDetector(sd);
	gAnalysisManager->SaveSDName(sdName);
	return sd;
}

G4VPhysicalVolume* DetectorConstruction::Construct() {
	
	//------------------------------------------------------ materials
	setVacuumPressure(fVacuumPressure);
	
	////////////////////////////////////////
	// user step limits
	////////////////////////////////////////	
	G4UserLimits* UserCoarseLimits = new G4UserLimits();
	UserCoarseLimits->SetMaxAllowedStep(10*m);
	G4UserLimits* UserGasLimits = new G4UserLimits();
	UserGasLimits->SetMaxAllowedStep(1*cm);
	G4UserLimits* UserSolidLimits = new G4UserLimits();
	UserSolidLimits->SetMaxAllowedStep(fScintStepLimit);
		
	///////////////////////////////////////
	//experimental Hall
	///////////////////////////////////////
	const G4double expHall_halfx=1.0*m;
	const G4double expHall_halfy=1.0*m;
	const G4double expHall_halfz=4.0*m;  	
	G4Box* experimentalHall_box = new G4Box("expHall_box",expHall_halfx,expHall_halfy,expHall_halfz);
	experimentalHall_log = new G4LogicalVolume(experimentalHall_box,Vacuum,"World_Log");  
	experimentalHall_log->SetVisAttributes(G4VisAttributes::Invisible);
	experimentalHall_log->SetUserLimits(UserCoarseLimits);
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
	trap.fCoatingThick[EAST]=trap.fCoatingThick[WEST]=0.150*um;
        trap.fWindowThick[EAST]=trap.fWindowThick[WEST]=0.50*um;
        if(sGeometry=="A"){
		// thick MWPC windows
		dets[EAST].mwpc.fWindowThick=dets[WEST].mwpc.fWindowThick=25*um;
	} else if(sGeometry=="B") {
		// thick MWPC, decay trap windows
		dets[EAST].mwpc.fWindowThick=dets[WEST].mwpc.fWindowThick=25*um;
		trap.fWindowThick[EAST]=trap.fWindowThick[WEST]=(0.7+12.5)*um;
	} else if(sGeometry=="C") {
		// "default" thin-windows configuration
	} else if(sGeometry=="D") {
		// open-ended decay trap
		trap.fWindowMat=Vacuum;
		trap.fCoatingMat=Vacuum;
	} else if(sGeometry=="2007") {
		// thick windows
		dets[EAST].mwpc.fWindowThick=dets[WEST].mwpc.fWindowThick=25*um;
		trap.fWindowThick[EAST]=trap.fWindowThick[WEST]=2.5*um;
	} else if(sGeometry=="siDet") {
	
	} else if(sGeometry=="thinFoil") {
		trap.fWindowThick[EAST]=0.130*um;
                trap.fWindowThick[WEST]=0.180*um;
		trap.fWindowMat=SixFSixF;
	} else if(sGeometry=="2011/2012") {
	        trap.fWindowThick[EAST]=trap.fWindowThick[WEST]=0.500*um;
		trap.fWindowMat=Mylar;
		trap.fIRcollimator=2.3*inch;
		trap.fColThickness=0.7*inch;
		dets[EAST].mwpc.activeRegion.anode_R=dets[WEST].mwpc.activeRegion.anode_R=5*um;
		dets[EAST].mwpc.activeRegion.cathode_R=dets[WEST].mwpc.activeRegion.cathode_R=39.1*um;
	} else if(sGeometry=="2012/2013") {
		trap.fWindowThick[EAST]=0.180*um;
                trap.fWindowThick[WEST]=0.130*um;
		trap.fWindowMat=SixFSixF;
		trap.fIRcollimator=2.25*inch;
		trap.fColThickness=0.75*inch; 
		trap.fColLength=0.25*inch;
		dets[EAST].mwpc.activeRegion.anode_R=dets[WEST].mwpc.activeRegion.anode_R=5.*um;
		dets[EAST].mwpc.activeRegion.cathode_R=dets[WEST].mwpc.activeRegion.cathode_R=39.1*um;
	}
	else if(sGeometry=="2012/2013_isobutane") {
		trap.fWindowThick[EAST]=0.180*um;
                trap.fWindowThick[WEST]=0.130*um;
		trap.fWindowMat=SixFSixF;
		trap.fIRcollimator=2.25*inch;
		trap.fColThickness=0.75*inch; 
		trap.fColLength=0.25*inch;
		dets[EAST].mwpc.activeRegion.anode_R=dets[WEST].mwpc.activeRegion.anode_R=5.*um;
		dets[EAST].mwpc.activeRegion.cathode_R=dets[WEST].mwpc.activeRegion.cathode_R=39.1*um;
		dets[EAST].mwpc.fMWPCGas=dets[WEST].mwpc.fMWPCGas=WCButane;
	}
	else {
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
		for(Side sd = EAST; sd <= WEST; ++sd) {
			dets[sd].mwpc.entranceToCathodes += fMWPCBowing;
			//dets[sd].mwpc.exitToCathodes += fMWPCBowing/2.;
			
			dets[sd].Construct(sd);
			G4RotationMatrix* sideFlip = new G4RotationMatrix();
			sideFlip->rotateZ(fDetRot*ssign(sd)*rad);
			if(sd==EAST)
				sideFlip->rotateY(M_PI*rad);
			G4ThreeVector sideTrans = G4ThreeVector(0.,0.,ssign(sd)*(2.2*m-dets[sd].getScintFacePos()))+fDetOffset*ssign(sd);
			
			detPackage_phys[sd] = new G4PVPlacement(sideFlip,sideTrans,
												   dets[sd].container_log,sideSubst("detPackage_phys%c",sd),experimentalHall_log,false,0);
			
			dets[sd].mwpc.myRotation = sideFlip;
			dets[sd].mwpc.myTranslation = (*sideFlip)(dets[sd].mwpc.myTranslation);
			dets[sd].mwpc.myTranslation += sideTrans;
			dets[sd].mwpc.setPotential(2700*volt);
			
			trap.trap_win_log[sd]->SetUserLimits(UserSolidLimits);
			dets[sd].mwpc.container_log->SetUserLimits(UserGasLimits);
			dets[sd].mwpc.winIn_log->SetUserLimits(UserSolidLimits);
			dets[sd].mwpc.winOut_log->SetUserLimits(UserSolidLimits);
			dets[sd].mwpc.kevStrip_log->SetUserLimits(UserSolidLimits);
			dets[sd].scint.container_log->SetUserLimits(UserSolidLimits);
		}
		
		for(Side sd = EAST; sd <= WEST; ++sd ) {
			
			scint_SD[sd] = registerSD(sideSubst("scint_SD%c",sd));
			dets[sd].scint.scint_log->SetSensitiveDetector(scint_SD[sd]);
			
			Dscint_SD[sd] = registerSD(sideSubst("Dscint_SD%c",sd));
			dets[sd].scint.Dscint_log->SetSensitiveDetector(Dscint_SD[sd]);
			dets[sd].scint.container_log->SetSensitiveDetector(Dscint_SD[sd]);
			dets[sd].mwpc_exit_N2_log->SetSensitiveDetector(Dscint_SD[sd]);		// include N2 volume here
			dets[sd].scint.lightguide_log->SetSensitiveDetector(Dscint_SD[sd]);	// and also light guides

			backing_SD[sd] = registerSD(sideSubst("backing_SD%c",sd));
			dets[sd].scint.backing_log->SetSensitiveDetector(backing_SD[sd]);
			
			winOut_SD[sd] = registerSD(sideSubst("winOut_SD%c",sd));
			dets[sd].mwpc.winOut_log->SetSensitiveDetector(winOut_SD[sd]);
			
			winIn_SD[sd] = registerSD(sideSubst("winIn_SD%c",sd));
			dets[sd].mwpc.winIn_log->SetSensitiveDetector(winIn_SD[sd]);
			
			trap_win_SD[sd] = registerSD(sideSubst("trap_win_SD%c",sd));
			trap.mylar_win_log[sd]->SetSensitiveDetector(trap_win_SD[sd]);
			trap.be_win_log[sd]->SetSensitiveDetector(trap_win_SD[sd]);
			trap.wigglefoils[sd].SetSensitiveDetector(trap_win_SD[sd]);
			
			mwpc_SD[sd] = registerSD(sideSubst("mwpc_SD%c",sd));
			dets[sd].mwpc.activeRegion.gas_log->SetSensitiveDetector(mwpc_SD[sd]);
			dets[sd].mwpc.activeRegion.anodeSeg_log->SetSensitiveDetector(mwpc_SD[sd]);
			dets[sd].mwpc.activeRegion.cathSeg_log->SetSensitiveDetector(mwpc_SD[sd]);
			
			mwpc_planes_SD[sd] = registerSD(sideSubst("mwpc_planes_SD%c",sd));
			dets[sd].mwpc.activeRegion.cathode_wire_log->SetSensitiveDetector(mwpc_planes_SD[sd]);
			dets[sd].mwpc.activeRegion.cath_plate_log->SetSensitiveDetector(mwpc_planes_SD[sd]);
			dets[sd].mwpc.activeRegion.anode_wire_log->SetSensitiveDetector(mwpc_planes_SD[sd]);
			
			mwpcDead_SD[sd] = registerSD(sideSubst("mwpcDead_SD%c",sd));
			dets[sd].mwpc.container_log->SetSensitiveDetector(mwpcDead_SD[sd]);
			
			kevlar_SD[sd] = registerSD(sideSubst("kevlar_SD%c",sd));
			dets[sd].mwpc.kevStrip_log->SetSensitiveDetector(kevlar_SD[sd]);
			
		}
		
		// source holder
		source_SD = registerSD("source_SD");
		source.window_log->SetSensitiveDetector(source_SD);
		for(Side sd = EAST; sd <= WEST; ++sd)
			source.coating_log[sd]->SetSensitiveDetector(source_SD);
		
		// decay trap monitor volumes
		for(Side sd = EAST; sd <= WEST; ++sd ) {
			trap_monitor_SD[sd] = registerSD(sideSubst("trap_monitor_SD%c",sd));
			trap.trap_monitor_log[sd]->SetSensitiveDetector(trap_monitor_SD[sd]);
		}
		
		// experimental hall vacuum, decay tube, other inert parts
		hall_SD = registerSD("hall_SD");
		experimentalHall_log->SetSensitiveDetector(hall_SD);
		trap.decayTube_log->SetSensitiveDetector(hall_SD);
		for(Side sd = EAST; sd <= WEST; ++sd ) {
			dets[sd].mwpc_entrance_log->SetSensitiveDetector(hall_SD);
			dets[sd].mwpc_exit_log->SetSensitiveDetector(hall_SD);
			dets[sd].container_log->SetSensitiveDetector(hall_SD);
			trap.collimator_log[sd]->SetSensitiveDetector(hall_SD);
			trap.collimatorBack_log[sd]->SetSensitiveDetector(hall_SD);
		}
		
		ConstructField();
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

void DetectorConstruction::ConstructField() {
	
	if(!fpMagField) {
		cout << "##### Constructing Field #####" << endl;
		
		// get magnetic field profile
		fpMagField = new Field();
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
				
		for(Side sd = EAST; sd <= WEST; ++sd) {
			dets[sd].mwpc.myBField = fpMagField;
			dets[sd].mwpc.ConstructField();
		}
	}
}
