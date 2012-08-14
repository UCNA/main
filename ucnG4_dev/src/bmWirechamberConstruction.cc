#include "bmWirechamberConstruction.hh"
#include "G4PVReplica.hh"
#include "G4FieldManager.hh"
#include "G4ChordFinder.hh"
#include "G4EqMagElectricField.hh"
#include "G4ClassicalRK4.hh"
#include <math.h>

void bmWirechamberConstruction::Construct(Side s) {
	
	///////////////////////////////////////////////////
	// main volume
	///////////////////////////////////////////////////
	
	// construct active gas volume
	activeRegion.fMWPCGas = fMWPCGas;
	activeRegion.Construct(s);
	
	mwpcContainer_halfZ = 0.5*(entranceToCathodes+exitToCathodes+activeRegion.GetWidth());
	const G4double mwpc_volume_width=8.0*inch;		// MWPC gas box width
	
	// container volume for all MWPC
	G4Box* mwpcContainer_box = new G4Box("mwpcContainer_box",mwpc_volume_width/2,mwpc_volume_width/2,mwpcContainer_halfZ);
	container_log = new G4LogicalVolume(mwpcContainer_box, fMWPCGas,sideSubst("mwpcContainer_log%c",s));
	container_log->SetVisAttributes(G4VisAttributes::Invisible);
	
	// MWPC active gas volume placement with wireplane, relative to MWPC container volume
	myTranslation = G4ThreeVector(0.,0.,(entranceToCathodes-exitToCathodes)/2);
	new G4PVPlacement(NULL,myTranslation,
					  activeRegion.gas_log,sideSubst("mwpc_phys%c",s),container_log,false,0);
	
	d = activeRegion.spacing;
	L = activeRegion.planeSpacing;
	r = activeRegion.anode_R;
	
	///////////////////////////////////////////////////
	// kevlar strings
	///////////////////////////////////////////////////
	
	// rectangular cross section strings with equal volume to nominal 140um cylinders
	const G4double kevlar_R=0.07*mm;
	const G4double kevlar_spacing=5.*mm;
	const G4int NbOfKevWires=32;
	const G4double kevLength=15.*cm;
	const G4double kev_AR = 16.;						// aspect ratio, width:depth
	const G4double kev_area = PI*kevlar_R*kevlar_R;		// total cross section area
	const G4double kev_eff_w = sqrt(kev_area*kev_AR);	// effective width
	const G4double kev_eff_t = sqrt(kev_area/kev_AR);	// effective thickness
	
	const G4double kevlarPosZ = -mwpcContainer_halfZ+kev_eff_t/2.;
	
	G4Box* kevContainer_box = new G4Box("kevContainer_box",NbOfKevWires*kevlar_spacing/2.,kevLength/2.,kev_eff_t/2.);
	G4Box* kevSeg_box = new G4Box("kevSeg_box",kevlar_spacing/2.,kevLength/2.,kev_eff_t/2);
	G4Box* kevStrip_box = new G4Box("kevStrip_box",kev_eff_w/2.,kevLength/2.,kev_eff_t/2.);
	
	kevContainer_log = new G4LogicalVolume(kevContainer_box,Vacuum,sideSubst("kevContainer_log%c",s));
	kevSeg_log = new G4LogicalVolume(kevSeg_box,Vacuum,"kevSeg_log");
	kevStrip_log = new G4LogicalVolume(kevStrip_box,Kevlar,"kevStrip_log");
	
	G4VisAttributes* visKevlar = new G4VisAttributes(G4Colour(1,1.0,0,1.0));
	kevStrip_log->SetVisAttributes(visKevlar);
	
	// place components and replicate array
	new G4PVPlacement(NULL,G4ThreeVector(0.,0.,kevlarPosZ),
					  kevContainer_log,sideSubst("kevContainer_phys%c",s),container_log,false,0);
	new G4PVPlacement(NULL,G4ThreeVector(0.,0.,0.),
					  kevStrip_log,"kevStrip_phys",kevSeg_log,false,0);	
	new G4PVReplica(sideSubst("kevlar_plane_%c",s),
					kevSeg_log,
					kevContainer_log,
					kXAxis,
					NbOfKevWires,
					kevlar_spacing);

	///////////////////////////////////////////////////
	// mylar windows
	///////////////////////////////////////////////////
	
	G4Tubs* winInner_tube = new G4Tubs("winInner_tube",0.,mwpc_entrance_R,fWindowThick/2,0.,2*M_PI);  
	G4Tubs* winOuter_tube = new G4Tubs("winOuter_tube",0.,mwpc_exit_R,fWindowThick/2,0.,2*M_PI); 
	G4VisAttributes* visWindow = new G4VisAttributes(G4Colour(0,1.0,0,1));
	winIn_log = new G4LogicalVolume(winInner_tube,Mylar,sideSubst("winIn_log%c",s));
	winIn_log->SetVisAttributes(visWindow);
	new G4PVPlacement(NULL,G4ThreeVector(0.,0.,-mwpcContainer_halfZ+kev_eff_t+fWindowThick/2),winIn_log,sideSubst("winIn%c",s),
					  container_log,false,0);
	winOut_log = new G4LogicalVolume(winOuter_tube,Mylar,sideSubst("winOut_log%c",s));
	winOut_log->SetVisAttributes(visWindow);
	new G4PVPlacement(NULL,G4ThreeVector(0.,0.,mwpcContainer_halfZ-fWindowThick/2),winOut_log,sideSubst("winOut%c",s),
					  container_log,false,0);
	
}

void bmWirechamberConstruction::GetFieldValue(const G4double Point[4], G4double* Bfield) const {
	// set magnetic field
	if(myBField) myBField->GetFieldValue(Point,Bfield);
	else Bfield[0]=Bfield[1]=Bfield[2]=0;
	if(!E0) { Bfield[3]=Bfield[4]=Bfield[5]=0; return; }
	
	// local position
	G4ThreeVector localPos = G4ThreeVector(Point[0],Point[1],Point[2])-myTranslation;
	if(myRotation) localPos = (*myRotation)(localPos);

	// electric field components
	G4ThreeVector E(0,0,0);	
	double l = localPos[2];
	if(fabs(l)<L) {
		double a = localPos[0]/d;
		a = (a-floor(a)-0.5)*d;
		if(a*a+l*l > r*r) {
			double denom = cosh(2*M_PI*l/d)-cos(2*M_PI*a/d);
			E[2] = E0*sinh(2*M_PI*l/d)/denom;
			E[0] = E0*sin(2*M_PI*a/d)/denom;
		}
	}
	
	// return to global coordinates
	if(myRotation) E = myRotation->inverse()(E);
	for(unsigned int i=0; i<3; i++) Bfield[3+i] = E[i];
}

void bmWirechamberConstruction::setPotential(G4double Vanode) {
	E0 = M_PI*Vanode/d/log(sinh(M_PI*L/d)/sinh(M_PI*r/d));
	G4cout << "Wirechamber voltage set to " << Vanode/volt <<" V => E0 = " << E0/(volt/cm) << " V/cm" << G4endl;
}

void bmWirechamberConstruction::ConstructField() {
	G4cout << "Setting up wirechamber electromagnetic field...";
	
	// local field manager
	G4FieldManager* localFieldMgr = new G4FieldManager();
	localFieldMgr->SetDetectorField(this);
	
	// equation of motion, stepper for field
	G4EqMagElectricField* pEquation = new G4EqMagElectricField(this);
	G4ClassicalRK4* pStepper = new G4ClassicalRK4(pEquation,8);
	G4MagInt_Driver* pIntgrDriver = new G4MagInt_Driver(0.01*um,pStepper,pStepper->GetNumberOfVariables());
	G4ChordFinder* pChordFinder = new G4ChordFinder(pIntgrDriver);
	localFieldMgr->SetChordFinder(pChordFinder);
	
	// accuracy settings
	localFieldMgr->GetChordFinder()->SetDeltaChord(10*um);
	localFieldMgr->SetMinimumEpsilonStep(1e-6);
	localFieldMgr->SetMaximumEpsilonStep(1e-5);
	localFieldMgr->SetDeltaOneStep(0.1*um);
	
	// apply field manager to wirechamber and all daughter volumes
	container_log->SetFieldManager(localFieldMgr,true);
	
	G4cout << " Done." << G4endl;
}
