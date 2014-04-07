#include "DecayTrapConstruction.hh"

#include <G4PVPlacement.hh>
#include <G4Tubs.hh>
#include <G4RotationMatrix.hh>

void DecayTrapConstruction::Construct(G4LogicalVolume* world, double crinkleAngle) {
	
	////////////////////////////////////////
	// decay tube
	////////////////////////////////////////
	const G4double decayTube_OR = fIRtrap+decayTube_Wall;
	const G4double decayTube_Length = 3.0*m;
	
	G4Tubs* decayTube_tube = new G4Tubs("decayTube_tube",fIRtrap, decayTube_OR, decayTube_Length/2.,0.,2*M_PI);        
	decayTube_log = new G4LogicalVolume(decayTube_tube,Cu,"decayTube_log");
	decayTube_log->SetVisAttributes(new G4VisAttributes(G4Colour(1,1,0,0.5)));
	new G4PVPlacement(NULL,G4ThreeVector(),decayTube_log,"decayTube",world,false,0);
	
	////////////////////////////////////////
	// "plug"
	////////////////////////////////////////
	if(false) {
		G4Tubs* plug_tube = new G4Tubs("plug_tube",fIRtrap-1*mm, fIRtrap, 2.*cm, -2*cm/fIRtrap, 2*cm/fIRtrap);
		plug_log = new G4LogicalVolume(plug_tube,Cu,"plug_log");
		plug_log->SetVisAttributes(new G4VisAttributes(G4Colour(1,1,0,0.5)));
		new G4PVPlacement(NULL,G4ThreeVector(),plug_log,"plug",world,false,0);
	}
	
	////////////////////////////////////////
	// Trap windows, collimator, monitors
	////////////////////////////////////////
	const G4double thicknessOfTrapWindow = fWindowThick+fCoatingThick;
	const G4double collimator_thick = 0.8*inch;
	
	G4Tubs* trap_win_tube = new G4Tubs("trap_win_tube",0.,decayTube_OR,thicknessOfTrapWindow/2.,0.,2*M_PI);  
	G4Tubs* mylarTube = new G4Tubs("mylarTube",0.,decayTube_OR,fWindowThick/2.,0.,2*M_PI);  
	G4Tubs* beTube = new G4Tubs("beTube",0.,decayTube_OR,fCoatingThick/2.,0.,2*M_PI);
	G4Tubs* collimatorTube = new G4Tubs("collimatorTube",fIRcollimator,fIRcollimator+collimator_thick,collimator_thick/2.,0.,2*M_PI);
	G4Tubs* collimatorBackTube = new G4Tubs("collimatorBackTube",decayTube_OR+1.*mm,fIRcollimator+collimator_thick,collimator_thick,0.,2*M_PI);
	
	const G4double beWinPosZ = -thicknessOfTrapWindow/2.+fCoatingThick/2.;
	const G4double mylarWinPosZ = thicknessOfTrapWindow/2.-fWindowThick/2.;
	const G4double trap_winPosZ=(decayTube_Length+thicknessOfTrapWindow)/2.;
	const G4double trap_monitor_thickness = 1.0*mm;
	const G4double trap_monitor_posZ = 0.5*m;
	G4Tubs* trap_monitor_tube = new G4Tubs("trap_monitor_tube",0.,fIRtrap,trap_monitor_thickness/2,0.,2*M_PI);	
	G4VisAttributes* visWindow = new G4VisAttributes(G4Colour(0,1.0,0,1));
	
	G4RotationMatrix* wigRot = new G4RotationMatrix;  
	wigRot->rotateY(M_PI/2.*rad);
	wigRot->rotateX(M_PI/2.*rad);
	
	for(Side sd = EAST; sd <= WEST; ++sd) {
		
		if(crinkleAngle) wigglefoils[sd].thetaMax = crinkleAngle;
		wigglefoils[sd].containerMat = Vacuum;
		wigglefoils[sd].length = 2*decayTube_OR;
		wigglefoils[sd].nseg = int(4*decayTube_OR/wigglefoils[sd].period)+1;
		wigglefoils[sd].addLayer(fWindowMat,fWindowThick,true,visWindow);
		wigglefoils[sd].addLayer(fCoatingMat,fCoatingThick,true,visWindow);
		wigglefoils[sd].Construct();
		
		trap_win_log[sd] = new G4LogicalVolume(trap_win_tube,Vacuum,sideSubst("trap_win_log%c",sd));
		trap_win_log[sd]->SetVisAttributes(visWindow);
		
		if(crinkleAngle) {
			new G4PVPlacement(wigRot,G4ThreeVector(0.,0.,ssign(sd)*(decayTube_Length+wigglefoils[sd].getContainerThick())/2),
							  wigglefoils[sd].container_log,sideSubst("trap_wigglefoil%c",sd),world,false,0);
		} else {
			new G4PVPlacement(NULL,G4ThreeVector(0.,0.,ssign(sd)*trap_winPosZ),trap_win_log[sd],
							  sideSubst("trap_win%c",sd),world,false,0);
		}
		
		mylar_win_log[sd] = new G4LogicalVolume(mylarTube, fWindowMat, sideSubst("mylar_win_log%c",sd));
		new G4PVPlacement(NULL,G4ThreeVector(0.,0.,ssign(sd)*mylarWinPosZ),mylar_win_log[sd],
						  sideSubst("mylar_win%c",sd),trap_win_log[sd],false,0);
		
		be_win_log[sd] = new G4LogicalVolume(beTube, fCoatingMat,sideSubst("be_win_log%c",sd));
		new G4PVPlacement(NULL,G4ThreeVector(0.,0.,ssign(sd)*beWinPosZ),be_win_log[sd],
						  sideSubst("be_win%c",sd),trap_win_log[sd],false,0);
		
		G4double collimatorPosZ = (decayTube_Length+collimator_thick)/2.;
		collimatorPosZ += (crinkleAngle?wigglefoils[sd].getContainerThick():thicknessOfTrapWindow)/2.;
		collimator_log[sd] = new G4LogicalVolume(collimatorTube, fCollimatorMat, sideSubst("collimator_log%c",sd));
		G4double collimatorBackZ = decayTube_Length/2.-collimator_thick;
		collimatorBack_log[sd] = new G4LogicalVolume(collimatorBackTube, fCollimatorMat, sideSubst("collimatorBack_log%c",sd));
		new G4PVPlacement(NULL,G4ThreeVector(0.,0.,ssign(sd)*collimatorPosZ),collimator_log[sd],
						  sideSubst("collimator%c",sd),world,false,0);
		new G4PVPlacement(NULL,G4ThreeVector(0.,0.,ssign(sd)*collimatorBackZ),collimatorBack_log[sd],
						  sideSubst("collimatorBack%c",sd),world,false,0);
		
		trap_monitor_log[sd] = new G4LogicalVolume(trap_monitor_tube,Vacuum,sideSubst("trap_monitor_log%c",sd));
		new G4PVPlacement(NULL,G4ThreeVector(0.,0.,ssign(sd)*trap_monitor_posZ),
						  trap_monitor_log[sd],sideSubst("trap_monitor%c",sd),world,false,0);
	}
}
