#include "bmDecayTrapConstruction.hh"

void bmDecayTrapConstruction::Construct(G4LogicalVolume* world) {
	
	////////////////////////////////////////
	// decay tube
	////////////////////////////////////////
	const G4double decayTube_IR = 7.5*cm;
	const G4double decayTube_Wall = 0.2*cm; 
	const G4double decayTube_Length = 3.0*m;
	
	G4Tubs* decayTube_tube = new G4Tubs("decayTube_tube",decayTube_IR, decayTube_IR+decayTube_Wall,decayTube_Length/2.,0.,2*M_PI);        
	decayTube_log = new G4LogicalVolume(decayTube_tube,Cu,"decayTube_log");
	decayTube_log->SetVisAttributes(new G4VisAttributes(G4Colour(1,1,0,0.5)));
	decayTube_phys = new G4PVPlacement(NULL,G4ThreeVector(),decayTube_log,"decayTube",world,false,0);
	
	////////////////////////////////////////
	// Trap windows, monitors
	////////////////////////////////////////
	const G4double thicknessOfTrapWindow = fWindowThick+fCoatingThick;
	const G4double radiusOfTrapWindow = decayTube_IR+decayTube_Wall;
	G4Tubs* trap_winInner_tube = new G4Tubs("trap_winInner_tube",0.,radiusOfTrapWindow,thicknessOfTrapWindow/2.,0.,2*M_PI);  
	const G4double trap_winPosZ=(decayTube_Length+thicknessOfTrapWindow)/2.;
	G4Tubs* mylarTube = new G4Tubs("mylarTube",0.,radiusOfTrapWindow,fWindowThick/2.,0.,2*M_PI);  
	G4Tubs* beTube = new G4Tubs("beTube",0.,radiusOfTrapWindow,fCoatingThick/2.,0.,2*M_PI);
	const G4double beWinPosZ = -thicknessOfTrapWindow/2.+fCoatingThick/2.;
	const G4double mylarWinPosZ = thicknessOfTrapWindow/2.-fWindowThick/2.;
	const G4double trap_monitor_thickness = 1.0*mm;
	const G4double trap_monitor_posZ = 0.5*m;
	G4Tubs* trap_monitor_tube = new G4Tubs("trap_monitor_tube",0.,decayTube_IR,trap_monitor_thickness/2,0.,2*M_PI);	
	G4VisAttributes* visWindow = new G4VisAttributes(G4Colour(0,1.0,0,1));
	
	for(Side s = EAST; s <= WEST; ++s) {
		
		trap_win_log[s] = new G4LogicalVolume(trap_winInner_tube,Vacuum,sideSubst("trap_win_log%c",s));
		trap_win_log[s]->SetVisAttributes(visWindow);
		trap_win_phys[s] = new G4PVPlacement(NULL,G4ThreeVector(0.,0.,ssign(s)*trap_winPosZ),trap_win_log[s],
											 sideSubst("trap_win%c",s),world,false,0);
		
		mylar_win_log[s] = new G4LogicalVolume(mylarTube, fWindowMat,sideSubst("mylar_win_log%c",s));
		mylar_win_phys[s] = new G4PVPlacement(NULL,G4ThreeVector(0.,0.,ssign(s)*mylarWinPosZ),mylar_win_log[s],
											  sideSubst("mylar_win%c",s),trap_win_log[s],false,0);
		
		be_win_log[s] = new G4LogicalVolume(beTube, fCoatingMat,sideSubst("be_win_log%c",s));
		be_win_phys[s] = new G4PVPlacement(NULL,G4ThreeVector(0.,0.,ssign(s)*beWinPosZ),be_win_log[s],
										   sideSubst("be_win%c",s),trap_win_log[s],false,0);
		
		trap_monitor_log[s] = new G4LogicalVolume(trap_monitor_tube,Vacuum,sideSubst("trap_monitor_log%c",s));
		trap_monitor_phys[s] = new G4PVPlacement(NULL,G4ThreeVector(0.,0.,ssign(s)*trap_monitor_posZ),
												 trap_monitor_log[s],sideSubst("trap_monitor%c",s),world,false,0);
	}
}
