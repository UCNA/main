#include "bmScintillatorConstruction.hh"

void bmScintillatorConstruction::Construct(Side s) {
	const G4double scint_Radius = 7.5*cm;	// scintillator disc radius
	const G4double scint_thick = 3.5*mm;	// scintillator disc thickness
	const G4double windowToScint = 1*cm;	// nitrogen volume gap between MWPC window and scintillator
	const G4double dead_thick = 3.0*um;		// dead scintillator thickness, 3 um according to Junhua's thesis
	const G4double scintToBacking = 0.5*cm;	// scintillator-to-backing-veto gap
	const G4double backing_thick = 2.54*cm;	//just a guess
	N2_volume_Z = windowToScint+dead_thick+scint_thick+scintToBacking+backing_thick;	// length of N2 volume
	
	scintFacePos = -N2_volume_Z/2+windowToScint;
	
	const G4double DscintPosZ = scintFacePos+dead_thick/2;
	const G4double scintPosZ = DscintPosZ+dead_thick/2+scint_thick/2;
	const G4double backingPosZ = (N2_volume_Z-backing_thick)/2;
	
	G4Tubs* N2_vol_tube = new G4Tubs(sideSubst("N2_vol_tube%c",s),0.,scint_Radius,N2_volume_Z/2,0.,2*M_PI);
	G4Tubs* Dscint_tube = new G4Tubs(sideSubst("Dscint_tube%c",s),0.,scint_Radius,dead_thick/2,0.,2*M_PI);
	G4VisAttributes* visDScint= new G4VisAttributes(G4Colour(1.0,0.0,1.0,0.5));
	G4Tubs* scint_tube = new G4Tubs("scint_tube",0.,scint_Radius,scint_thick/2,0.,2*M_PI);
	G4VisAttributes* visScint= new G4VisAttributes(G4Colour(0.0,1.0,1.0,0.2));
	G4Tubs* backing_tube = new G4Tubs("backing_tube",0.,scint_Radius,backing_thick/2,0.,2*M_PI);
	G4VisAttributes* visBacking= new G4VisAttributes(G4Colour(0.0,0.0,1,0.2));
	
	container_log = new G4LogicalVolume(N2_vol_tube,WCNitrogen,sideSubst("N2_vol_log%c",s));
		
	Dscint_log = new G4LogicalVolume(Dscint_tube,WCNitrogen,sideSubst("Dscint_log%c",s));
	Dscint_log->SetVisAttributes(visDScint);
	Dscint_phys = new G4PVPlacement(NULL,G4ThreeVector(0.,0.,DscintPosZ),
									   Dscint_log,sideSubst("Dscint_phys%c",s),container_log,false,0);
	
	scint_log = new G4LogicalVolume(scint_tube,Sci,sideSubst("scint_log%c",s));
	scint_log->SetVisAttributes(visScint);
	scint_phys = new G4PVPlacement(NULL,G4ThreeVector(0.,0.,scintPosZ),
									  scint_log,sideSubst("scint_phys%c",s),container_log,false,0);
	
	backing_log = new G4LogicalVolume(backing_tube,Sci,sideSubst("backing_log%c",s));
	backing_log->SetVisAttributes(visBacking);
	backing_phys = new G4PVPlacement(NULL,G4ThreeVector(0.,0.,backingPosZ),
										backing_log,sideSubst("backing_phys%c",s),container_log,false,0);
	
}
