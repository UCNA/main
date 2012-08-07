#include "bmScintillatorConstruction.hh"
#include <G4Polycone.hh>
#include <cassert>

void bmScintillatorConstruction::Construct(Side s) {
	
	assert(backing_Radius > scint_Radius);
	assert(lightguide_thick >= scint_thick);
	N2_volume_Z = lightguide_thick+backing_thick;	// length of N2 volume
	scintFacePos = -N2_volume_Z/2;					// scintillator face position, at edge of N2 volume
		
	// overall container
	G4Tubs* N2_vol_tube = new G4Tubs(sideSubst("N2_vol_tube%c",s),0.,backing_Radius,N2_volume_Z/2,0.,2*M_PI);
	container_log = new G4LogicalVolume(N2_vol_tube,WCNitrogen,sideSubst("N2_vol_log%c",s));
	container_log->SetVisAttributes(G4VisAttributes::Invisible);
	
	// dead layer
	G4Tubs* Dscint_tube = new G4Tubs(sideSubst("Dscint_tube%c",s),0.,scint_Radius,dead_thick/2,0.,2*M_PI);
	Dscint_log = new G4LogicalVolume(Dscint_tube,Sci,sideSubst("Dscint_log%c",s));
	G4VisAttributes* visDScint= new G4VisAttributes(G4Colour(1.0,0.0,1.0,0.5));
	Dscint_log->SetVisAttributes(visDScint);
	Dscint_phys = new G4PVPlacement(NULL,G4ThreeVector(0.,0.,-(N2_volume_Z-dead_thick)/2.),
									Dscint_log,sideSubst("Dscint_phys%c",s),container_log,false,0);

	// scintillator
	G4Tubs* scint_tube = new G4Tubs("scint_tube",0.,scint_Radius,(scint_thick-dead_thick)/2.,0.,2*M_PI);
	scint_log = new G4LogicalVolume(scint_tube,Sci,sideSubst("scint_log%c",s));
	G4VisAttributes* visScint= new G4VisAttributes(G4Colour(0.0,1.0,1.0,0.2));
	scint_log->SetVisAttributes(visScint);
	scint_phys = new G4PVPlacement(NULL,G4ThreeVector(0.,0.,-N2_volume_Z/2+dead_thick+(scint_thick-dead_thick)/2.),
								   scint_log,sideSubst("scint_phys%c",s),container_log,false,0);

	// light guides around/behind scintillator
	const G4double zPlane[] = {0.,			scint_thick,	scint_thick,	lightguide_thick};
	G4double lg_rad = scint_Radius-(lightguide_thick-scint_thick);
	const G4double rInner[] = {scint_Radius,	scint_Radius,	lg_rad,			lg_rad};
	const G4double rOuter[] = {backing_Radius,backing_Radius,	backing_Radius,	backing_Radius};
	G4Polycone* lightguide_polycone = new G4Polycone("lightguide_polycone",0,2*M_PI,4,zPlane,rInner,rOuter);
	lightguide_log = new G4LogicalVolume(lightguide_polycone,Sci,sideSubst("lightguide_log%c",s));
	G4VisAttributes* visLG = new G4VisAttributes(G4Colour(0.0,1.0,0.5,0.2));
	lightguide_log->SetVisAttributes(visLG);
	lightguide_phys = new G4PVPlacement(NULL,G4ThreeVector(0.,0.,-N2_volume_Z/2.),
										lightguide_log,sideSubst("scint_phys%c",s),container_log,false,0);
	
	// backing veto
	G4Tubs* backing_tube = new G4Tubs("backing_tube",0.,backing_Radius,backing_thick/2.,0.,2*M_PI);
	backing_log = new G4LogicalVolume(backing_tube,Sci,sideSubst("backing_log%c",s));
	G4VisAttributes* visBacking= new G4VisAttributes(G4Colour(0.0,0.0,1,0.2));
	backing_log->SetVisAttributes(visBacking);
	backing_phys = new G4PVPlacement(NULL,G4ThreeVector(0.,0.,(N2_volume_Z-backing_thick)/2.),
										backing_log,sideSubst("backing_phys%c",s),container_log,false,0);
	
}
