#include "bmDetectorPackageConstruction.hh"

void bmDetectorPackageConstruction::Construct(Side s) {	
	scint.Construct(s);
	mwpc.Construct(s);
	
	const G4double detPackageRadius = 6.0*inch;
	const G4double detPackageHalfZ = 5.0*inch;
	const G4double mwpc_pos = -(mwpc.GetWidth()+scint.GetWidth())/2-scint.getScintFacePos();
	
	G4Tubs* detPackage_tube = new G4Tubs(sideSubst("detPackage_tube%c",s),0,detPackageRadius,detPackageHalfZ,0.,2*M_PI);   
	container_log = new G4LogicalVolume(detPackage_tube,Vacuum,sideSubst("container_log%c",s));
	
	scint_phys = new G4PVPlacement(NULL,G4ThreeVector(0.,0.,-scint.getScintFacePos()),
								   scint.container_log,sideSubst("N2_vol_phys%c",s),container_log,false,0);
	mwpc_phys = new G4PVPlacement(NULL,G4ThreeVector(0.,0.,mwpc_pos),
								  mwpc.container_log,sideSubst("mwpcContainer%c",s),container_log,false,0);
	
	// aluminum entrance, exit windows around container
	const G4double mwpc_entrance_thickness = 2.0*inch;
	const G4double mwpc_exit_thickness = 0.5*inch;
	G4Tubs* mwpc_entrance_tube = new G4Tubs("mwpc_entrance_tube",mwpc.mwpc_entrance_R+0.1*mm,detPackageRadius,0.5*mwpc_entrance_thickness,0.,2*M_PI); 
	G4Tubs* mwpc_exit_tube = new G4Tubs("mwpc_exit_tube",mwpc.mwpc_exit_R+0.1*mm,detPackageRadius,0.5*mwpc_exit_thickness,0.,2*M_PI);
	G4VisAttributes* visMWPCEntrance = new G4VisAttributes(G4Colour(0.7,0.7,0.7,0.8));
	G4VisAttributes* visMWPCExit = new G4VisAttributes(G4Colour(0.3,0.3,0.3,0.8));
	mwpc_entrance_log = new G4LogicalVolume(mwpc_entrance_tube, Al, sideSubst("mwpc_entrance_log%c",s));
	mwpc_entrance_log->SetVisAttributes(visMWPCEntrance);
	mwpc_entrance_phys = new G4PVPlacement(NULL,G4ThreeVector(0.,0.,mwpc_pos-(mwpc.GetWidth()+mwpc_entrance_thickness)/2),
										   mwpc_entrance_log,sideSubst("mwpc_entrance%c",s),container_log,false,0);
	mwpc_exit_log = new G4LogicalVolume(mwpc_exit_tube, Al, sideSubst("mwpc_exit_log%c",s));
	mwpc_exit_log->SetVisAttributes(visMWPCExit);
	mwpc_exit_phys = new G4PVPlacement(NULL,G4ThreeVector(0.,0.,mwpc_pos+(mwpc.GetWidth()+mwpc_exit_thickness)/2),
									   mwpc_exit_log,sideSubst("mwpc_exit%c",s),container_log,false,0);	
}
