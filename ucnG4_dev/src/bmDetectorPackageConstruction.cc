#include "bmDetectorPackageConstruction.hh"

void bmDetectorPackageConstruction::Construct(Side s) {	
	scint.Construct(s);
	mwpc.Construct(s);
	
	//////////////////////
	// aluminum entrance collimator to detector package
	//////////////////////
	
	const G4double entrance_section_length = mwpc_entrance_depth+frontwin_frame_thick;
	
	G4Tubs* mwpc_entrance_tube = new G4Tubs("mwpc_entrance_tube",0,detPackageRadius,0.5*entrance_section_length,0.,2*M_PI);
	G4Tubs* entrance_front_tube = new G4Tubs("entrance_front_tube",mwpc_entrance_r+mwpc_entrance_thickness,
											 detPackageRadius,0.5*mwpc_entrance_thickness,0.,2*M_PI);
	G4Tubs* entrance_mid_tube = new G4Tubs("entrance_mid_tube",mwpc_entrance_r,
										   mwpc_entrance_r+mwpc_entrance_thickness,0.5*mwpc_entrance_depth,0.,2*M_PI);
	G4Tubs* entrance_back_tube = new G4Tubs("entrance_back_tube",mwpc.mwpc_entrance_R,
											detPackageRadius,0.5*frontwin_frame_thick,0.,2*M_PI);
	
	G4VisAttributes* visMWPCEntrance = new G4VisAttributes(G4Colour(0.7,0.7,0.7,0.8));
	mwpc_entrance_log = new G4LogicalVolume(mwpc_entrance_tube, Vacuum, sideSubst("mwpc_entrance_log%c",s));
	mwpc_entrance_log->SetVisAttributes(G4VisAttributes::Invisible);
	entrance_front_log = new G4LogicalVolume(entrance_front_tube, Al, sideSubst("entrance_front_log%c",s));
	entrance_mid_log = new G4LogicalVolume(entrance_mid_tube, Al, sideSubst("entrance_mid_log%c",s));
	entrance_back_log = new G4LogicalVolume(entrance_back_tube, Al, sideSubst("entrance_back_log%c",s));
	entrance_front_log->SetVisAttributes(visMWPCEntrance);
	entrance_mid_log->SetVisAttributes(visMWPCEntrance);
	entrance_back_log->SetVisAttributes(visMWPCEntrance);
	
	entrance_front_phys = new G4PVPlacement(NULL,G4ThreeVector(0.,0.,-0.5*(entrance_section_length-mwpc_entrance_thickness)),
											entrance_front_log,sideSubst("entrance_front%c",s),mwpc_entrance_log,false,0);
	entrance_mid_phys = new G4PVPlacement(NULL,G4ThreeVector(0.,0.,-0.5*frontwin_frame_thick),
										  entrance_mid_log,sideSubst("entrance_mid%c",s),mwpc_entrance_log,false,0);
	entrance_back_phys = new G4PVPlacement(NULL,G4ThreeVector(0.,0.,0.5*(entrance_section_length-frontwin_frame_thick)),
										   entrance_back_log,sideSubst("entrance_back%c",s),mwpc_entrance_log,false,0);
	
	
	//////////////////////
	// overall detector package
	//////////////////////
	
	const G4double detPackageHalfZ = mwpc_entrance_depth+mwpc.GetWidth()+1.0*inch;
	G4Tubs* detPackage_tube = new G4Tubs(sideSubst("detPackage_tube%c",s),0,detPackageRadius,detPackageHalfZ,0.,2*M_PI);   
	container_log = new G4LogicalVolume(detPackage_tube,Vacuum,sideSubst("container_log%c",s));
	container_log->SetVisAttributes(G4VisAttributes::Invisible);
	
	////////////////////////
	// place components relative to scintillator face at 0
	////////////////////////
	
	scint_phys = new G4PVPlacement(NULL,G4ThreeVector(0.,0.,-scint.getScintFacePos()),
								   scint.container_log,sideSubst("N2_vol_phys%c",s),container_log,false,0);
	
	const G4double mwpc_pos = -mwpc.GetWidth()/2.-backwin_frame_thick-(scint.GetWidth()/2.+scint.getScintFacePos());
	mwpc_phys = new G4PVPlacement(NULL,G4ThreeVector(0.,0.,mwpc_pos),
								  mwpc.container_log,sideSubst("mwpcContainer%c",s),container_log,false,0);
	mwpc.myTranslation[2] += mwpc_pos;

	const G4double entrance_pos = mwpc_pos-(mwpc.GetWidth()+entrance_section_length)/2;
	entrance_face_pos = entrance_pos - 0.5*entrance_section_length;
	entrance_win_pos = entrance_pos + 0.5*entrance_section_length;
	mwpc_entrance_phys = new G4PVPlacement(NULL,G4ThreeVector(0.,0.,entrance_pos),
										   mwpc_entrance_log,sideSubst("mwpc_entrance%c",s),container_log,false,0);
	
	//////////////////////
	// aluminum exit window and N2 volume at back of gas box
	//////////////////////
	
	G4Tubs* mwpc_exit_tube = new G4Tubs("mwpc_exit_tube",mwpc.mwpc_exit_R,detPackageRadius,0.5*backwin_frame_thick,0.,2*M_PI);
	G4VisAttributes* visMWPCExit = new G4VisAttributes(G4Colour(0.3,0.3,0.3,0.8));
	mwpc_exit_log = new G4LogicalVolume(mwpc_exit_tube, Al, sideSubst("mwpc_exit_log%c",s));
	mwpc_exit_log->SetVisAttributes(visMWPCExit);
	exit_frame_pos = mwpc_pos+(mwpc.GetWidth()+backwin_frame_thick)/2;
	mwpc_exit_phys = new G4PVPlacement(NULL,G4ThreeVector(0.,0.,exit_frame_pos),
									   mwpc_exit_log,sideSubst("mwpc_exit%c",s),container_log,false,0);	
	G4Tubs* mwpc_exit_N2_tube = new G4Tubs("mwpc_exit_N2_tube",0,mwpc.mwpc_exit_R,0.5*backwin_frame_thick,0.,2*M_PI);
	mwpc_exit_N2_log = new G4LogicalVolume(mwpc_exit_N2_tube, WCNitrogen, sideSubst("mwpc_exit_N2_log%c",s));
	mwpc_exit_N2_log->SetVisAttributes(G4VisAttributes::Invisible);
	mwpc_exit_N2_phys = new G4PVPlacement(NULL,G4ThreeVector(0.,0.,exit_frame_pos),
										  mwpc_exit_N2_log,sideSubst("mwpc_exit%c",s),container_log,false,0);

	
	/////////////////////
	// material behind detector
	/////////////////////
	
	const G4double backstuff_thick = 1.*inch;
	G4Tubs* backstuff_tube = new G4Tubs(sideSubst("backstuff_tube%c",s),0,detPackageRadius,backstuff_thick*0.5,0.,2*M_PI); 
	backstuff_log = new G4LogicalVolume(backstuff_tube,SS304,sideSubst("backstuff_log%c",s));
	backstuff_phys = new G4PVPlacement(NULL,G4ThreeVector(0.,0.,detPackageHalfZ-0.5*backstuff_thick),
									   backstuff_log,sideSubst("backstuff%c",s),container_log,false,0);
}
