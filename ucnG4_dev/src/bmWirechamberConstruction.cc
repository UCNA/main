#include "bmWirechamberConstruction.hh"
#include "G4PVReplica.hh"

void bmWirechamberConstruction::Construct(Side s) {
	
	///////////////////////////////////////////////////
	// main volume
	///////////////////////////////////////////////////
	
	// construct active gas volume
	activeRegion.fMWPCGas = fMWPCGas;
	activeRegion.Construct(s);
	
	//assume 8 mm between front window and front cathode, 20 mm between cathodes,
	//5 mm from back cath to back window, 2 mm N2 volume.
	//see Junhua's thesis
	const G4double entranceToCathodes = 5.0*mm;		// entrance-window-to-cathode distance
	const G4double exitToCathodes = 5.0*mm;			// exit-window-to-cathode distance
	mwpcContainer_halfZ = 0.5*(entranceToCathodes+exitToCathodes+activeRegion.GetWidth());
	mwpc_entrance_R = 7.0*cm;						// entrance window radius to MWPC, windowed by Kevlar holder
	mwpc_exit_R = 7.5*cm;							// exit window radius to MWPC
	const G4double mwpc_volume_width=8.0*inch;		// MWPC gas box width
	
	// container volume for all MWPC
	G4Box* mwpcContainer_box = new G4Box("mwpcContainer_box",mwpc_volume_width/2,mwpc_volume_width/2,mwpcContainer_halfZ);
	container_log = new G4LogicalVolume(mwpcContainer_box, fMWPCGas,sideSubst("mwpcContainer_log%c",s));
	
	// MWPC active gas volume placement with wireplane, relative to MWPC container volume
	new G4PVPlacement(NULL,G4ThreeVector(0.,0.,(entranceToCathodes-exitToCathodes)/2),
					  activeRegion.gas_log,sideSubst("mwpc_phys%c",s),container_log,false,0);
	
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
