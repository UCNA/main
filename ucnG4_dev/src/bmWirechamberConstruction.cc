#include "bmWirechamberConstruction.hh"
#include "G4PVParameterised.hh"
#include "ucnaWireParameterisation.hh"

void bmWirechamberConstruction::Construct(Side s) {
	
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
	mwpc_phys = new G4PVPlacement(NULL,G4ThreeVector(0.,0.,(entranceToCathodes-exitToCathodes)/2),
								  activeRegion.gas_log,sideSubst("mwpc_phys%c",s),container_log,false,0);
	
	// wirechamber mylar windows
	G4Tubs* winInner_tube = new G4Tubs("winInner_tube",0.,mwpc_entrance_R,fWindowThick/2,0.,2*M_PI);  
	G4Tubs* winOuter_tube = new G4Tubs("winOuter_tube",0.,mwpc_exit_R,fWindowThick/2,0.,2*M_PI); 
	G4VisAttributes* visWindow = new G4VisAttributes(G4Colour(0,1.0,0,1));
	winIn_log = new G4LogicalVolume(winInner_tube,Mylar,sideSubst("winIn_log%c",s));
	winIn_log->SetVisAttributes(visWindow);
	winIn_phys = new G4PVPlacement(NULL,G4ThreeVector(0.,0.,-mwpcContainer_halfZ+fWindowThick/2),winIn_log,sideSubst("winIn%c",s),
								   container_log,false,0);
	winOut_log = new G4LogicalVolume(winOuter_tube,Mylar,sideSubst("winOut_log%c",s));
	winOut_log->SetVisAttributes(visWindow);
	winOut_phys = new G4PVPlacement(NULL,G4ThreeVector(0.,0.,mwpcContainer_halfZ-fWindowThick/2),winOut_log,sideSubst("winOut%c",s),
									container_log,false,0);
	
	
	///////////////////////////////////////////////////
	//kevlar wires:
	///////////////////////////////////////////////////
	
	const G4double kevlar_R=0.07*mm;
	const G4int NbOfKevWires=28;
	const G4double kevlar_spacing=5.*mm;
	const G4double kevlarPosZ = -mwpcContainer_halfZ+fWindowThick+3*kevlar_R;
	G4double firstPos=0.-13.5*kevlar_spacing;
	G4double maxLength=7.5*cm;
	G4Tubs* kevContainer_tube = new G4Tubs("kevContainer_tube",0.,mwpc_entrance_R,kevlar_R,0.,2*M_PI);
	G4Tubs* kevlar_tube = new G4Tubs("kevlar_tube",0,kevlar_R,1*cm,0.,2*M_PI);
	G4VPhysicalVolume* kevContainer_phys;
	
	G4VPVParameterisation* kevlarParam;
	G4VisAttributes* visKevlar = new G4VisAttributes(G4Colour(1,1.0,0,1.0));
	
	kevContainer_log = new G4LogicalVolume(kevContainer_tube,Vacuum,sideSubst("kevContainer_log%c",s));
	kevContainer_phys = new G4PVPlacement(NULL,G4ThreeVector(0.,0.,kevlarPosZ),
										  kevContainer_log,sideSubst("kevContainer_phys%c",s),container_log,false,0);
	
	kevlar_log = new G4LogicalVolume(kevlar_tube,Kevlar,sideSubst("kevlar_log%c",s));
	kevlar_log->SetVisAttributes(visKevlar);
	
	kevlarParam = new ucnaWireParameterisation(NbOfKevWires,	// Number of Wires 
											   0.0,				// Z position
											   firstPos,		// first position On X
											   kevlar_spacing,	// spacing of centers
											   kevlar_R,		// tube Radius
											   maxLength,		// Disk Radius
											   1);				// placement along X
	
	kevlar_phys = new G4PVParameterised(sideSubst("kevlar%c",s),	// their name
										kevlar_log,					// their logical volume
										kevContainer_phys,			// Mother  volume
										kXAxis,						// Are placed along this axis 
										NbOfKevWires,				// Number of wires
										kevlarParam);				// The parametrisation
}
