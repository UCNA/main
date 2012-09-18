#include "bmWireVolumeConstruction.hh"
#include <cassert>

void bmWireVolumeConstruction::Construct(Side s) {
	
	G4Material* fCathodeWireMaterial = Al;			//< cathode wire material
	G4Material* fAnodeWireMaterial = Wu;			//< anode wire material
	G4Material* fCathodePlateMaterial = Au;			//< cathode plating material
	
	assert(fMWPCGas);
	assert(fCathodeWireMaterial);
	assert(fAnodeWireMaterial);	
	
	// wireplane box size
	const G4double wireplane_half_width = NbOfWires*spacing/2;
	
	//effective mwpc gas volume containing the cathodes and anode
	G4Box* mwpc_box = new G4Box("mwpc_box",wireplane_half_width,wireplane_half_width,planeSpacing);
	gas_log = new G4LogicalVolume(mwpc_box,fMWPCGas,sideSubst("mwpc_log%c",s));
	gas_log->SetVisAttributes(G4VisAttributes::Invisible);
	
	// anode, cathode wire containers... note these are "on their side" to allow wireplane parametrization,
	// and will be rotated later
	G4Box* cathContainer_box = new G4Box("cathContainer_box",wireplane_half_width,cathode_R,wireplane_half_width);
	G4Box* anodeContainer_box = new G4Box("anodeContainer_box",wireplane_half_width,anode_R,wireplane_half_width);
	// anode, cathode wires and surrounding gas
	G4Tubs* cathplate_tube = new G4Tubs("cathplate_tube",cathode_R-plating_thick,cathode_R,wireplane_half_width,0,2*M_PI);
	G4Tubs* cathode_tube = new G4Tubs("cathode_tube",0,cathode_R-plating_thick,wireplane_half_width,0.,2*M_PI);
	G4Tubs* anode_tube = new G4Tubs("anode_tube",0,anode_R,wireplane_half_width,0.,2*M_PI);
	G4Box* cathodeSegmentBox = new G4Box("cathodeSegmentBox",spacing/2,cathode_R,wireplane_half_width);
	G4Box* anodeSegmentBox = new G4Box("anodeSegmentBox",spacing/2,anode_R,wireplane_half_width);
	
	// Rotate 90 degrees around X axis
	G4RotationMatrix* xRot90 = new G4RotationMatrix;  
	xRot90->rotateX(M_PI/2.*rad);
	// Rotate 90 degrees around X then Z axis (Y axis in object coordinates)
	G4RotationMatrix* xzRot90 = new G4RotationMatrix;  
	xzRot90->rotateX(M_PI/2.*rad);
	xzRot90->rotateY(M_PI/2.*rad);
	
	G4VisAttributes* visCathWires = new G4VisAttributes(G4Colour(1,0.7,0,0.8));
	G4VisAttributes* visAnodeWires = new G4VisAttributes(G4Colour(1,0.3,0,0.8));
	
	// anode, cathode segments
	cathSeg_log = new G4LogicalVolume(cathodeSegmentBox,fMWPCGas,sideSubst("cathSeg_log%c",s));
	cathSeg_log->SetVisAttributes(G4VisAttributes::Invisible);
	anodeSeg_log = new G4LogicalVolume(anodeSegmentBox,fMWPCGas,sideSubst("anodeSeg_log%c",s));
	anodeSeg_log->SetVisAttributes(G4VisAttributes::Invisible);
	cathode_wire_log = new G4LogicalVolume(cathode_tube,fCathodeWireMaterial,sideSubst("cathode_log%c",s));
	cathode_wire_log->SetVisAttributes(visCathWires);
	cath_plate_log = new G4LogicalVolume(cathplate_tube,fCathodePlateMaterial,sideSubst("cathplate_log%c",s));
	cath_plate_log->SetVisAttributes(visCathWires);
	anode_wire_log = new G4LogicalVolume(anode_tube,fAnodeWireMaterial,sideSubst("anode_log%c",s));				
	anode_wire_log->SetVisAttributes(visAnodeWires);
	
	new G4PVPlacement(NULL,G4ThreeVector(),cathode_wire_log,
					  sideSubst("cathode_wire_phys%c",s),cathSeg_log,true,0);
	new G4PVPlacement(NULL,G4ThreeVector(),cath_plate_log,
					  sideSubst("cath_plate_phys%c",s),cathSeg_log,true,0);
	new G4PVPlacement(NULL,G4ThreeVector(),anode_wire_log,
					  sideSubst("anode_wire_phys%c",s),anodeSeg_log,true,0);
	
	
	// anode, cathode plane container volumes
	G4LogicalVolume* cathContainer1_log = new G4LogicalVolume(cathContainer_box,fMWPCGas,sideSubst("cathContainer1_log%c",s));
	G4LogicalVolume* cathContainer2_log = new G4LogicalVolume(cathContainer_box,fMWPCGas,sideSubst("cathContainer2_log%c",s));
	G4LogicalVolume* anodContainer_log = new G4LogicalVolume(anodeContainer_box,fMWPCGas,sideSubst("anodContainer_log%c",s));
	
	new G4PVPlacement(xRot90,G4ThreeVector(0.,0.,cathode_R-planeSpacing),
					  cathContainer1_log,sideSubst("cathContainer1_phys%c",s),gas_log,false,0);
	new G4PVPlacement(xzRot90,G4ThreeVector(0.,0.,planeSpacing-cathode_R),
					  cathContainer2_log,sideSubst("cathContainer2_phys%c",s),gas_log,false,0);
	new G4PVPlacement(xRot90,G4ThreeVector(0.,0.,0.),
					  anodContainer_log,sideSubst("anodContainer_phys%c",s),gas_log,false,0);
	
	// replicate segments into cathode, anode arrays
	new G4PVReplica(sideSubst("cathode_array_1_%c",s),
					cathSeg_log,
					cathContainer1_log,
					kXAxis,
					NbOfWires,
					spacing);
	
	new G4PVReplica(sideSubst("cathode_array_2_%c",s),
					cathSeg_log,
					cathContainer2_log,
					kXAxis,
					NbOfWires,
					spacing);
	
	new G4PVReplica(sideSubst("anode_array_%c",s),
					anodeSeg_log,
					anodContainer_log,
					kXAxis,
					NbOfWires,
					spacing);
}
