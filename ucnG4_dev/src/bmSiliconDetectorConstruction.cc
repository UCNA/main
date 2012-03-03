#include "bmSiliconDetectorConstruction.hh"
#include <cassert>

void bmSiliconDetectorConstruction::Construct() {
	
	assert(pDetMat);
	assert(pHolderMat);
	assert(pBackingMat);	
	assert(fHolderThick >= fDetThick + fBackingThick);
	assert(fHolderRadius > fDetRadius);
	
	// container volume for whole assembly
	G4Tubs* container_tube = new G4Tubs("si_container_tube",0,fHolderRadius,fHolderThick/2.,0.,2*M_PI);
	container_log = new G4LogicalVolume(container_tube, Vacuum, "si_container_log");
	
	// detector holder
	G4Tubs* holder_tube = new G4Tubs("si_holder_tube",fDetRadius,fHolderRadius,fHolderThick/2.,0.,2*M_PI);
	holder_log = new G4LogicalVolume(holder_tube, pHolderMat, "si_holder_log");
	holder_phys = new G4PVPlacement(NULL,G4ThreeVector(),holder_log,"si_holder_phys",container_log,true,0);
	G4VisAttributes* visHolder = new G4VisAttributes(G4Colour(1.0,0.7,0.0,1.0));
	holder_log->SetVisAttributes(visHolder);
	
	// active detector
	G4Tubs* det_tube = new G4Tubs("si_det_tube",0,fDetRadius,fDetThick/2.,0.,2*M_PI);
	det_log = new G4LogicalVolume(det_tube, pDetMat, "si_det_log");
	det_phys = new G4PVPlacement(NULL,G4ThreeVector(0,0,fHolderThick*0.5-fBackingThick-fDetThick*0.5),
								 det_log,"si_det_phys",container_log,true,0);
	G4VisAttributes* visDet = new G4VisAttributes(G4Colour(0.7,0.7,0.7,1.0));
	det_log->SetVisAttributes(visDet);

	// detector backing
	G4Tubs* backing_tube = new G4Tubs("si_backing_tube",0,fDetRadius,fBackingThick/2.,0.,2*M_PI);
	backing_log = new G4LogicalVolume(backing_tube, pBackingMat, "si_backing_log");
	backing_phys = new G4PVPlacement(NULL,G4ThreeVector(0,0,fHolderThick*0.5-fBackingThick*0.5),
									 backing_log,"si_backing_phys",container_log,true,0);
	G4VisAttributes* visBack = new G4VisAttributes(G4Colour(0.2,0.2,0.2,1.0));
	backing_log->SetVisAttributes(visBack);
}
