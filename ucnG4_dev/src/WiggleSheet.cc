#include "WiggleSheet.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4PVPlacement.hh"
#include "SMExcept.hh"
#include <cmath>

void WiggleSheet::addLayer(G4Material* m, double t, bool stretched, G4VisAttributes* va) {
	if(!t) return;
	assert(m);
	if(stretched) t /= getStretch();
	layerMat.push_back(m);
	layerThick.push_back(t);
	layerVis.push_back(va);
	ttotal += t;
}

//-----------------
//    \\	.
// ^   ||	zig
// |  //	.
// | ||		zag
// y  \\	.
// zx---->
//-----------------

void WiggleSheet::Construct() {
	double rcenter = period/(4.*sin(thetaMax));	// center radius
	double rmax = rcenter+ttotal/2.;			// outer radius
	double cmax = rmax-sqrt(rcenter*rcenter-period*period/16.);	// cord distance to outer radius
	containerThick = 2*cmax;
	
	if(!containerMat) throw(SMExcept("ContainerMaterialUnspecified"));
	
	G4Box* container_box = new G4Box("container_box",cmax,(nseg*period/2+ttotal)/2.,length/2.);
	container_log = new G4LogicalVolume(container_box, containerMat, "wigglesheet_container_log");
	container_log->SetVisAttributes(G4VisAttributes::Invisible);

	// construct logical volumes
	double rout = rmax;
	double rin = rmax-ttotal;
	for(unsigned int i=0; i<layerMat.size(); i++) {
		double t = layerThick[i];
		G4Tubs* zig_tube = new G4Tubs("zig_tube",rout-t,rout,length/2.,-thetaMax,2*thetaMax);
		G4Tubs* zag_tube = new G4Tubs("zag_tube",rin,rin+t,length/2.,M_PI-thetaMax,2*thetaMax);
		zigs_log.push_back(new G4LogicalVolume(zig_tube, layerMat[i], "zig_log"));
		zags_log.push_back(new G4LogicalVolume(zag_tube, layerMat[i], "zag_log"));
		if(layerVis[i]) {
			zigs_log.back()->SetVisAttributes(layerVis[i]);
			zags_log.back()->SetVisAttributes(layerVis[i]);
		}
		rout -= t;
		rin += t;
	}
	
	//G4cout << "WiggleSheet for theta=" << thetaMax << " pd=" << period << " nseg=" << nseg <<
	//" length=" << length << " cthick=" << getContainerThick() << G4endl;
	
	// position physical volumes
	for(unsigned int n=0; n<nseg; n++) {
		std::vector<G4LogicalVolume*>& logs = (n%2)?zags_log:zigs_log;
		double ypos = (int(nseg)-1-2*int(n))*period/4.;
		double xpos = (rmax-cmax)*(n%2?1:-1);
		for(unsigned int i=0; i<logs.size(); i++)
			new G4PVPlacement(NULL,G4ThreeVector(xpos,ypos,0),logs[i],"seg_phys",container_log,true,n/2);
	}
}

double WiggleSheet::getStretch() const {
	return thetaMax/fabs(sin(thetaMax));
}

void WiggleSheet::SetSensitiveDetector(G4VSensitiveDetector* SD) {
	for(unsigned int i=0; i<zigs_log.size(); i++) {
		zigs_log[i]->SetSensitiveDetector(SD);
		zags_log[i]->SetSensitiveDetector(SD);
	}
}
