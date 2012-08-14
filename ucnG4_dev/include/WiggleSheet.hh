#ifndef WIGGLESHEET_HH
#define WIGGLESHEET_HH 1

#include "G4LogicalVolume.hh"
#include "G4VSensitiveDetector.hh"
#include "G4VisAttributes.hh"
#include <vector>

/// build a wavy sheet of materials
class WiggleSheet {
public:
	/// constructor
	WiggleSheet(double tm = M_PI/4., double pd = 2*mm, unsigned int n = 75, double l = 75*mm):
	thetaMax(tm), period(pd), nseg(n), length(l), containerMat(NULL), ttotal(0) {}
	/// add material layer of given thickness (optionally "stretched" to preserve equivalent flat thickness)
	void addLayer(G4Material* m, double t, bool stretched = false, G4VisAttributes* va = NULL);
	/// build geometry
	void Construct();
	/// set sensitive detector for all volumes
	void SetSensitiveDetector(G4VSensitiveDetector* SD);
	/// get ratio of sheet length to straight line length
	double getStretch() const;
	/// get constructed thickness of container volume
	double getContainerThick() const { return containerThick; }
	
	double thetaMax;	//< maximum wiggle angle, radians
	double period;		//< wiggle period
	unsigned int nseg;	//< number of wiggle segments
	double length;		//< sheet length
	G4Material* containerMat;		//< material for container volume
	G4LogicalVolume* container_log; 		//< overall container volume
	std::vector<G4LogicalVolume*> zigs_log;	//< segments for each material
	std::vector<G4LogicalVolume*> zags_log;	//< segments for each material
	
protected:
	std::vector<G4Material*> layerMat;	//< layer materials
	std::vector<G4VisAttributes*> layerVis;	//< layer visualization attributes
	std::vector<double> layerThick;		//< layer thicknesses
	double ttotal;						//< total thickness
	double containerThick;				//< thickness of container volume
};


#endif
