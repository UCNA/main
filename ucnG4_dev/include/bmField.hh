//2-19-02 J. Yuan: magnet field for UCNA experiment
////////////////////////////////////////////////////

#ifndef bmField_H
#define bmField_H 1

#include "globals.hh"
#include "G4MagneticField.hh"
#include <TString.h>
#include "G4ThreeVector.hh"
#include <vector>

using namespace std;

/// class for describing magnetic field
class bmField: public G4MagneticField {
public:
	/// constructor
	bmField();
	/// constructor with fieldmap
	bmField(const TString filename);
	/// get field at given point
	void GetFieldValue( const  G4double Point[3], G4double *Bfield ) const;
	/// zero out the field (via scaling factor)
	void SetFieldToZero() { fieldScale=0; };
	/// set fieldmap scaling factor
	void SetFieldScale(const G4double val) { fieldScale=val; }
	/// load fieldmap from file
	void LoadFieldMap(const TString filename);
	
	bool addAFP;	//< whether to add the AFP fringe field
	
private:
	/// add point to field profile
	void addPoint(G4double z, G4double B) { Zpoints.push_back(z); Bpoints.push_back(B); }
	vector<G4double> Bpoints;	//< field profile B values
	vector<G4double> Zpoints;	//< field profile z positions
	G4double rmax2;				//< max radius squared (position in world volume) to apply field
	G4double fieldScale;		//< scaling factor for field strength
};

#endif

