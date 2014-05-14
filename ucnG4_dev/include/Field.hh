#ifndef FIELD_HH
#define FIELD_HH

#include <vector>

#include <globals.hh>
#include <G4MagneticField.hh>
#include <G4ThreeVector.hh>

using namespace std;
class FieldMessenger;

/// class for describing magnetic field
class Field: public G4MagneticField {
public:
	/// constructor
	Field(const G4String& filename = "");
	/// get field at given point
	void GetFieldValue( const  G4double Point[3], G4double *Bfield ) const;
	/// set fieldmap scaling factor
	void SetFieldScale(G4double val) { fieldScale = val; }
	/// set AFP dipole fringe
	void SetAFPDipole(G4double val) { afp_m = val; }
	/// load fieldmap from file
	void LoadFieldMap(const G4String& filename);
	
private:
	
	FieldMessenger* myMessenger;	///< UI messenger
	
	/// add point to field profile
	void addPoint(G4double z, G4double B) { Zpoints.push_back(z); Bpoints.push_back(B); }
	vector<G4double> Bpoints;		///< field profile B values
	vector<G4double> Zpoints;		///< field profile z positions
	G4double rmax2;					///< max radius squared (position in world volume) to apply field
        
	/// add AFP fringe field contribution
	void addAFPFringeField(const G4double Point[3], G4double *Bfield) const;
	
	G4double fieldScale;			///< scaling factor for field strength
        
	G4double afp_m;					///< magnitude of AFP fringe dipole contribution, in A * m^2 (~14600 A*m^2)
        bool RobbyField;          ///< for implementing Robby's Fortran code
};

#endif

