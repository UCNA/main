#ifndef bmWirechamberEMField_HH
#define bmWirechamberEMField_HH 1

#include "globals.hh"
#include "G4MagneticField.hh"
#include "G4ElectroMagneticField.hh"
#include "G4ThreeVector.hh"

using namespace std;

/// class for describing magnetic field
class bmWirechamberEMField: public G4ElectroMagneticField {
public:
	/// constructor
	bmWirechamberEMField(G4MagneticField* b);
	/// get field at given point
	void GetFieldValue(const G4double Point[4], G4double* Bfield) const;
	/// whether the field changes particle energy (true for electric field)
	virtual G4bool DoesFieldChangeEnergy() { return anodeVoltage != 0; }
private:
	G4double anodeVoltage;		//< anode plane voltage
	G4MagneticField* bField;	//< magnetic field in volume
};

#endif
