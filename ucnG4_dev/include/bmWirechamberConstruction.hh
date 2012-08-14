#ifndef bmWirechamberConstruction_HH
#define bmWirechamberConstruction_HH 1

#include "G4ElectroMagneticField.hh"
#include "G4MagneticField.hh"
#include "G4RotationMatrix.hh"

#include "bmDetectorConstructionUtils.hh"
#include "bmWireVolumeConstruction.hh"

/// class for constructing entire wirechamber
class bmWirechamberConstruction: public MaterialUser, G4ElectroMagneticField {
public:
	/// constructor
	bmWirechamberConstruction(): fWindowThick(6*um), mwpc_entrance_R(7.0*cm), mwpc_exit_R(7.5*cm),
	fMWPCGas(WCPentane), entranceToCathodes(5.0*mm), exitToCathodes(5.0*mm), myBField(NULL), E0(0) {}
	
	/// get constructed width
	G4double GetWidth() const { return 2*mwpcContainer_halfZ; }
	
	G4double fWindowThick;			//< mylar window thickness
	G4double mwpc_entrance_R;		//< entrance window radius
	G4double mwpc_exit_R;			//< exit window radius
	G4Material* fMWPCGas;			//< MWPC fill gas
	G4double entranceToCathodes;	//< entrance-window-to-cathode distance
	G4double exitToCathodes;		//< exit-window-to-cathode distance
	
	bmWireVolumeConstruction activeRegion;	//< active gas region with wireplanes
	
	G4LogicalVolume* container_log;	//< overall gas box
	G4LogicalVolume* winIn_log;		//< inner window
	G4LogicalVolume* winOut_log;	//< outer window
		
	G4LogicalVolume* kevContainer_log;	//< container volume for kevlar strip array
	G4LogicalVolume* kevSeg_log;		//< one segment of kevlar strip array
	G4LogicalVolume* kevStrip_log;		//< kevlar strip in one segment
	
	/// construct logical container volume
	void Construct(Side s);
	
	/// electromagnetic field
	virtual void  GetFieldValue(const G4double Point[4], G4double* Bfield) const;
	/// whether the field changes particle energy
	virtual G4bool DoesFieldChangeEnergy() const { return E0 != 0; }
	/// set up tracking in field
	void ConstructField();
	/// set anode voltage
	void setPotential(G4double Vanode);
	
	G4MagneticField* myBField;		//< Magnetic field pointer
	G4RotationMatrix* myRotation;	//< rotation from global frame to local coordinates
	G4ThreeVector myTranslation;	//< translation from global coordinates to center of anode plane	
	
protected:
	G4double mwpcContainer_halfZ;	//< half-width of wirechamber
	G4double E0;					//< field scaling constant
	double d;						//< spacing between anode wires
	double L;						//< spacing between planes
	double r;						//< anode wire radius
	
};

#endif
