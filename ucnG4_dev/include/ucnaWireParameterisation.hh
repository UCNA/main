//Follow the example of novice/N02/ExN02ChamberParameterisation.
//--J.Yuan
///////////////////////////////////////////////////////////////

#ifndef ucnaWireParameterisation_H
#define ucnaWireParameterisation_H 1

#include "globals.hh"
#include "G4VPVParameterisation.hh"
#include "G4RotationMatrix.hh"

class G4VPhysicalVolume;
class G4Tubs;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class ucnaWireParameterisation : public G4VPVParameterisation
{ 
  public:
  
  ucnaWireParameterisation(G4int    NoWires, 
			   G4double Zpos,
			   G4double startPos,
			   G4double spacing,
			   G4double tubeR, 
			   G4double maxLength,
			   G4int xyz);
  
  virtual ~ucnaWireParameterisation();
  
  void ComputeTransformation (const G4int copyNo,
			      G4VPhysicalVolume* physVol) const;
  
  void ComputeDimensions (G4Tubs & wireTube, const G4int copyNo,
				  const G4VPhysicalVolume* physVol) const;
  // Dummy declarations to get rid of warnings ...
  
  void ComputeDimensions (G4Box&,const G4int,const G4VPhysicalVolume*) const {}
  
  void ComputeDimensions (G4Trd&,const G4int,const G4VPhysicalVolume*) const {}
  void ComputeDimensions (G4Trap&,const G4int,const G4VPhysicalVolume*) const 
  {}
  void ComputeDimensions(G4Cons&, const G4int, const G4VPhysicalVolume*) const 
  {}
  void ComputeDimensions (G4Orb&,const G4int,const G4VPhysicalVolume*) const {}
  void ComputeDimensions (G4Sphere&,const G4int,const G4VPhysicalVolume*) const {}
  void ComputeDimensions (G4Torus&,const G4int,const G4VPhysicalVolume*) const {}
  void ComputeDimensions (G4Para&,const G4int,const G4VPhysicalVolume*) const
  {}
  void ComputeDimensions (G4Hype&,const G4int,const G4VPhysicalVolume*) const
  {}
  void ComputeDimensions (G4Polycone&,const G4int,const G4VPhysicalVolume*) const {}
  void ComputeDimensions (G4Polyhedra&,const G4int,const G4VPhysicalVolume*) const {}
  
  
private:
  
  G4int    fNoWires;
  G4double fStartPos;  //start position on X or Y.
  G4double fZpos;  //Z positon.
  G4double fRadius;    //  Radius of the tube
  G4double fSpacing;      //  The distance between the wires' center
  G4double fMaxLength;   //  The Increment for the half-length 
  G4RotationMatrix* fRM;
  G4int fXYZ;   //rotation on X or Y?
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
