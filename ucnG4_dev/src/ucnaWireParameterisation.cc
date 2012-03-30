
#include "ucnaWireParameterisation.hh"

#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4Tubs.hh"
//#include "G4RotationMatrix.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ucnaWireParameterisation::ucnaWireParameterisation(
        G4int    NoWires,
	G4double Zpos,     //Z position
        G4double startPos,          // start Positon on X or Y
        G4double spacing,        //  spacing of Wire centers
        G4double tubeR,
        G4double maxLength,
	G4int xyz) //xyz = 1/2: place replication along x/y 
{
   fNoWires =  NoWires;
   fStartPos     =  startPos;
   fZpos=Zpos;
   fRadius  =  tubeR;
   fSpacing    =  spacing;
   fMaxLength = maxLength;
   fRM=new G4RotationMatrix();
   fXYZ=xyz;
   if(xyz==1) fRM->rotateX(90*deg);
   else if(xyz==2) fRM->rotateY(90*deg);
   //else G4Exception("WireParameterisation: rotation matrix needed!");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ucnaWireParameterisation::~ucnaWireParameterisation()
{
  delete fRM;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
   

void ucnaWireParameterisation::ComputeTransformation
(const G4int copyNo, G4VPhysicalVolume* physVol) const
{
  G4double      position= fStartPos + copyNo * fSpacing;
  G4ThreeVector origin;
  if(fXYZ==1)   origin.set(position,0,fZpos);  //pos in the new frame.
  else if(fXYZ==2)  origin.set(0,position,fZpos); //pos in the new frame.
  //else G4Exception("WireParameterisation:rotation?");
  physVol->SetTranslation(origin);
  physVol->SetRotation(fRM);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ucnaWireParameterisation::ComputeDimensions
(G4Tubs& wireTube, const G4int copyNo,
 const G4VPhysicalVolume* /*physVol*/) const
{
  G4double      position= fStartPos + copyNo * fSpacing;
  G4double  halfLength= sqrt(fMaxLength*fMaxLength-position*position);
  if(halfLength<0) halfLength=0;
  wireTube.SetZHalfLength(halfLength);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
