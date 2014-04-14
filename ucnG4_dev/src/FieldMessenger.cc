#include "FieldMessenger.hh"

#include <G4SystemOfUnits.hh>

FieldMessenger::FieldMessenger(Field* f): myField(f) {
	
	fFieldDir = new G4UIdirectory("/field/");
	fFieldDir->SetGuidance("Magnetic field settings");
	
	fFieldScaleCmd = new G4UIcmdWithADouble("/field/scale",this);
	fFieldScaleCmd->SetGuidance("Scale factor for main magnetic field strength");
	fFieldScaleCmd->SetDefaultValue(1.0);
	fFieldScaleCmd->AvailableForStates(G4State_Idle);

	fFieldMapFileCmd = new G4UIcmdWithAString("/field/mapfile",this);
	fFieldMapFileCmd->SetParameterName("fname",true);
	fFieldMapFileCmd->SetGuidance("Set B field map file");
	fFieldMapFileCmd->SetDefaultValue("");
	fFieldMapFileCmd->AvailableForStates(G4State_Idle);
	
	fAfpDipoleCmd = new G4UIcmdWithADouble("/field/AFP_dipole",this);
	fAfpDipoleCmd->SetGuidance("AFP dipole strength in A*m^2");
	fAfpDipoleCmd->SetDefaultValue(0.);
	fAfpDipoleCmd->AvailableForStates(G4State_Idle);
	
	fCheckFieldCmd = new G4UIcmdWith3VectorAndUnit("/field/check",this);
	fCheckFieldCmd->SetGuidance("Check field value at given position");
	fCheckFieldCmd->AvailableForStates(G4State_Idle);
}

FieldMessenger::~FieldMessenger() {
	delete fFieldScaleCmd;
	delete fFieldMapFileCmd;
	delete fAfpDipoleCmd;
	delete fFieldDir;
}

void FieldMessenger::SetNewValue(G4UIcommand * command, G4String newValue) {
	if(command == fFieldScaleCmd) myField->SetFieldScale(fFieldScaleCmd->GetNewDoubleValue(newValue));
	else if(command == fFieldMapFileCmd) myField->LoadFieldMap(newValue);
	else if(command == fAfpDipoleCmd) myField->SetAFPDipole(fAfpDipoleCmd->GetNewDoubleValue(newValue));
	else if(command == fCheckFieldCmd) {
		G4ThreeVector x = fCheckFieldCmd->GetNew3VectorValue(newValue);
		G4ThreeVector B;
		myField->GetFieldValue(&x[0], &B[0]);
		G4cout << " B = " << B[0]/tesla << " " << B[1]/tesla << " " << B[2]/tesla << " T" << G4endl;
	}
}
