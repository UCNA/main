#include "PrimaryGeneratorMessenger.hh"
#include "PrimaryGeneratorAction.hh"

#include <G4UIdirectory.hh>
#include <G4UIcmdWithAString.hh>

PrimaryGeneratorMessenger::PrimaryGeneratorMessenger(PrimaryGeneratorAction* Gun): Action(Gun) {
	gunDir = new G4UIdirectory("/generator/");
	gunDir->SetGuidance("PrimaryGenerator control");
	
	eventFileCmd = new G4UIcmdWithAString("/generator/evtfile",this);
	eventFileCmd->SetGuidance("Set input file for events");
	eventFileCmd->SetDefaultValue("");
	eventFileCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
	
	srcrelCmd = new G4UIcmdWithABool("/generator/relholder",this);
	srcrelCmd->SetGuidance("Set event positioning relative to source holder location");
	srcrelCmd->SetDefaultValue(false);
	
	sourceRadiusCmd = new G4UIcmdWithADoubleAndUnit("/generator/sourceRadius",this);
	sourceRadiusCmd->SetGuidance("Radius for SourceDrop generator");
	sourceRadiusCmd->SetDefaultValue(1.5*mm);
	sourceRadiusCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger() {
	delete eventFileCmd;
	delete srcrelCmd;
	delete sourceRadiusCmd;
	delete gunDir;
}

void PrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command, G4String newValue) {
	if( command == sourceRadiusCmd )
		Action->SetSourceRadius(sourceRadiusCmd->GetNewDoubleValue(newValue));
	if( command == srcrelCmd )
		Action->SetPosRelHolder(srcrelCmd->GetNewBoolValue(newValue));
	if( command == eventFileCmd )
		Action->SetEventFile(newValue);
}
