#include "PrimaryGeneratorMessenger.hh"
#include "PrimaryGeneratorAction.hh"

#include <G4UIdirectory.hh>
#include <G4UIcmdWithAString.hh>

PrimaryGeneratorMessenger::PrimaryGeneratorMessenger(PrimaryGeneratorAction* Gun): Action(Gun) {
	gunDir = new G4UIdirectory("/benchmark/gun/");
	gunDir->SetGuidance("PrimaryGenerator control");
	
	gunTypeCmd = new G4UIcmdWithAString("/benchmark/gun/type",this);
	gunTypeCmd->SetGuidance("Set the generator gun type.");
	gunTypeCmd->SetGuidance(" Choices: eGun eGunRandMomentum");
	gunTypeCmd->SetDefaultValue("eGun");
	gunTypeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
	
	gunPtclCmd = new G4UIcmdWithAString("/benchmark/gun/particle",this);
	gunPtclCmd->SetGuidance("Set the gun particle thrown.");
	gunPtclCmd->SetDefaultValue("e-");
	gunPtclCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
	
	positionerCmd = new G4UIcmdWithAString("/benchmark/gun/positioner",this);
	positionerCmd->SetGuidance("Set the generator gun positioner.");
	positionerCmd->SetGuidance(" Choice : Fixed, SourceDrop, DecayTrapUniform, DecayTrapFiducial");
	positionerCmd->SetDefaultValue("Fixed");
	positionerCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
	
	srcrelCmd = new G4UIcmdWithABool("/benchmark/gun/relholder",this);
	srcrelCmd->SetGuidance("Set event positioning relative to source holder location");
	srcrelCmd->SetDefaultValue(false);
	
	eventFileCmd = new G4UIcmdWithAString("/benchmark/gun/evtfile",this);
	eventFileCmd->SetGuidance("Set input file for events");
	eventFileCmd->SetDefaultValue("");
	eventFileCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

	sourceRadiusCmd = new G4UIcmdWithADoubleAndUnit("/benchmark/gun/sourceRadius",this);
	sourceRadiusCmd->SetGuidance("Radius for SourceDrop generator");
	sourceRadiusCmd->SetDefaultValue(1.5*mm);
	sourceRadiusCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger() {
	delete gunTypeCmd;
	delete positionerCmd;
	delete sourceRadiusCmd;
	delete gunDir;
}

void PrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command, G4String newValue) { 
	if( command == gunTypeCmd )
		Action->SetGunType(newValue);
	if( command == gunPtclCmd )
		Action->SetParticleType(newValue);
	if( command == positionerCmd )
		Action->SetPositioner(newValue);
	if( command == sourceRadiusCmd )
		Action->SetSourceRadius(sourceRadiusCmd->GetNewDoubleValue(newValue));
	if( command == srcrelCmd )
		Action->SetPosRelHolder(srcrelCmd->GetNewBoolValue(newValue));
	if( command == eventFileCmd )
		Action->SetEventFile(newValue);
}
