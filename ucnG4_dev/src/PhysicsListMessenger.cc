#include "PhysicsListMessenger.hh"

PhysicsListMessenger::PhysicsListMessenger(PhysList495* pPhys) :fPhysList(pPhys) {
  fPhysDir = new G4UIdirectory("/phys/");
  fPhysDir->SetGuidance("physics list commands");
   
  fListCmd = new G4UIcmdWithAString("/phys/physlist",this);
  fListCmd->SetGuidance("Set EM physics list.");
  fListCmd->SetDefaultValue("Livermore");
  fListCmd->AvailableForStates(G4State_PreInit);
}

PhysicsListMessenger::~PhysicsListMessenger() {
  delete fListCmd;
  delete fPhysDir;
}

void PhysicsListMessenger::SetNewValue(G4UIcommand* command, G4String newValue) {
  if( command == fListCmd ) fPhysList->setPhysicsList(newValue);
}
