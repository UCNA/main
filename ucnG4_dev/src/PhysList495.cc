#include "PhysList495.hh"

#include <G4SystemOfUnits.hh>
#include <G4EmLivermorePhysics.hh>
#include <G4EmPenelopePhysics.hh>

#include <G4Gamma.hh>
#include <G4Electron.hh>
#include <G4Positron.hh>

#include <G4LossTableManager.hh>
#include <G4EmConfigurator.hh>
#include <G4UnitsTable.hh>

#include <G4ProcessManager.hh>

PhysList495::PhysList495(bool usePenelope) : G4VModularPhysicsList() {
	
	G4LossTableManager::Instance();
	defaultCutValue = 1.*um;
	cutForGamma     = defaultCutValue;
	cutForElectron  = 0.5*defaultCutValue;
	cutForPositron  = defaultCutValue;
	
	SetVerboseLevel(1);
	
	// EM physics
	if(usePenelope) {
		emName = G4String("Penelope");  
		emPhysicsList = new G4EmPenelopePhysics();
	} else {
		emName = G4String("Livermore");  
		emPhysicsList = new G4EmLivermorePhysics();
	}
}

PhysList495::~PhysList495() {
	delete emPhysicsList;
}

////////////////////////////////////////////////////////////////////////////
// Construct Particles /////////////////////////////////////////////////////

void PhysList495::ConstructParticle() {	
	emPhysicsList->ConstructParticle();
}

void PhysList495::ConstructProcess() {
	// transportation process
	AddTransportation();
	// electromagnetic physics list
	emPhysicsList->ConstructProcess();
}

/*
void Livermore_EMPhysList::AddStepMax()
{
	// Step limitation seen as a process
	stepMaxProcess = new StepMax();
	
	theParticleIterator->reset();
	while ((*theParticleIterator)()){
		G4ParticleDefinition* particle = theParticleIterator->value();
		G4ProcessManager* pmanager = particle->GetProcessManager();
		if (stepMaxProcess->IsApplicable(*particle) && pmanager) {
			pmanager->AddDiscreteProcess(stepMaxProcess);
		}
	}
}
*/

void PhysList495::SetCuts() {
	
	if (verboseLevel >0) {
		G4cout << "PhysicsList::SetCuts:";
		G4cout << "CutLength : " << G4BestUnit(defaultCutValue,"Length") << G4endl;
	}
	
	if(emName == "Livermore")
		G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(250*eV, 1*GeV);
	else if(emName == "Penelope")
		G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(100*eV, 1*GeV);

	// set cut values for gamma at first and for e- second and next for e+,
	// because some processes for e+/e- need cut values for gamma
	SetCutValue(cutForGamma, "gamma");
	SetCutValue(cutForElectron, "e-");
	SetCutValue(cutForPositron, "e+");
	
	if (verboseLevel>0) DumpCutValuesTable();
}

void PhysList495::SetCutForGamma(G4double cut) {
	cutForGamma = cut;
	SetParticleCuts(cutForGamma, G4Gamma::Gamma());
}

void PhysList495::SetCutForElectron(G4double cut) {
	cutForElectron = cut;
	SetParticleCuts(cutForElectron, G4Electron::Electron());
}

void PhysList495::SetCutForPositron(G4double cut) {
	cutForPositron = cut;
	SetParticleCuts(cutForPositron, G4Positron::Positron());
}
