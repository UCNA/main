#include "bmLivermore_EMPhysList.hh"

#include "G4ProcessManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4ProductionCutsTable.hh"

#include "globals.hh"
#include "G4UnitsTable.hh"

////////////////////////////////////////////////////////////////////////////
// Construct Particles /////////////////////////////////////////////////////

#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Proton.hh"
#include "G4AntiProton.hh"
#include "G4KaonMinus.hh"
#include "G4IonConstructor.hh"

void bmLivermore_EMPhysList::ConstructParticle() {	
	G4Gamma::GammaDefinition();
	G4Electron::ElectronDefinition();
	G4Positron::PositronDefinition();
	G4Proton::ProtonDefinition();
	G4AntiProton::AntiProtonDefinition();
	G4KaonMinus::KaonMinusDefinition();
	G4IonConstructor ions;
	ions.ConstructParticle();
}

/////////////////////////////////////////////////////////////////////////////
// Construct Processes //////////////////////////////////////////////////////

void bmLivermore_EMPhysList::ConstructProcess() {
	AddTransportation();
	StandardEM();		
}


// production cuts
void bmLivermore_EMPhysList::SetCuts() {
	G4double defaultCutValue = 0.1*um;
	if (GetVerboseLevel() >0){
		G4cout << "PhysicsList::SetCuts:";
		G4cout << "CutLength : " << G4BestUnit(defaultCutValue,"Length") << G4endl;
	}
	SetCutValue(defaultCutValue, "gamma");
	SetCutValue(defaultCutValue, "e-");
	SetCutValue(defaultCutValue, "e+");
	G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(250*eV, 1*GeV);
	if (GetVerboseLevel()>0) DumpCutValuesTable();
}


/////////////////////////////////////////////////////////////////////////////
// Electromagnetic Processes, Penelope model ////////////////////////////////

#include "G4StepLimiter.hh"

// gamma
#include "G4PhotoElectricEffect.hh" 
#include "G4LivermorePhotoElectricModel.hh" 
#include "G4ComptonScattering.hh" 
#include "G4LivermoreComptonModel.hh" 
#include "G4GammaConversion.hh"
#include "G4LivermoreGammaConversionModel.hh"
#include "G4RayleighScattering.hh"
#include "G4LivermoreRayleighModel.hh"

// e-,e+
#include "G4eMultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4LivermoreIonisationModel.hh"
#include "G4eBremsstrahlung.hh" 
#include "G4LivermoreBremsstrahlungModel.hh"
#include "G4eplusAnnihilation.hh"

// other charged particles
#include "G4hMultipleScattering.hh"
#include "G4hIonisation.hh"
#include "G4ionIonisation.hh"
#include "G4IonParametrisedLossModel.hh"

void bmLivermore_EMPhysList::StandardEM() {
	
	theParticleIterator->reset();
	
	while( (*theParticleIterator)() ){
		
		G4ParticleDefinition* particle = theParticleIterator->value();
		G4ProcessManager* pmanager = particle->GetProcessManager();
		G4String particleName = particle->GetParticleName();
		
		if (particleName == "gamma") {
			// gammas
			
			G4PhotoElectricEffect* thePhotoElectricEffect = new G4PhotoElectricEffect();
			thePhotoElectricEffect->SetModel(new G4LivermorePhotoElectricModel());
			pmanager->AddDiscreteProcess(thePhotoElectricEffect);
			
			G4ComptonScattering* theComptonScattering = new G4ComptonScattering();
			theComptonScattering->SetModel(new G4LivermoreComptonModel());
			pmanager->AddDiscreteProcess(theComptonScattering);
			
			G4GammaConversion* theGammaConversion = new G4GammaConversion();
			theGammaConversion->SetModel(new G4LivermoreGammaConversionModel());
			pmanager->AddDiscreteProcess(theGammaConversion);
			
			G4RayleighScattering* theRayleigh = new G4RayleighScattering();
			theRayleigh->SetModel(new G4LivermoreRayleighModel());
			pmanager->AddDiscreteProcess(theRayleigh);
			
			pmanager->AddProcess(new G4StepLimiter,       -1,-1,3);
		} else if (particleName == "e-" || particleName == "e+") {
			// electrons and positrons
			pmanager->AddProcess(new G4eMultipleScattering(), -1, 1,1);
			
			G4eIonisation* theIonisation = new G4eIonisation();
			if(particleName=="e-")
				theIonisation->SetEmModel(new G4LivermoreIonisationModel());
			pmanager->AddProcess(theIonisation, -1, 2, 2);
			
			G4eBremsstrahlung* theBremsstrahlung = new G4eBremsstrahlung();
			if(particleName=="e-")
				theBremsstrahlung->SetEmModel(new G4LivermoreBremsstrahlungModel());
			pmanager->AddProcess(theBremsstrahlung,  -1,-3, 3);
			
			if(particleName=="e+")
				pmanager->AddProcess(new G4eplusAnnihilation(),   0,-1,4);
			
			pmanager->AddProcess(new G4StepLimiter,       -1,-1,3);
		} else if (particleName == "alpha" || particleName == "GenericIon" ) { 
			// ions
			pmanager->AddProcess(new G4hMultipleScattering, -1, 1,1);
			
			G4ionIonisation* theIonisation = new G4ionIonisation();
			if(particleName == "GenericIon")
				theIonisation->SetEmModel(new G4IonParametrisedLossModel());
			pmanager->AddProcess(theIonisation, -1, 2,2);
			
			pmanager->AddProcess(new G4StepLimiter, -1,-1,3);
		} else if ((!particle->IsShortLived()) &&
				   (particle->GetPDGCharge() != 0.0) && 
				   (particle->GetParticleName() != "chargedgeantino")) {
			// misc. charged
			pmanager->AddProcess(new G4hMultipleScattering, -1,1,1);
			pmanager->AddProcess(new G4hIonisation,        -1,2,2);
			pmanager->AddProcess(new G4StepLimiter,       -1,-1,3);
		}
	}
}
