//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************

/*
#include "bmPhysicsList.hh"

#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"

#include "globals.hh"
#include "G4ios.hh"
#include <iomanip>


// Constructor /////////////////////////////////////////////////////////////
bmPhysicsList::bmPhysicsList() : G4VUserPhysicsList(), defaultCutValue(1.0*mm) {
	
	SetVerboseLevel(1);
		
	physDir = new G4UIdirectory("/bm/phys/");
	physDir->SetGuidance("physics list commands");
	
	allCutCmd = new G4UIcmdWithADoubleAndUnit("/bm/phys/setCuts",this);  
	allCutCmd->SetGuidance("Set cut for all.");
	allCutCmd->SetParameterName("cut",false);
	allCutCmd->SetUnitCategory("Length");
	allCutCmd->SetRange("cut>0.0");
	allCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  
	
}


////////////////////////////////////////////////////////////////////////////
// Construct Particles /////////////////////////////////////////////////////

#include "G4LeptonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4BosonConstructor.hh"

void bmPhysicsList::ConstructParticle() {
	
	G4LeptonConstructor lepton;
	lepton.ConstructParticle();
	
	G4BosonConstructor boson;
	boson.ConstructParticle();

	G4BaryonConstructor baryon;
	baryon.ConstructParticle();
}

/////////////////////////////////////////////////////////////////////////////
// Construct Processes //////////////////////////////////////////////////////

void bmPhysicsList::ConstructProcess() {
	AddTransportation();
	StandardEM();		
}

/////////////////////////////////////////////////////////////////////////////
// Electromagnetic Processes ////////////////////////////////////////////////

// all charged particles
#include "G4eMultipleScattering.hh"
#include "G4MuMultipleScattering.hh"
#include "G4hMultipleScattering.hh"

// gamma
#include "G4LowEnergyRayleigh.hh" 
#include "G4LowEnergyPhotoElectric.hh"
#include "G4LowEnergyCompton.hh"  
#include "G4LowEnergyGammaConversion.hh" 

// e-
#include "G4LowEnergyIonisation.hh" 
#include "G4LowEnergyBremsstrahlung.hh" 

// e+
#include "G4eIonisation.hh" 
#include "G4eBremsstrahlung.hh" 
#include "G4eplusAnnihilation.hh"

// muons
#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"
#include "G4MuonMinusCaptureAtRest.hh"

// hadrons and ions
#include "G4hIonisation.hh"
#include "G4ionIonisation.hh"

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4StepLimiter.hh"

void bmPhysicsList::StandardEM()
{
	// Add standard EM Processes
	
	theParticleIterator->reset();
	while( (*theParticleIterator)() ){
		G4ParticleDefinition* particle = theParticleIterator->value();
		G4ProcessManager* pmanager = particle->GetProcessManager();
		G4String particleName = particle->GetParticleName();
		
		if (particleName == "gamma") {
			// gamma
			pmanager->AddDiscreteProcess(new G4PhotoElectricEffect);
			pmanager->AddDiscreteProcess(new G4ComptonScattering);
			pmanager->AddDiscreteProcess(new G4GammaConversion);
			pmanager->AddProcess(new G4StepLimiter,       -1,-1,3);
		} else if (particleName == "e-") {
			//electron
			G4eMultipleScattering *msc = new G4eMultipleScattering();
			//msc->SetFacrange(0.001);
			pmanager->AddProcess(msc, -1, 1,1);
			pmanager->AddProcess(new G4eIonisation,        -1, 2,2);
			pmanager->AddProcess(new G4eBremsstrahlung,    -1, 3,3);
			pmanager->AddProcess(new G4StepLimiter,       -1,-1,3);
		} else if (particleName == "e+") {
			//positron
			pmanager->AddProcess(new G4eMultipleScattering, -1, 1,1);
			pmanager->AddProcess(new G4eIonisation,        -1, 2,2);
			pmanager->AddProcess(new G4eBremsstrahlung,    -1, 3,3);
			pmanager->AddProcess(new G4eplusAnnihilation,   0,-1,4);
			pmanager->AddProcess(new G4StepLimiter,       -1,-1,3);
			
		} else if (particleName == "mu+" || 
				   particleName == "mu-"    ) {
			//muon  
			pmanager->AddProcess(new G4MuMultipleScattering, -1, 1,1);
			pmanager->AddProcess(new G4MuIonisation,       -1, 2,2);
			pmanager->AddProcess(new G4MuBremsstrahlung,   -1, 3,3);
			pmanager->AddProcess(new G4MuPairProduction,   -1, 4,4);
			
		} else if (particleName == "alpha" || particleName == "GenericIon" ) { 
			pmanager->AddProcess(new G4hMultipleScattering, -1, 1,1);
			pmanager->AddProcess(new G4ionIonisation,      -1, 2,2);
			
		} else if ((!particle->IsShortLived()) &&
				   (particle->GetPDGCharge() != 0.0) && 
				   (particle->GetParticleName() != "chargedgeantino")) {
			//all others charged particles except geantino
			pmanager->AddProcess(new G4hMultipleScattering, -1,1,1);
			pmanager->AddProcess(new G4hIonisation,        -1,2,2);
		}
	}
}


void bmPhysicsList::SetCuts() {
	
	if (GetVerboseLevel() >0){
		G4cout << "PhysicsList::SetCuts:";
		G4cout << "CutLength : " << G4BestUnit(defaultCutValue,"Length") << G4endl;
	}
	
	// set cut values for gamma at first and for e- second and next for e+,
	// because some processes for e+/e- need cut values for gamma
	cutForGamma     = defaultCutValue;
	cutForElectron  = defaultCutValue;
	cutForPositron  = defaultCutValue;
	SetCutValue(cutForGamma, "gamma");
	SetCutValue(cutForElectron, "e-");
	SetCutValue(cutForPositron, "e+");
	
	if (GetVerboseLevel()>0) DumpCutValuesTable();
}

void bmPhysicsList::SetNewValue(G4UIcommand* command, G4String newValue) { 
	if( command == allCutCmd ) {
		defaultCutValue = allCutCmd->GetNewDoubleValue(newValue);
		SetCuts();
    }
}
*/
