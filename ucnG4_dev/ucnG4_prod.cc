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
//
//
// $Id: ucnG4_prod.cc,v 1.9 2011-10-12 19:31:16 mmendenhall Exp $
// GEANT4 tag $Name:  $
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "bmDetectorConstruction.hh"
#include "bmPhysicsList.hh"
#include "bmPenelope2008_EMPhysList.hh"
#include "bmLivermore_EMPhysList.hh"
#include "bmPhysList495.hh"
#include "bmPrimaryGeneratorAction.hh"
#include "bmRunAction.hh"
#include "bmEventAction.hh"
#include "bmSteppingAction.hh"
#include "bmSteppingVerbose.hh"
#include "bmAnalysisManager.hh"
#include "G4UnitsTable.hh"


#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIExecutive.hh"
#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc, char** argv) {
	
	if(argc < 2) {
		G4cout << "Usage:" << G4endl << "\t" << argv[0] << " <macro filename> [physics list]" << G4endl;
		return 0;
	}
	std::string physlist = (argc >= 3)?argv[2]:"livermore";
	
	// User Verbose stepping output class
	G4VSteppingVerbose::SetInstance(new bmSteppingVerbose());
	
	// Run manager
	G4RunManager* runManager = new G4RunManager;
	
	// User Initialization classes (mandatory)
	bmDetectorConstruction* detector = new bmDetectorConstruction();
	runManager->SetUserInitialization(detector);
	
	if(physlist=="livermore") {
		//runManager->SetUserInitialization(new bmLivermore_EMPhysList());
		runManager->SetUserInitialization(new bmPhysList495(false));
		//else if(physlist=="g4default")
		//	runManager->SetUserInitialization(new bmPhysicsList());
	} else if(physlist=="penelope") {
		runManager->SetUserInitialization(new bmPhysList495(true));
		//runManager->SetUserInitialization(new bmPenelope2008_EMPhysList());
	} else {
		G4cout << "***ERROR*** Unknown physics list: " << physlist << G4endl;
		exit(-1);
	}
	G4cout << "Using physics list: " << physlist << G4endl;
	
	new G4UnitDefinition("torr","torr","Pressure",atmosphere/760.);
	//G4UnitDefinition::PrintUnitsTable();
	//G4ParticleTable::GetParticleTable()->DumpTable();
	
#ifdef G4VIS_USE
	// Visualization, if you choose to have it!
	G4VisManager* visManager = new G4VisExecutive;
	visManager->Initialize();
#endif
	
	// User Action classes
	runManager->SetUserAction(new bmPrimaryGeneratorAction(detector));
	runManager->SetUserAction(new bmRunAction);
	runManager->SetUserAction(new bmEventAction);
	runManager->SetUserAction(new bmSteppingAction);  
	
	//create global analysis manager for histograms and trees
	gbmAnalysisManager = new bmAnalysisManager();
	
	// Execute input macro file
	G4UImanager * UI = G4UImanager::GetUIpointer(); 
	G4String command = "/control/execute ";
	G4String fileName = argv[1];
	UI->ApplyCommand(command+fileName);
	
	// interactive UI session
	if(argc >= 3 && std::string(argv[argc-1]) == "ui") {
		G4UIExecutive* UIuser = new G4UIExecutive(argc, argv);
		UIuser->SessionStart();
		delete UIuser;
	}
	
#ifdef G4VIS_USE
	delete visManager;
#endif
	delete runManager;
	
	return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

