#include "AnalysisManager.hh"
#include "PhysList495.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"
#include "SteppingVerbose.hh"
#include "DetectorConstruction.hh"

#include <G4UnitsTable.hh>
#include <G4RunManager.hh>
#include <G4UImanager.hh>
#include <G4UIExecutive.hh>
#ifdef G4VIS_USE
#include <G4VisExecutive.hh>
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc, char** argv) {
	
	if(argc < 2) {
		G4cout << "Usage:" << G4endl << "\t" << argv[0] << " <macro filename> [physics list]" << G4endl;
		return 0;
	}
	std::string physlist = (argc >= 3)?argv[2]:"livermore";
	
	// User Verbose stepping output class
	G4VSteppingVerbose::SetInstance(new SteppingVerbose());
	
	// Run manager
	G4RunManager* runManager = new G4RunManager;
	
	// User Initialization classes (mandatory)
	DetectorConstruction* detector = new DetectorConstruction();
	runManager->SetUserInitialization(detector);
	
	if(physlist=="livermore") {
		runManager->SetUserInitialization(new PhysList495(false));
	} else if(physlist=="penelope") {
		runManager->SetUserInitialization(new PhysList495(true));
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
	runManager->SetUserAction(new PrimaryGeneratorAction(detector));
	runManager->SetUserAction(new RunAction);
	runManager->SetUserAction(new EventAction);
	runManager->SetUserAction(new SteppingAction);  
	
	//create global analysis manager for histograms and trees
	gAnalysisManager = new AnalysisManager();
	
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

