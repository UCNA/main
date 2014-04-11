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

int main(int argc, char** argv) {
	
	// User Verbose stepping output class
	G4VSteppingVerbose::SetInstance(new SteppingVerbose());

	// Run manager
	G4RunManager* runManager = new G4RunManager;
	
	// User Initialization classes
	DetectorConstruction* detector = new DetectorConstruction();
	runManager->SetUserInitialization(detector);
	runManager->SetUserInitialization(new PhysList495());
	// User Action classes
	runManager->SetUserAction(new PrimaryGeneratorAction(detector));
	runManager->SetUserAction(new RunAction);
	runManager->SetUserAction(new EventAction);
	runManager->SetUserAction(new SteppingAction);
	
	new G4UnitDefinition("torr","torr","Pressure",atmosphere/760.);
	
#ifdef G4VIS_USE
	// Visualization, if you choose to have it!
	G4VisManager* visManager = new G4VisExecutive;
	visManager->Initialize();
#endif
		
	//create global analysis manager for histograms and trees
	gAnalysisManager = new AnalysisManager();
	
	// Execute input macro file, or enter interactive mode
	if(argc >= 2) {
		G4UImanager * UI = G4UImanager::GetUIpointer();
		G4String command = "/control/execute ";
		G4String fileName = argv[1];
		UI->ApplyCommand(command+fileName);
	} else {
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
