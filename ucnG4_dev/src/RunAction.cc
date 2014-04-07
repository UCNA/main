#include "RunAction.hh"
#include "AnalysisManager.hh"

#include <unistd.h>

#include <G4Run.hh>
#include <G4UImanager.hh>
#include <G4UIcommand.hh>
#include <G4UIcmdWithAString.hh>

RunAction::RunAction() {
	//too lazy to create a messenger class; include inside this class
	fFileDir = new G4UIdirectory("/files/");
	fFileDir->SetGuidance("File IO control");
	fFileCmd= new G4UIcommand("/files/output",this);
	fFileCmd->SetGuidance("Set the file name for event output.");
	fFileCmd->SetParameter( new G4UIparameter("filename", 's', true) ); 
	
	//GEANT4 has an internal variable runID. Suprisingly it is not a "UI"
	//variable, i.e. it can not be set in G4RunMessenger.
	//GEANT4 also suggests not to change G4Run object in Begin/EndRunActions.
	//So create a run number UI command, and a variable in
	//dywAnalysisManager::fRunNumber.
	//Jianglai 10-01-2006
	fRunDir = new G4UIdirectory("//run/");
	fRunDir->SetGuidance(" customized run control");
	fRunNumberCmd = new G4UIcommand("//run/runNumber",this);
	fRunNumberCmd->SetGuidance("Set the run Number.");
	fRunNumberCmd->SetParameter( new G4UIparameter("run number", 'i', true) ); 
}

RunAction::~RunAction() {
	delete fFileCmd;
	delete fFileDir;  
	delete fRunDir;
	delete fRunNumberCmd;
}

void RunAction::BeginOfRunAction(const G4Run* aRun) {
	G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
	
	//This goes into the SD manager, and get the collection ID for 
	//individual hit collections (grouped by physical detectors)
	G4cout<<"Start storing hit collection IDs ..."<<G4endl;
	gAnalysisManager->StoreHitCollectionIDs();
	
}

void RunAction::EndOfRunAction(const G4Run* aRun) { 
	G4cout << "number of event = " << aRun->GetNumberOfEvent()<< G4endl;
	//write to root file
	if(gAnalysisManager) gAnalysisManager->CloseFile();
}

void RunAction::SetNewValue(G4UIcommand * command,G4String newValue) {
	if(command->GetCommandName() == "output") {
		// always call CloseFile() here, in case file is open
		// (user code should ignore CloseFile() if file is not open)
		if(gAnalysisManager) gAnalysisManager->CloseFile();
		// open new file, if new filename given
		// (it is okay not to give a name, in case user just wants to close file)
		if (newValue.length() <= 0) {
			G4cerr << "Null Output File Name" << G4endl;
			return;
        } else {
			gAnalysisManager->OpenFile(newValue);
		}
    }
	
	//runNumber UI
	else if (command->GetCommandName() == "runNumber") {
		int run_num = atoi((const char *)newValue);
		if (run_num>=0) {
			fRunNumber = run_num;
			gAnalysisManager->SetRunNumber(run_num);
		} else {
			G4cerr << "Run Number can not be negative!" << G4endl;
		}
	} else {
		G4cerr << "Unknown command:" << command->GetCommandName()
		<< " passed to EventAction::SetNewValue\n";
    }
}
