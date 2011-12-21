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
// $Id: bmRunAction.cc,v 1.2 2011-10-01 18:48:38 mmendenhall Exp $
// GEANT4 tag $Name:  $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "bmRunAction.hh"
#include "bmAnalysisManager.hh"

#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithAString.hh"
#include <unistd.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

bmRunAction::bmRunAction() {
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
	fRunDir = new G4UIdirectory("/bm/run/");
	fRunDir->SetGuidance("bm customized run control");
	fRunNumberCmd = new G4UIcommand("/bm/run/runNumber",this);
	fRunNumberCmd->SetGuidance("Set the run Number.");
	fRunNumberCmd->SetParameter( new G4UIparameter("run number", 'i', true) ); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

bmRunAction::~bmRunAction() {
	delete fFileCmd;
	delete fFileDir;  
	delete fRunDir;
	delete fRunNumberCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void bmRunAction::BeginOfRunAction(const G4Run* aRun)
{
	G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
	
	//This goes into the SD manager, and get the collection ID for 
	//individual hit collections (grouped by physical detectors)
	G4cout<<"Start storing hit collection IDs ..."<<G4endl;
	gbmAnalysisManager->StoreHitCollectionIDs();
	
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void bmRunAction::EndOfRunAction(const G4Run* aRun) { 
	G4cout << "number of event = " << aRun->GetNumberOfEvent()<< G4endl;
	//write to root file
	if(gbmAnalysisManager) gbmAnalysisManager->CloseFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void bmRunAction::SetNewValue(G4UIcommand * command,G4String newValue) {
	if(command->GetCommandName() == "output") {
		// always call CloseFile() here, in case file is open
		// (user code should ignore CloseFile() if file is not open)
		if(gbmAnalysisManager) gbmAnalysisManager->CloseFile();
		// open new file, if new filename given
		// (it is okay not to give a name, in case user just wants to close file)
		if (newValue.length() <= 0) {
			G4cerr << "Null Output File Name" << G4endl;
			return;
        } else {
			gbmAnalysisManager->OpenFile(newValue);
		}
    }
	
	//runNumber UI
	else if (command->GetCommandName() == "runNumber") {
		int run_num = atoi((const char *)newValue);
		if (run_num>=0) {
			fRunNumber = run_num;
			gbmAnalysisManager->SetRunNumber(run_num);
		} else {
			G4cerr << "Run Number can not be negative!" << G4endl;
		}
	} else {
		G4cerr << "Unknown command:" << command->GetCommandName()
		<< " passed to bmEventAction::SetNewValue\n";
    }
}


