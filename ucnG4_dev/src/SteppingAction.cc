#include <TMath.h>

#include "SteppingAction.hh"
#include "AnalysisManager.hh"

#include <G4SteppingManager.hh>
#include <G4String.hh>
#include <G4EventManager.hh>
#include <G4Event.hh>

SteppingAction::SteppingAction() { 
	Reset();
}

void SteppingAction::UserSteppingAction(const G4Step* aStep) { 
	//if(aStep->GetTrack()->GetTrackStatus() == fStopAndKill && aStep->GetPostStepPoint()->GetKineticEnergy() > 0) {
	//	G4cout << "***** abnormal event abortion at " << aStep->GetPostStepPoint()->GetKineticEnergy()/keV << "keV ******" << G4endl;
	//}
	// check that computation limit is not exceeded (trapped events)
	G4int StepNo = aStep->GetTrack()->GetCurrentStepNumber();
	timeSpentSoFar = ((EventAction*)G4EventManager::GetEventManager()->GetUserEventAction())->getCPUTime();
	if(StepNo >= 200 || timeSpentSoFar > 60) {
		G4cout << "Tracking killed by computation time limit" << G4endl;
		aStep->GetTrack()->SetTrackStatus(fStopAndKill);
		SetTrappedFlag(1);  
	}
}
