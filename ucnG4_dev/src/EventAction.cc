#include <cassert>
#include "SMExcept.hh"

#include "EventAction.hh"
#include "AnalysisManager.hh"
#include "SteppingAction.hh"

#include <G4Event.hh>
#include <G4EventManager.hh>
#include <G4TrajectoryContainer.hh>
#include <G4Trajectory.hh>
#include <G4ios.hh>

void EventAction::BeginOfEventAction(const G4Event* evt) {
	timer.Start();	
	G4cout<<"Beginning of event "<<evt->GetEventID()<<G4endl;
	smassert(!system("date"));
	((SteppingAction*)fpEventManager->GetUserSteppingAction())->Reset();
}

void EventAction::EndOfEventAction(const G4Event* evt) {
	timer.Stop();  
	if(evt->IsAborted())
		G4cout << "** Event aborted. **" << G4endl;
	G4cout<<"End of event "<<evt->GetEventID()<<G4endl;
	gAnalysisManager->FillTrackerData(evt);
	gAnalysisManager->FillEventTree();
}

double EventAction::getCPUTime() {
	timer.Stop();
	double t = timer.CpuTime();
	timer.Continue();
	return t;
}

