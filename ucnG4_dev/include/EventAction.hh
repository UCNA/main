#ifndef EVENTACTION_HH
#define EVENTACTION_HH

#include <TStopwatch.h>

#include <G4UserEventAction.hh>
#include <G4Event.hh>

/// user actions for each event (mainly gets computation time)
class EventAction : public G4UserEventAction {
public:
	/// constructor
    EventAction(): timer() {}
	/// perform at start of event simulation
    void BeginOfEventAction(const G4Event*);
	/// perform at end of event simulation
    void EndOfEventAction(const G4Event*);
	/// get computation time spent so far
	double getCPUTime();
protected:
	TStopwatch timer;	//< event computation time timer
};

#endif
