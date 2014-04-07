#ifndef SteppingAction_h
#define SteppingAction_h

#include <globals.hh>
#include <G4UserSteppingAction.hh>

/// user stepping action to check for and abort "trapped" events
class SteppingAction : public G4UserSteppingAction {
public:
	/// constructor
    SteppingAction();
	
	/// custom per-step action: checks computation time not exceeded
    void UserSteppingAction(const G4Step*);
	
    int GetTrappedFlag() const { return fTrappedFlag; }
    void SetTrappedFlag(int flag) { fTrappedFlag = flag; }
	double GetTimeSpent() const { return timeSpentSoFar; }
	
	/// reset trapping flags
    void Reset(){ fTrappedFlag = timeSpentSoFar = 0.; }
	
private:
    int fTrappedFlag;		//< whether current event is "trapped"
	double timeSpentSoFar;	//< CPU time spent on current event
};

#endif
