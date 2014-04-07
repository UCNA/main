#ifndef RunAction_h
#define RunAction_h

#include <globals.hh>
#include <G4Run.hh>
#include <G4UserRunAction.hh>
#include <G4UImessenger.hh>
#include <G4UIdirectory.hh>
#include <G4UIcommand.hh>

/// user actions for whole run
class RunAction : public G4UserRunAction, public G4UImessenger {
public:
	/// constructor
    RunAction();
	/// destructor
	~RunAction();
	
	/// at begin of run
    void BeginOfRunAction(const G4Run*);
	/// at end of run (write/close files)
    void EndOfRunAction(const G4Run*);
	/// UI controls
	void SetNewValue(G4UIcommand * command,G4String newValue);

protected:
	G4UIdirectory  *fFileDir;		//< UI directory for opening output file
	G4UIcommand    *fFileCmd;		//< UI command for opening output file
	G4UIdirectory  *fRunDir;		//< UI directory for setting run number
	G4UIcommand    *fRunNumberCmd;	//< UI command for setting run number
	G4int          fRunNumber;		//< run number
};

#endif
