#ifndef PrimaryGeneratorMessenger_h
#define PrimaryGeneratorMessenger_h

#include "PrimaryGeneratorAction.hh"

#include <globals.hh>
#include <G4UImessenger.hh>
#include <G4UIcmdWithABool.hh>
#include <G4UIcmdWithADoubleAndUnit.hh>
#include <G4UIcmdWithAString.hh>
#include <G4UIdirectory.hh>

class PrimaryGeneratorAction;

/// messenger UI for primary event generator
class PrimaryGeneratorMessenger: public G4UImessenger {
public:
	/// constructor
    PrimaryGeneratorMessenger(PrimaryGeneratorAction*);
	/// destructor
	~PrimaryGeneratorMessenger();
    
	/// receive command
    void SetNewValue(G4UIcommand*, G4String);
    
private:
    PrimaryGeneratorAction* Action;					///< primary generator action being controlled
    G4UIdirectory*			gunDir;					///< '/gun/' commands directory
	G4UIcmdWithABool*		srcrelCmd;				///< control positions relative to source holder
	G4UIcmdWithAString*		eventFileCmd;			///< control file to read events from
	G4UIcmdWithADoubleAndUnit* sourceRadiusCmd;		///< spread vertices out by specified radius in x-y plane (for calibration source spots)
};

#endif
