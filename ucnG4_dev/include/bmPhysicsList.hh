#ifndef bmPhysicsList_h
#define bmPhysicsList_h 1

#include "globals.hh"
#include "G4VUserPhysicsList.hh"
#include "G4UImessenger.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

/// UCNA electromagnetic physics list
class bmPhysicsList: public G4VUserPhysicsList, public G4UImessenger{
public:
	/// constructor
    bmPhysicsList();
	/// set particle cuts
	virtual void SetCuts();
	
protected:
    /// list particles to consider
    virtual void ConstructParticle();
	/// list physics processes to consider
    virtual void ConstructProcess();
    /// electromagnetic physics processes
	virtual void StandardEM();
	
	/// UI commands response
	virtual void SetNewValue(G4UIcommand* command, G4String newValue);
	
	G4double defaultCutValue;	//< length cut for everything
	G4double cutForGamma;		//< length cut for photons
	G4double cutForElectron;	//< length cut for electrons
	G4double cutForPositron;	//< length cut for positrons
	
	G4UIdirectory *physDir;					//< UI directory for commands
	G4UIcmdWithADoubleAndUnit *allCutCmd;	//< UI command for setting length cuts
};

#endif
