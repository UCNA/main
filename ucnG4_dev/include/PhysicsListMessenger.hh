#ifndef PHYSICSLISTMESSENGER_HH
#define PHYSICSLISTMESSENGER_HH

#include "PhysList495.hh"

#include <G4UImessenger.hh>
#include <G4UIdirectory.hh>
#include <G4UIcmdWithAString.hh>

/// UI for selecting physics list properties
class PhysicsListMessenger: public G4UImessenger {
public:
	/// constructor
	PhysicsListMessenger(PhysList495*);
	/// destructor
	~PhysicsListMessenger();
	/// respond to UI input
	virtual void SetNewValue(G4UIcommand*, G4String);

private:
	PhysList495* fPhysList;			///< physics list being managed
	
	G4UIdirectory* fPhysDir;		///< UI directory for physics list commands
	G4UIcmdWithAString* fListCmd;	///< command for selecting physics list
};

#endif
