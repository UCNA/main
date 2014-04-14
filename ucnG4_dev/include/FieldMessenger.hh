#ifndef FIELDMESSENGER_HH
#define FIELDMESSENGER_HH

#include "Field.hh"

#include <G4UImessenger.hh>
#include <G4UIdirectory.hh>
#include <G4UIcmdWithAString.hh>
#include <G4UIcmdWithADouble.hh>

/// UI for controlling magnetic field settings
class FieldMessenger: public G4UImessenger {
public:
	/// constructor
	FieldMessenger(Field*);
	/// destructor
	virtual ~FieldMessenger();
	
	/// respond to UI input
	virtual void SetNewValue(G4UIcommand*, G4String);

private:

	Field* myField;								///< field being controlled
	
	G4UIdirectory* fFieldDir;					///< UI directory for field commands
	
	G4UIcmdWithADouble* fFieldScaleCmd;			///< command for setting magnetic field scale factor
	G4UIcmdWithADouble* fAfpDipoleCmd;			///< command for setting AFP dipole strength
	G4UIcmdWithAString* fFieldMapFileCmd;		///< which field map to use
};

#endif
