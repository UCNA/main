#ifndef Livermore_EMPhysList_hh
#define Livermore_EMPhysList_hh

#include "G4VUserPhysicsList.hh"

/// UCNA electromagnetic physics list based on Livermore low-energy routines
class Livermore_EMPhysList: public G4VUserPhysicsList {
public:
	/// constructor
    Livermore_EMPhysList() {}
	/// set particle cuts
	virtual void SetCuts();
	
protected:
    /// list particles to consider
    virtual void ConstructParticle();
	/// list physics processes to consider
    virtual void ConstructProcess();
    /// electromagnetic physics processes
	virtual void StandardEM();
};

#endif
