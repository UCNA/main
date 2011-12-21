#ifndef bmLivermore_EMPhysList_hh
#define bmLivermore_EMPhysList_hh 1

#include "G4VUserPhysicsList.hh"

/// UCNA electromagnetic physics list based on Livermore low-energy routines
class bmLivermore_EMPhysList: public G4VUserPhysicsList {
public:
	/// constructor
    bmLivermore_EMPhysList() {}
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
