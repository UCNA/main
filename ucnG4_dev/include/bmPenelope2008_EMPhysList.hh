#ifndef bmPenelope2008_EMPhysList_hh
#define bmPenelope2008_EMPhysList_hh 1

#include "G4VUserPhysicsList.hh"

/// UCNA electromagnetic physics list
class bmPenelope2008_EMPhysList: public G4VUserPhysicsList {
public:
	/// constructor
    bmPenelope2008_EMPhysList() {}
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
