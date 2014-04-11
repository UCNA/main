#ifndef PhysList495_h
#define PhysList495_h

#include <G4VModularPhysicsList.hh>
#include <G4VPhysicsConstructor.hh>

class PhysList495: public G4VModularPhysicsList {
public:
	
	PhysList495(bool usePenelope);
	virtual ~PhysList495();
	
	void ConstructParticle();
    
	void SetCuts();
	void SetCutForGamma(G4double);
	void SetCutForElectron(G4double);
	void SetCutForPositron(G4double);
	
	void ConstructProcess();
	
private:
		
	G4double cutForGamma;
	G4double cutForElectron;
	G4double cutForPositron;
    
	G4String emName;
	G4VPhysicsConstructor* emPhysicsList;
};

#endif
