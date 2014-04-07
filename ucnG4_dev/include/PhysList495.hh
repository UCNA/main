#ifndef PhysList495_h
#define PhysList495_h

#include "G4VModularPhysicsList.hh"

class G4VPhysicsConstructor;
//class StepMax;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

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
    
	//void AddStepMax();       
	//StepMax* GetStepMaxProcess() {return stepMaxProcess;};
	
private:
		
	G4double cutForGamma;
	G4double cutForElectron;
	G4double cutForPositron;
    
	G4String emName;
	G4VPhysicsConstructor* emPhysicsList;
    
	//StepMax* stepMaxProcess;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
