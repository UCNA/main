#include "bmDetectorConstructionUtils.hh"

G4Material* MaterialUser::Be = NULL;
G4Material* MaterialUser::Al= NULL;
G4Material* MaterialUser::Si= NULL;
G4Material* MaterialUser::Cu= NULL;
G4Material* MaterialUser::Wu= NULL;

G4Material* MaterialUser::Vacuum = NULL;
G4Material* MaterialUser::Brass = NULL;
G4Material* MaterialUser::Kevlar = NULL;
G4Material* MaterialUser::Mylar = NULL;
G4Material* MaterialUser::WCPentane = NULL;
G4Material* MaterialUser::WCNitrogen = NULL;
G4Material* MaterialUser::Sci = NULL;

MaterialUser::MaterialUser() {
	static bool isConstructed = false;
	if(isConstructed) return;
	
	Vacuum = NULL; // be sure to set this later
	
	std::string name,symbol;
	int z;
	G4double a;
	new G4Element(name="H",		symbol="H", z=1,  a=1.0079*g/mole);
	new G4Element(name="C",		symbol="C", z=6,  a=12.0107*g/mole);
	new G4Element(name="N",		symbol="N", z=7,  a=14.0067*g/mole);
	new G4Element(name="O",		symbol="O", z=8,  a=15.9994*g/mole);
	new G4Element(name="Al",	symbol="Al",z=13, a=26.9815*g/mole);
	new G4Element(name="Cu",	symbol="Cu",z=29, a=63.55*g/mole);
	new G4Element(name="Zn",	symbol="Zn",z=30, a=65.39*g/mole);
	
	Be = new G4Material("Beryllium",4.,9.01*g/mole,1.848*g/cm3);
	Al = new G4Material("Aluminum",13.,26.98*g/mole,2.7*g/cm3);
	Si = new G4Material("Silicon",14.,28.09*g/mole,2.33*g/cm3);
	Cu = new G4Material("Copper", 29., 63.55*g/mole, 8.96*g/cm3);
	Wu = new G4Material("Tungsten",74.,183.84*g/mole,19.3*g/cm3);
	
	Brass = new G4Material("Brass", 8.5*g/cm3,2);
	Brass->AddElement(G4Element::GetElement("Cu"),7);
	Brass->AddElement(G4Element::GetElement("Zn"),3);
	
	Kevlar=new G4Material("Kevlar",1.44*g/cm3,4);
	Kevlar->AddElement(G4Element::GetElement("N"),2);
	Kevlar->AddElement(G4Element::GetElement("C"),14);
	Kevlar->AddElement(G4Element::GetElement("H"),10);
	Kevlar->AddElement(G4Element::GetElement("O"),2);
	
	Mylar=new G4Material("Mylar",1.4*g/cm3,3);
	Mylar->AddElement(G4Element::GetElement("C"),5);
	Mylar->AddElement(G4Element::GetElement("H"),4);
	Mylar->AddElement(G4Element::GetElement("O"),2);
	
	//Wirechamber fill: pentane @ 100torr
	WCPentane=new G4Material("Pentane",0.389*mg/cm3,2);
	WCPentane->AddElement(G4Element::GetElement("C"),5);
	WCPentane->AddElement(G4Element::GetElement("H"),12);
	//Wirechamber fill: Nitrogen @ 100torr
	WCNitrogen = new G4Material("Nitrogen",7.,14.01*g/mole,1.251*100/760*mg/cm3);
	
	//Scintillator
	Sci=new G4Material("Scintillator",1.032*g/cm3,2);
	Sci->AddElement(G4Element::GetElement("C"),9);
	Sci->AddElement(G4Element::GetElement("H"),10);
	
	isConstructed = true;
}

void MaterialUser::setVacuumPressure(G4double pressure) {
	// our slightly crappy vacuum: low-pressure air (density @20c; 1.290*mg/cm3 @STP)
	G4cout<<"+++++++++++++++++ Detector vacuum is set at "<<pressure/(atmosphere/760)<<" Torr"<<G4endl;
	Vacuum = new G4Material("Vacuum",1.2048*mg/cm3*pressure/atmosphere,2);
	Vacuum->AddElement(G4Element::GetElement("N"),0.7);
	Vacuum->AddElement(G4Element::GetElement("O"),0.3);	
}
