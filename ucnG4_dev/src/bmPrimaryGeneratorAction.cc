//// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: bmPrimaryGeneratorAction.cc,v 1.34 2011-12-12 01:29:37 mmendenhall Exp $
// GEANT4 tag $Name:  $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "bmPrimaryGeneratorAction.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "BetaSpectrum.hh"
#include "globals.hh"
#include "Randomize.hh"
#include "bmAnalysisManager.hh"
#include <iostream>
#include <fstream>
#include <vector>
#include <TString.h>
#include <TRandom.h>
#include <TF2.h>

#include <unistd.h>
#include <vector>
#include <algorithm>
#include <cassert>
#include <bitset>  // djaffe 2aug06

namespace CLHEP {}
using namespace CLHEP; 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//for beta spectrum
TF1* funcBetaSpectrum = NULL;
TF1* funcHeavyBeta = NULL;
TF2* funcNeutronBeta = NULL;

double rand_outof_list(const double* line, const double* weight, const int npoints, int& selected) {
	// build cumulative list
	std::vector<double> division(npoints+1,0.0);
	for(int i=0; i<npoints; i++)
		division[i+1] = division[i]+weight[i];
	// pick selected line
	std::vector<double>::iterator itsel = std::upper_bound(division.begin(),division.end(),division.back()*G4UniformRand());
	selected = int(itsel-division.begin())-1;	
	assert(0<=selected && selected<npoints);
	return line[selected];
}

void RandomizeMomentum(G4ThreeVector& mom) {
	// isotropic emission
	G4double phi = 2.0*M_PI*G4UniformRand();
	G4double costheta = 2.0*G4UniformRand()-1.0;
	G4double sintheta = sqrt(1.0-costheta*costheta);
	mom.set(cos(phi)*sintheta,sin(phi)*sintheta,costheta);
}  


/// generate a random position in a tube
void randomTubePosition(const G4ThreeVector centerpos, const G4double radius, const G4double halfz, G4ThreeVector& pos) {
	G4double x0,y0;
	while(true) {
		x0 = (2.0*G4UniformRand()-1.)*radius;
		y0 = (2.0*G4UniformRand()-1.)*radius;
		if(x0*x0+y0*y0<=radius*radius) break;
	}
	G4double z0=(2.0*G4UniformRand()-1.0)*halfz;
	pos = centerpos+G4ThreeVector(x0,y0,z0);
}

double genericBetaSpectrum(double* x, double *par) {
	double KE = x[0];	// particle kinetic energy
	double Q = par[0];	// spectrum endpoint
	if(KE>=0 && KE<=Q) return (Q-KE)*(Q-KE)*(KE+m_e)*sqrt(KE*(KE+2*m_e));
	else return 0;
}

/// beta decay spectrum from heavy nucleus
double heavyBetaSpectrum(double* x, double* par) {
	double KE = x[0];	// beta kinetic energy
	double Q = par[0];	// spectrum endpoint
	double A = par[1];	// nucleus A
	double Z = par[2];	// nucleus Z
	double W = (KE+m_e)/m_e;
	double W0 = (Q+m_e)/m_e;
	double R = pow(A,1./3.)*neutron_R0;
	
	//G4cout << "HeavyBeta: KE=" << KE << " Q=" << Q << " Z=" << Z << " W=" << W << " W0=" << W0 << " R=" << R << " s=" << plainPhaseSpace(W,W0) << " F0=" << WilkinsonF0(Z,W,R) << G4endl;
	
	if(0<KE && KE<Q)
		return plainPhaseSpace(W,W0)*WilkinsonF0(Z,W,R)*WilkinsonL0(Z,W,R)*(1.+Wilkinson_g(W,W0));
	return 0;
}

double asymNeutronBetaSpectrum(double* x, double*) {
	double e0 = x[0];		// kinetic energy
	double costh = x[1];	// cos theta
	
	double nbe = neutronBetaEp;
	double beta = sqrt(e0*e0+2*m_e*e0)/(m_e+e0);
	return genericBetaSpectrum(x,&nbe)*(1+A0_PDG*beta*costh);
}

void bmPrimaryGeneratorAction::polarizedNeutronBetaDecayGenerator(G4Event* anEvent, bool flipper) {
	
	// choose kinetic energy and cos theta
	double costheta, beta_energy;
	funcNeutronBeta->GetRandom2(beta_energy,costheta);
	if(flipper)
		costheta *= -1.0;
	
	// choose transverse momentum components
	double sintheta=sqrt(1.0-costheta*costheta);
	double phi=6.28318531*G4UniformRand();
	G4ThreeVector direction(cos(phi)*sintheta,sin(phi)*sintheta,costheta);
	
	// throw event
	G4cout << "Throwing neutron beta decay with KE = " << beta_energy << "keV, cos(theta)=" << costheta << G4endl;
	particleGun->SetParticleEnergy(beta_energy*keV);
	particleGun->SetParticleMomentumDirection(direction);
	particleGun->GeneratePrimaryVertex(anEvent);
}


void bmPrimaryGeneratorAction::
throwElectronsAndGammas(const std::vector<G4double>& electrons,
						const std::vector<G4double>& gammas,
						G4Event* anEvent) {
	
	
	G4cout << "Throwing " << electrons.size() << " electrons and " << gammas.size() << " gammas." << G4endl;
	
	G4ThreeVector direction;
	
	// generate electrons
	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	particleGun->SetParticleDefinition(particleTable->FindParticle("e-"));
	for(std::vector<G4double>::const_iterator it = electrons.begin(); it != electrons.end(); it++) {
		G4cout << "\te- = " << *it/keV << "keV" << G4endl;
		RandomizeMomentum(direction);
		particleGun->SetParticleEnergy(*it);
		particleGun->SetParticleMomentumDirection(direction);
		particleGun->GeneratePrimaryVertex(anEvent);
	}
	
	// generate photons
	particleGun->SetParticleDefinition(particleTable->FindParticle("gamma"));
	for(std::vector<G4double>::const_iterator it = gammas.begin(); it != gammas.end(); it++) {	
		G4cout << "\tgamma = " << *it/keV << "keV" << G4endl;
		RandomizeMomentum(direction);
		particleGun->SetParticleEnergy(*it);
		particleGun->SetParticleMomentumDirection(direction);
		particleGun->GeneratePrimaryVertex(anEvent);
	}
}

//position is in units of mm, mom is dimensionless, energy is in units of keV
void bmPrimaryGeneratorAction::
Sn113SourceGenerator(G4Event* anEvent) {
	
	G4ThreeVector direction;
	G4double gunEnergy;
	RandomizeMomentum(direction);
	
	const double line[4] = {363.758, 387.461, 390.872, 391.576};
	const double branching[4] = {28.4, 5.57, 1.108, 0.2449};
	
	int selected(-1);
	gunEnergy = rand_outof_list(line, branching, 4, selected)*keV;
	
	particleGun->SetParticleEnergy(gunEnergy);
	particleGun->SetParticleMomentumDirection(direction);
	particleGun->GeneratePrimaryVertex(anEvent);
	
	//now add in correlated Auger electrons
	if(selected==0) {
		double ran = G4UniformRand();
		if(ran<0.1405) { //Auger coefficient
			gunEnergy = 20.1*keV;
			RandomizeMomentum(direction);
			particleGun->SetParticleEnergy(gunEnergy);
			particleGun->SetParticleMomentumDirection(direction);
			particleGun->GeneratePrimaryVertex(anEvent);
		}
	}
	return;
}

// Sn113 with gammas, data based on NuDat 2.0 www.nndc.bnl.gov/nudat2/decaysearchdirect.jsp?nuc=113SN&unc=nds
//position is in units of mm, mom is dimensionless, energy is in units of keV
void bmPrimaryGeneratorAction::
Sn113GSourceGenerator(G4Event* anEvent) {
	G4double gunEnergy;    	
	// choose gamma decay path
	std::vector<G4double> electrons;
	std::vector<G4double> gammas;
	int selected(-1);
	unsigned int nklines = 0;
	
	const double p646 = 2.11 + 0.082+0.114+0.022+0.0004;
	if(100.0*G4UniformRand() < p646) {
		/// 255keV gamma and conversions
		const double line255[] =	{255.134,	227.194,	250.896,	254.308,	255.012};
		const double branch255[] =	{2.11,		0.082,		0.0114,		0.0022,		0.0004};
		gunEnergy = rand_outof_list(line255, branch255, 5, selected)*keV;
		if(!selected)
			gammas.push_back(gunEnergy);
		else
			electrons.push_back(gunEnergy);
		if(selected == 1)
			nklines++;
	}
	
	const double line392[] =	{391.698,	363.758,	387.461,	390.872,	391.576,	391.697};
	const double branch392[] =	{64.97,			28.8,		5.60,		1.137,		0.205,		0.01260};
	gunEnergy = rand_outof_list(line392, branch392, 6, selected)*keV;
	if(!selected)
		gammas.push_back(gunEnergy);
	else
		electrons.push_back(gunEnergy);
	if(selected == 1)
		nklines++;
	
	// uncorrelated K auger, without other K lines
	if(!nklines) {
		if(G4UniformRand() < 0.170-0.1405*(0.00082+0.288))
			electrons.push_back(20.1*keV);
	}
	// correlated K auger
	while(nklines-- > 0) {
		if(G4UniformRand() < 0.1405)
			electrons.push_back(20.1*keV);
	}
	
	throwElectronsAndGammas(electrons,gammas,anEvent);
}


//position is in units of mm, mom is dimensionless, energy is in units of keV
void bmPrimaryGeneratorAction::
Bi207SourceGenerator(G4Event* anEvent) {
	
	G4ThreeVector direction;
	G4double gunEnergy;
	
	RandomizeMomentum(direction);
	
	//note, the total Auger K intensity is 2.9% per decay. 
	//The uncorrelated Auger K is computed as 
	//2.9%-(1.515%+7.03%)*3.66% where 3.66% is the Auger coefficent estimated for
	//Pb207
	
	const double line[9] = {56.7, 481.6935, 553.8372, 565.8473, 809.80, 881.94,
		975.651, 1047.795, 1059.805};
	const double branching[9] = {2.587, 1.515, 0.438, 0.147, 0.00263, 4.45E-4,
		7.03, 1.84, 0.54};
	
	int selected(-1);
	gunEnergy = rand_outof_list(line, branching, 9, selected)*keV;
	
	particleGun->SetParticleEnergy(gunEnergy);
	particleGun->SetParticleMomentumDirection(direction);
	particleGun->GeneratePrimaryVertex(anEvent);
	
	//now add in correlated Auger electrons 
	if(selected==1||selected==4||selected==6) {
		double ran = G4UniformRand();
		if(ran<0.0366) { //Auger electron, 3.66% Auger coefficent
			gunEnergy = 56.7*keV;                            
			RandomizeMomentum(direction);
			particleGun->SetParticleEnergy(gunEnergy);
			particleGun->SetParticleMomentumDirection(direction);
			particleGun->GeneratePrimaryVertex(anEvent);
		}
	}
	
	return;
}

//more complicated Bi207 including gammas, correlated decays
//data from NuDat 2.5, www.nndc.bnl.gov/nudat2
//position is in units of mm, mom is dimensionless, energy is in units of keV
void bmPrimaryGeneratorAction::Bi207GSourceGenerator(G4Event* anEvent) {
	
	G4double gunEnergy;
	
	// 2350 energy level branching
	const double p1442 = 0.1310+0.00035;
	const double p1770 = 6.87+0.0241+0.0034;
	
	// select which level we initially land on:
	//	0: 560keV
	//	1: 1630keV
	//	2: 2350keV
	const double p1630 = 74.6 + 0.54+1.84+7.03;			// observed gammas+electrons percent
	const double p2350 = 6.87+0.1310 + 0.0034+0.0241;	// observed gammas+electrons percent
	unsigned int gammaChoice = 0;
	double randsel = 100.0*G4UniformRand();
	if(randsel < p1630+p2350)
		gammaChoice = 2;
	if(randsel < p1630)
		gammaChoice = 1;
	
	std::vector<G4double> electrons;
	std::vector<G4double> gammas;
	int selected(-1);
	int nKlines = 0;	// number of K decays with possible correlated auger
	
	if(gammaChoice == 0 || gammaChoice == 1 || (gammaChoice == 2 && G4UniformRand()<p1770/(p1770+p1442))) {
		// 569keV conversions
		const double line569[] =		{569.698,	481.6935,	553.8372,	565.8473};
		const double branching569[] =	{97.76,		1.515,		0.438,		0.147};
		gunEnergy = rand_outof_list(line569, branching569, 4, selected)*keV;
		if(selected)
			electrons.push_back(gunEnergy);
		else
			gammas.push_back(gunEnergy);
		if(selected == 1 || selected == 3)
			nKlines++;
		
		if(gammaChoice == 2) {
			// 1770keV conversions
			const double line1770[] =		{1770.228,	1682.224,	1754.367};
			const double branching1770[] =	{6.87,		0.0241,		0.0034};
			gunEnergy = rand_outof_list(line1770, branching1770, 3, selected)*keV;
			if(selected)
				electrons.push_back(gunEnergy);
			else
				gammas.push_back(gunEnergy);
			if(selected == 1)
				nKlines++;
		}
	} else if(gammaChoice == 2) {
		// 1442keV conversions
		const double line1442[] =		{1442.2,	1354.20,	1426.34};
		const double branching1442[] =	{0.1310,	3.54E-4,	5.50E-5};
		gunEnergy = rand_outof_list(line1442, branching1442, 3, selected)*keV;
		if(selected)
			electrons.push_back(gunEnergy);
		else
			gammas.push_back(gunEnergy);
		if(selected == 1)
			nKlines++;
		// 897keV conversions
		const double line897[] =		{897.8,		 881.94,	809.80};
		const double branching897[] =	{0.131,		4.45E-4,	0.00263};
		gunEnergy = rand_outof_list(line897, branching897, 3, selected)*keV;
		if(selected)
			electrons.push_back(gunEnergy);
		else
			gammas.push_back(gunEnergy);	
		if(selected == 1)
			nKlines++;
	}
	
	// 1063keV conversions
	if(gammaChoice == 1) {
		const double line1063[] =		{1063.656,	975.651,	1047.795,	1059.805};
		const double branching1063[] =	{74.6,		7.03,		1.84,		0.54};
		gunEnergy = rand_outof_list(line1063, branching1063, 4, selected)*keV;
		if(selected)
			electrons.push_back(gunEnergy);
		else
			gammas.push_back(gunEnergy);
		if(selected == 1)
			nKlines++;
	}
	
	// uncorrelated K auger, without other K lines
	if(!nKlines) {
		if(G4UniformRand() < 0.029-0.0366*(0.0703+0.01515))
			electrons.push_back(56.7*keV);
	}
	
	// correlated K auger, 3.66% chance with other K lines
	while(nKlines-- > 0) {
		if(G4UniformRand() < 0.0366)
			electrons.push_back(56.7*keV);
	}
	
	throwElectronsAndGammas(electrons,gammas,anEvent);
}


//position is in units of mm, mom is dimensionless, energy is in units of keV
void bmPrimaryGeneratorAction::
Ce139SourceGenerator(G4Event* anEvent) {
	
	G4ThreeVector direction;
	G4double gunEnergy;
	RandomizeMomentum(direction);
	
	const double line[4] = {126.9329, 159.5912, 164.4962, 165.5871};
	const double branching[4] = {17.1, 2.30, 0.48, 0.131};
	int selected(-1);
	gunEnergy = rand_outof_list(line, branching, 4, selected)*keV;
	
	particleGun->SetParticleEnergy(gunEnergy);
	particleGun->SetParticleMomentumDirection(direction);
	particleGun->GeneratePrimaryVertex(anEvent);
	
	//now add in correlated Auger electrons 
	if(selected==0) {
		double ran = G4UniformRand();
		if(ran<0.0953) { //Auger electron, 9.53% Auger coefficent
			gunEnergy = 28.4*keV;                            
			RandomizeMomentum(direction);
			particleGun->SetParticleEnergy(gunEnergy);
			particleGun->SetParticleMomentumDirection(direction);
			particleGun->GeneratePrimaryVertex(anEvent);
		}
	}
	
	return;
}

void bmPrimaryGeneratorAction::
Ce139GSourceGenerator(G4Event* anEvent)
{
	G4double gunEnergy;
	
	std::vector<G4double> electrons;
	std::vector<G4double> gammas;
	int selected(-1);
	
	// 166keV gamma and conversions
	const double line166[] =	{166.8575,	126.9329,	159.5912,	164.4962,	165.5871};
	const double branch166[] =	{80.0,		17.146,		2.2980,		0.4751,		0.13120};
	gunEnergy = rand_outof_list(line166, branch166, 5, selected)*keV;
	if(!selected)
		gammas.push_back(gunEnergy);
	else
		electrons.push_back(gunEnergy);
	// correlated and uncorrelated auger
	if(selected == 1) {
		if(G4UniformRand() < 0.0953)
			electrons.push_back(27.4*keV);
	} else if(G4UniformRand() < 0.083-0.0953*0.17146)
		electrons.push_back(27.4*keV);
	
	throwElectronsAndGammas(electrons,gammas,anEvent);
}


void bmPrimaryGeneratorAction::
Cd109GSourceGenerator(G4Event* anEvent) {
	
	G4double gunEnergy;
	
	std::vector<G4double> electrons;
	std::vector<G4double> gammas;
	int selected(-1);
	
	// 88keV gamma and conversions
	const double line88[] =		{88.03,		62.5196,	84.2278,	87.3161,	87.9384};
	const double branch88[] =	{3.7,		41.7,		44.0,		8.9,		1.6};
	gunEnergy = rand_outof_list(line88, branch88, 5, selected)*keV;
	if(!selected)
		gammas.push_back(gunEnergy);
	else
		electrons.push_back(gunEnergy);
	
	// auger
	if(G4UniformRand() < 0.208)
		electrons.push_back(18.5*keV);
	
	throwElectronsAndGammas(electrons,gammas,anEvent);
}

// lines based on NuDat 2.6
void bmPrimaryGeneratorAction::
Cd113mGSourceGenerator(G4Event* anEvent) {
	
	std::vector<G4double> electrons;
	std::vector<G4double> gammas;
	
	funcHeavyBeta->FixParameter(0,585.7);
	funcHeavyBeta->FixParameter(1,113.);
	funcHeavyBeta->FixParameter(2,49.);
	electrons.push_back(funcHeavyBeta->GetRandom()*keV);
	
	throwElectronsAndGammas(electrons,gammas,anEvent);
}

// lines based on NuDat 2.6
void bmPrimaryGeneratorAction::
Cs137GSourceGenerator(G4Event* anEvent) {
	
	G4double gunEnergy;
	
	std::vector<G4double> electrons;
	std::vector<G4double> gammas;
	int selected(-1);
	
	// Beta spectrum
	const double lineBeta[] =	{513.97,	1175.63};
	const double branchBeta[] =	{94.70,		5.30};
	gunEnergy = rand_outof_list(lineBeta, branchBeta, 2, selected);
	funcHeavyBeta->FixParameter(0,gunEnergy);
	funcHeavyBeta->FixParameter(1,137.);
	funcHeavyBeta->FixParameter(2,56.);
	electrons.push_back(funcHeavyBeta->GetRandom()*keV);
	
	// lines
	if(G4UniformRand()*100. < 85.10+7.79+1.402+0.300+0.0646) {
		// 662keV line:				gamma		CE K		CE L		CE M		CE N
		const double line662[] =	{661.657,	624.216,	665.668,	660.364,	661.404};
		const double branch662[] =	{85.10,		7.79,		1.402,		0.300,		0.0646};
		gunEnergy = rand_outof_list(line662, branch662, 5, selected)*keV;
		if(selected)
			electrons.push_back(gunEnergy);
		else
			gammas.push_back(gunEnergy);
	}
	
	throwElectronsAndGammas(electrons,gammas,anEvent);
}

// lines based on NuDat 2.6
void bmPrimaryGeneratorAction::In114GSourceGenerator(G4Event* anEvent) {
	
	G4double gunEnergy;
	
	std::vector<G4double> electrons;
	std::vector<G4double> gammas;
	int selected(-1);
	
	if(G4UniformRand() < 0.5) {
		// 5+ decay to:				114Cd ------------	or 1+ state ------------------
		const double line190[] =	{725.24,	558.43,	190.27,	189.44,	186.03,	162.33};
		const double branch190[] =	{3.2,		3.2,	15.56,	6.71,	31.9,	40.1};
		gunEnergy = rand_outof_list(line190, branch190, 4, selected)*keV;
		if(selected < 3)
			gammas.push_back(gunEnergy);
		else
			electrons.push_back(gunEnergy);
		// auger K
		if(selected==5 && G4UniformRand() < 0.0598/0.401)
			electrons.push_back(20.1*keV);
	} else {
		// 1+ state 1988.7keV beta decay to 114Sn
		if(G4UniformRand() < 0.9936*0.94) {
			funcHeavyBeta->FixParameter(0,1988.7);
			funcHeavyBeta->FixParameter(1,114.);
			funcHeavyBeta->FixParameter(2,50.);
			electrons.push_back(funcHeavyBeta->GetRandom()*keV);
		}
	}
	
	throwElectronsAndGammas(electrons,gammas,anEvent);
}

// lines based on NuDat 2.6
void bmPrimaryGeneratorAction::Sc46GSourceGenerator(G4Event* anEvent) {
	std::vector<G4double> electrons;
	std::vector<G4double> gammas;
	funcHeavyBeta->FixParameter(0,356.9);
	funcHeavyBeta->FixParameter(1,46.);
	funcHeavyBeta->FixParameter(2,22.);
	electrons.push_back(funcHeavyBeta->GetRandom()*keV);
	gammas.push_back(889.277*keV);
	gammas.push_back(1120.545*keV);
	throwElectronsAndGammas(electrons,gammas,anEvent);
}

// lines based on NuDat 2.6
void bmPrimaryGeneratorAction::Co60GSourceGenerator(G4Event* anEvent) {
	std::vector<G4double> electrons;
	std::vector<G4double> gammas;
	funcHeavyBeta->FixParameter(0,318.2);
	funcHeavyBeta->FixParameter(1,60.);
	funcHeavyBeta->FixParameter(2,28.);
	electrons.push_back(funcHeavyBeta->GetRandom()*keV);
	gammas.push_back(1173.228*keV);
	gammas.push_back(1332.492*keV);
	throwElectronsAndGammas(electrons,gammas,anEvent);
}

// lines based on NuDat 2.6
void bmPrimaryGeneratorAction::Ag110GSourceGenerator(G4Event* anEvent) {
	std::vector<G4double> electrons;
	std::vector<G4double> gammas;
	
	// 6+ beta decays
	int selected(-1);
	const double lineBetas[] =	{83.0,	529.9};
	const double branchBetas[] ={67.0,	30.2};
	G4double gunEnergy = rand_outof_list(lineBetas, branchBetas, 2, selected)*keV;
	
	funcHeavyBeta->FixParameter(0,gunEnergy);
	funcHeavyBeta->FixParameter(1,110.);
	funcHeavyBeta->FixParameter(2,48.);
	electrons.push_back(funcHeavyBeta->GetRandom()*keV);
	
	// 6+ gammas
	if(G4UniformRand() < 0.943)
		gammas.push_back(657.76*keV);
	if(G4UniformRand() < 0.1056)
		gammas.push_back(677.62*keV);
	if(G4UniformRand() < 0.1633)
		gammas.push_back(706.68*keV);
	if(G4UniformRand() < .2262)
		gammas.push_back(763.94*keV);
	if(G4UniformRand() < .727)
		gammas.push_back(884.68*keV);
	if(G4UniformRand() < .342)
		gammas.push_back(937.49*keV);
	if(G4UniformRand() < .249)
		gammas.push_back(1384.29*keV);
	if(G4UniformRand() < .136)
		gammas.push_back(1505.03*keV);
	
	throwElectronsAndGammas(electrons,gammas,anEvent);
}


//position is in units of mm, mom is dimensionless, energy is in units of keV
void bmPrimaryGeneratorAction::
In114SourceGenerator(G4Event* anEvent) {
	
	G4ThreeVector direction;
	G4double gunEnergy;
	
	RandomizeMomentum(direction);
	
	const double line[4] = {162.33, 186.03, 189.44, 190.15};
	const double branching[4] = {40.1, 31.9, 6.71, 1.351};
	
	double ran = G4UniformRand();
	int selected(-1);
	if(ran<0.5) gunEnergy = rand_outof_list(line, branching, 4, selected)*keV;
	else gunEnergy = funcBetaSpectrum->GetRandom()*keV;
	
	particleGun->SetParticleEnergy(gunEnergy);
	particleGun->SetParticleMomentumDirection(direction);
	particleGun->GeneratePrimaryVertex(anEvent);
	
	return;
}

//position is in units of mm, mom is dimensionless, energy is in units of keV
void bmPrimaryGeneratorAction::
Sr85SourceGenerator(G4Event* anEvent) {
	
	G4ThreeVector direction;
	G4double gunEnergy;
	
	RandomizeMomentum(direction);
	
	const double line[2] = {498.8070, 511.9416};
	const double branching[2] = {0.6, 0.068};  
	
	double ran = G4UniformRand();
	int selected(-1);
	if(ran<0.9)
		gunEnergy = rand_outof_list(line, branching, 2, selected)*keV;
	else gunEnergy = funcBetaSpectrum->GetRandom()*keV;//simulate Sr89 background
	
	particleGun->SetParticleEnergy(gunEnergy);
	particleGun->SetParticleMomentumDirection(direction);
	particleGun->GeneratePrimaryVertex(anEvent);
	
	//now add in correlated Auger electrons 
	if(selected==0) {
		ran = G4UniformRand();
		if(ran<0.3256) { //Auger electron, 32.56% Auger coefficent
			gunEnergy = 11.4*keV;                            
			RandomizeMomentum(direction);
			particleGun->SetParticleEnergy(gunEnergy);
			particleGun->SetParticleMomentumDirection(direction);
			particleGun->GeneratePrimaryVertex(anEvent);
		}
	}
	
	return;
}

//position is in units of mm, mom is dimensionless, energy is in units of keV
void bmPrimaryGeneratorAction::
Cd109SourceGenerator(G4Event* anEvent) {
	
	G4ThreeVector direction;
	G4double gunEnergy;
	
	RandomizeMomentum(direction);
	
	const double line[4] = {62.5196, 84.2278, 87.3161, 87.9384};
	const double branching[4] = {41.7, 44.0, 8.9, 1.60};
	
	double ran = G4UniformRand();
	int selected(-1);
	if(ran<0.85) gunEnergy = rand_outof_list(line, branching, 4, selected)*keV;
	//based on source spec sheet, 15% of Cd113m beta decay
	else gunEnergy = funcBetaSpectrum->GetRandom()*keV;
	
	particleGun->SetParticleEnergy(gunEnergy);
	particleGun->SetParticleMomentumDirection(direction);
	particleGun->GeneratePrimaryVertex(anEvent);
	
	return;
}

//position is in units of mm, mom is dimensionless, energy is in units of keV
void bmPrimaryGeneratorAction::
Mn54SourceGenerator(G4Event* anEvent) {
	
	G4ThreeVector direction;
	G4double gunEnergy; 
	
	RandomizeMomentum(direction);
	
	const double line[1] = {828.859};
	const double branching[1] = {0.0224};
	
	int selected(-1);
	gunEnergy = rand_outof_list(line, branching, 1, selected)*keV;
	
	particleGun->SetParticleEnergy(gunEnergy);
	particleGun->SetParticleMomentumDirection(direction);
	particleGun->GeneratePrimaryVertex(anEvent);
	
	return;
}





// Xe125 1/2+ (incorrect correlations)
void bmPrimaryGeneratorAction::Xe125_1_2p_SourceGenerator(G4Event* anEvent) {
	std::vector<G4double> electrons;
	std::vector<G4double> gammas;
	int selected(-1);
	const double lines[] =		{54.968, 188.418,	243.378,	453.80,	846.51,		21.799,	49.780,	155.249,	183.230,	210.209};
	const double branches[] =	{6.80,	54.0,		30.1,		4.69,	1.11,		24.8,	3.29,	6.264,		0.886,		1.97};
	double gun_energy = rand_outof_list(lines, branches, 10, selected)*keV;
	if(selected >= 5)
		electrons.push_back(gun_energy);
	else
		gammas.push_back(gun_energy);
	throwElectronsAndGammas(electrons,gammas,anEvent);
}

// Xe133 metastable 11/2-
void bmPrimaryGeneratorAction::Xe133_11_2m_SourceGenerator(G4Event* anEvent) {
	std::vector<G4double> electrons;
	std::vector<G4double> gammas;
	int selected(-1);
	const double line233[] =	{233.221,	198.660,	227.768,	232.079,	233.013};
	const double branch233[] =	{10.,		63.5,		20.68,		4.57,		1.225};
	double gun_energy = rand_outof_list(line233, branch233, 5, selected)*keV;
	if(selected)
		electrons.push_back(gun_energy);
	else
		gammas.push_back(gun_energy);
	throwElectronsAndGammas(electrons,gammas,anEvent);
}

// Xe133 3/2+ (incorrect correlations and chains)
void bmPrimaryGeneratorAction::Xe133_3_2p_SourceGenerator(G4Event* anEvent) {
	std::vector<G4double> electrons;
	std::vector<G4double> gammas;
	int selected(-1);
	
	// beta decay
	const double betaep[] =		{346.4,	266.8};
	const double branchbeta[] =	{99.,	0.81};
	funcBetaSpectrum->FixParameter(0,rand_outof_list(betaep, branchbeta, 2, selected));
	electrons.push_back(funcBetaSpectrum->GetRandom()*keV);
	
	// gammas and electrons
	const double eline[] =	{80.997,	45.012,	75.283,	79.780,	80.766,	124.628};
	const double ebranch[] ={38.0,		55.1,	8.21,	1.69,	0.437,	0.0156};
	double gun_energy = rand_outof_list(eline, ebranch, 6, selected)*keV;
	if(selected)
		electrons.push_back(gun_energy);
	else
		gammas.push_back(gun_energy);
	
	throwElectronsAndGammas(electrons,gammas,anEvent);
}

// Xe135 metastable 11/2-
void bmPrimaryGeneratorAction::Xe135_11_2m_SourceGenerator(G4Event* anEvent) {
	std::vector<G4double> electrons;
	std::vector<G4double> gammas;
	int selected(-1);
	const double line526[] =	{526.561,	492.000,	521.108,	525.419,	526.353,	526.545};
	const double branch526[] =	{80.4,		15.34,		2.922,		0.619,		0.1275,		0.01524};
	double gun_energy = rand_outof_list(line526, branch526, 5, selected)*keV;
	if(selected)
		electrons.push_back(gun_energy);
	else
		gammas.push_back(gun_energy);
	throwElectronsAndGammas(electrons,gammas,anEvent);
}

// Xe135 3/2+ (incorrect correlations and chains)
void bmPrimaryGeneratorAction::Xe135_3_2p_SourceGenerator(G4Event* anEvent) {
	std::vector<G4double> electrons;
	std::vector<G4double> gammas;
	int selected(-1);
	
	// beta decay
	const double betaep[] =		{103.,	184.,	557.,	757.,	915.};
	const double branchbeta[] =	{0.123,	0.075,	3.11,	0.59,	96.};
	funcBetaSpectrum->FixParameter(0,rand_outof_list(betaep, branchbeta, 5, selected));
	electrons.push_back(funcBetaSpectrum->GetRandom()*keV);
	
	// gammas and electrons
	const double eline[] =	{249.794,	608.185,	213.809,	244.080,	248.577,	249.563};
	const double ebranch[] ={90.0,		2.90,		5.61,		0.82,		0.169,		0.035};
	double gun_energy = rand_outof_list(eline, ebranch, 6, selected)*keV;
	if(selected > 1)
		electrons.push_back(gun_energy);
	else
		gammas.push_back(gun_energy);
	
	throwElectronsAndGammas(electrons,gammas,anEvent);
}

// Xe137 7/2- (incorrect correlations and chains)
void bmPrimaryGeneratorAction::Xe137_7_2m_SourceGenerator(G4Event* anEvent) {
	std::vector<G4double> electrons;
	std::vector<G4double> gammas;
	int selected(-1);
	// beta decay
	const double betaep[] =		{4173.,	3718.,	3324.,	1323.,	2390.};
	const double branchbeta[] =	{67.,	31.,	0.65,	0.72,	0.38};
	funcBetaSpectrum->FixParameter(0,rand_outof_list(betaep, branchbeta, 5, selected));
	electrons.push_back(funcBetaSpectrum->GetRandom()*keV);
	if(G4UniformRand() < 0.0037)
		electrons.push_back(419.505*keV);
	throwElectronsAndGammas(electrons,gammas,anEvent);
}




bmPrimaryGeneratorAction::bmPrimaryGeneratorAction(bmDetectorConstruction* myDC): myDetector(myDC) {
	//G4int n_particle = 1;
	particleGun = new G4ParticleGun();
	
	//create a messenger for this class
	gunMessenger = new bmPrimaryGeneratorMessenger(this);
	
	// default particle
	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	G4ParticleDefinition* particle = particleTable->FindParticle("e-");
	
	
	particleGun->SetParticleDefinition(particle);
	particleGun->SetParticleTime(0.0*ns);
	
	funcBetaSpectrum = new TF1("funcBetaSpectrum",genericBetaSpectrum,0,3000,1);
	funcBetaSpectrum->SetNpx(3000);
	
	funcHeavyBeta = new TF1("funcHeavyBeta",heavyBetaSpectrum,0,5000,3);
	funcHeavyBeta->SetNpx(500);
	
	funcNeutronBeta = new TF2("funcNeutronBeta",asymNeutronBetaSpectrum,0,800,-1.0,1.0,0);
	funcNeutronBeta->SetNpx(800);
	funcNeutronBeta->SetNpy(200);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

bmPrimaryGeneratorAction::~bmPrimaryGeneratorAction() {
	delete particleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void bmPrimaryGeneratorAction::displayGunStatus() {
	G4cout << "Electron gun from " << particleGun->GetParticlePosition()/m
	<< " towards " << particleGun->GetParticleMomentumDirection()
	<< ": " << particleGun->GetParticleEnergy()/keV << "keV" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void bmPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{	
	std::bitset<32> a ( (long) 0);
	//use the UI fRunNumber in bmAnalysisManager which got set 
	//in bmRunAction, instead of the GEANT4 default RunID
	std::bitset<32> run ( gbmAnalysisManager->GetRunNumber() );
	std::bitset<32> evt ( anEvent->GetEventID() );
	
	//   if one wants to use the hostname somehow, this is the way to do it
	//   char hostname[256];
	//   gethostname(hostname,255);
	//   cout<<hostname<<endl;
	
	//this is the "host id", unique for each host
	std::bitset<32> site ( (long) gethostid() ) ; 
	
	for ( int i = 32-1-1; i >= 0 ; i-- ) {
		if ( run.test(i) ) a.set(31-1-i) ;
	}
	
	// create seed = (a XOR site) XOR evt
	std::bitset<32> seed = (a^site)^evt ;
	
	// set highest bit to zero to avoid negative seed
	if ( seed.test(31) ) seed.reset(31) ; 
	
	long myseed = seed.to_ulong();
	
	HepRandom::setTheSeed(myseed);	// random seed for Geant
	gRandom->SetSeed(myseed);		// random seed for ROOT
	
	G4cout<<"run "<<gbmAnalysisManager->GetRunNumber()<<" "
	<<"evt "<<anEvent->GetEventID()<<" "
	<<"seed "<<myseed<<G4endl;
	
	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	G4String particleName;
	G4ParticleDefinition* particle = particleTable->FindParticle(particleName="e-");
	particleGun->SetParticleDefinition(particle);
	particleGun->SetParticleTime(0.0*ns);
	
	// set vertex position for event according to position generator
	// assumed position has been set in /gun/position!
	static G4ThreeVector position_saved = particleGun->GetParticlePosition();
	G4ThreeVector vertex_position = position_saved;
	if(positioner=="Fixed") {
		// fixed position set from macro file
		vertex_position = position_saved;
	} else if(positioner=="SourceDrop") {
		// nominally 3mm disc around pre-set position for calibration sources
		randomTubePosition(position_saved, sourceRadius, 0*um, vertex_position);  
	} else if(positioner=="DecayTrapUniform") {
		// uniform fill of decay trap volume
		randomTubePosition(G4ThreeVector(0,0,0), 7.5*cm, 1.5*m, vertex_position);
	} else if(positioner=="DecayTrapFiducial") {
		// uniform fill of 5.5cm radius volume from which events reach detectors
		randomTubePosition(G4ThreeVector(0,0,0), 5.5*cm, 1.5*m, vertex_position);
	} else if(positioner=="SpectrometerVolumeUniform") {
		// uniform fill of (visible) spectrometer volume, for Xe
		randomTubePosition(G4ThreeVector(0,0,0), 7.5*cm, 2.17*m, vertex_position);
	} else {
		G4cout << "********* WARNING: Undefined positioner type! Defaulting to 'Fixed'! **********" << G4endl;
	}
	particleGun->SetParticlePosition(vertex_position);
	
	if(gunType=="Sn113") {
		Sn113SourceGenerator(anEvent);
	} else if(gunType=="Sn113G") {
		Sn113GSourceGenerator(anEvent);
	} else if(gunType=="Bi207") {
		Bi207SourceGenerator(anEvent);
	} else if(gunType=="Bi207G") {
		Bi207GSourceGenerator(anEvent);
	} else if(gunType=="Cd109") {
		funcBetaSpectrum->FixParameter(0,579.54);
		Cd109SourceGenerator(anEvent);
	} else if(gunType=="Cd109G") {
		Cd109GSourceGenerator(anEvent);
	} else if(gunType=="Cd113mG") {
		Cd113mGSourceGenerator(anEvent);
	} else if(gunType=="Sr85") {
		funcBetaSpectrum->FixParameter(0,1495.1);//assume Sr89 impurity
		Sr85SourceGenerator(anEvent);
	} else if(gunType=="Mn54") {
		Mn54SourceGenerator(anEvent);
	} else if(gunType=="Ce139") {
		Ce139SourceGenerator(anEvent);
	} else if(gunType=="Ce139G") {
		Ce139GSourceGenerator(anEvent);
	} else if(gunType=="Xe125_1/2+") {
		Xe125_1_2p_SourceGenerator(anEvent);
	} else if(gunType=="Xe133_11/2-") {
		Xe133_11_2m_SourceGenerator(anEvent);
	} else if(gunType=="Xe133_3/2+") {
		Xe133_3_2p_SourceGenerator(anEvent);
	} else if(gunType=="Xe135_11/2-") {
		Xe135_11_2m_SourceGenerator(anEvent);
	} else if(gunType=="Xe135_3/2+") {
		Xe135_3_2p_SourceGenerator(anEvent);
	} else if(gunType=="Xe137_7/2-") {
		Xe137_7_2m_SourceGenerator(anEvent);
	} else if (gunType=="In114") {
		funcBetaSpectrum->FixParameter(0,1988.7);
		In114SourceGenerator(anEvent);
	} else if (gunType=="In114G") {
		In114GSourceGenerator(anEvent);
	} else if (gunType=="endpoint" || gunType=="neutronBetaUnpol") {
		double eOrig = particleGun->GetParticleEnergy();
		if(gunType=="neutronBetaUnpol") {
			funcBetaSpectrum->FixParameter(0,neutronBetaEp);
		} else {
			// assume endpoint is set by particle gun energy if not for neutrons
			funcBetaSpectrum->FixParameter(0,eOrig/keV);
		}
		particleGun->SetParticleEnergy(funcBetaSpectrum->GetRandom()*keV);
		G4ThreeVector direction;
		RandomizeMomentum(direction);
		particleGun->SetParticleMomentumDirection(direction);
		displayGunStatus();
		particleGun->GeneratePrimaryVertex(anEvent);
		particleGun->SetParticleEnergy(eOrig);
	} else if(gunType=="uniformRandMomentum") {
		double eOrig = particleGun->GetParticleEnergy();
		particleGun->SetParticleEnergy(G4UniformRand()*eOrig);
		G4ThreeVector direction;
		RandomizeMomentum(direction);
		particleGun->SetParticleMomentumDirection(direction);
		displayGunStatus();
		particleGun->GeneratePrimaryVertex(anEvent);
		particleGun->SetParticleEnergy(eOrig);
	} else if (gunType=="eGunRandMomentum") {
		//for both eGun generator, assume gun energy is set by the
		//standard gun energy command!!!
		G4ThreeVector direction;
		RandomizeMomentum(direction);
		particleGun->SetParticleMomentumDirection(direction);
		displayGunStatus();
		particleGun->GeneratePrimaryVertex(anEvent);
	} else if (gunType=="eGun") {
		//for both eGun generator, assume gun energy is set by the   
		//standard gun energy command!!!
		particleGun->SetParticleMomentumDirection(G4ThreeVector(0,0,1));
		displayGunStatus();
		particleGun->GeneratePrimaryVertex(anEvent);
	} else if (gunType=="neutronBetaOff"){
		polarizedNeutronBetaDecayGenerator(anEvent,false);
	} else if (gunType=="neutronBetaOn"){
		polarizedNeutronBetaDecayGenerator(anEvent,true);
	} else {
		G4cout << "********* WARNING: Undefined generator type! No events generated!! **********" << G4endl;
	}
	
	
	//now fill the primary tree here
	gbmAnalysisManager->FillPrimaryData(anEvent,myseed);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

