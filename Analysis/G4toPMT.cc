#include "G4toPMT.hh"
#include "PathUtils.hh"
#include "BetaSpectrum.hh"
#include <cassert>
#include <cmath>

TRandom3 mc_rnd_source;	

Sim2PMT::Sim2PMT(const std::string& treeName): ProcessedDataScanner(treeName,false),
reSimulate(true), afp(AFP_OTHER) {
	for(Side s = EAST; s <= WEST; ++s) {
		PGen[s].setSide(s);
		PGen[s].larmorField = 0;
	}
	totalTime = BlindTime(1.0);
	fPID = PID_BETA;	// only simulating beta events
}

Stringmap Sim2PMT::evtInfo() {
	Stringmap m = ProcessedDataScanner::evtInfo();
	m.insert("Eprim",ePrim);
	m.insert("costh",costheta);
	m.insert("weight",physicsWeight);
	return m;
}
	
void Sim2PMT::calcReweight() {
	physicsWeight = 1.0; //= spectrumCorrectionFactor(ePrim) for beta spectrum
	if(afp==AFP_ON||afp==AFP_OFF)
		physicsWeight *= 1.0+correctedAsymmetry(ePrim,costheta*(afp==AFP_OFF?1:-1));
}

void Sim2PMT::setCalibrator(PMTCalibrator& PCal) {
	for(Side s = EAST; s <= WEST; ++s)
		PGen[s].setCalibrator(&PCal);
	evtRun = PCal.rn;
	ActiveCal = &PCal;
}

bool Sim2PMT::nextPoint() {
	bool np = ProcessedDataScanner::nextPoint();
	reverseCalibrate();
	calcReweight();
	nSimmed++;
	return np;
}

void Sim2PMT::reverseCalibrate() {
	
	doUnits();
	
	bool passesMWPC[2];
	bool passesScint[2];
	bool is2fold[2];
	
	// simulate event on both sides
	for(Side s = EAST; s <= WEST; ++s) {
		
		// TODO efficiency/resolution fancier wirechamber simulation
		mwpcEnergy[s] = eW[s];
		passesMWPC[s] = (eW[s] > 0.02);		
		
		// simulate detector energy response, or use un-smeared original data 
		if(reSimulate) {
			PGen[s].setOffsets(scintPos[s][X_DIRECTION], scintPos[s][Y_DIRECTION], wires[s][X_DIRECTION].center, wires[s][Y_DIRECTION].center);
			scints[s] = PGen[s].generate(eQ[s]);
			passesScint[s] = PGen[s].triggered();
		} else {
			scints[s].energy = eQ[s];
			passesScint[s] = (eQ[s] > 0);
		}
		
		is2fold[s] = passesMWPC[s] && passesScint[s];
	}
	
	if(is2fold[EAST] || is2fold[WEST]) fPID = PID_BETA;
	else fPID = PID_SINGLE;
	
	fType = TYPE_IV_EVENT;
	fSide = NONE;
	for(Side s = EAST; s<=WEST; ++s) {
		if(is2fold[EAST] && is2fold[WEST]) fType = TYPE_I_EVENT;
		if(is2fold[s] && !passesScint[otherSide(s)]) fType = passesMWPC[otherSide(s)]?TYPE_II_EVENT:TYPE_0_EVENT;
		if(passesScint[s]&&!passesScint[otherSide(s)]) fSide = s;
	}
	if(passesScint[EAST] && passesScint[WEST])
		fSide = time[EAST]<time[WEST]?EAST:WEST;
}

float Sim2PMT::getEtrue() {
	if(fSide>WEST) return 0;
	return PGen[fSide].getCalibrator()->Etrue(fSide,fType,scints[EAST].energy.x,scints[WEST].energy.x);
}



void G4toPMT::setReadpoints() {
	Tch->SetBranchAddress("EdepQ",eQ);
	Tch->SetBranchAddress("Edep",eDep);
	Tch->SetBranchAddress("MWPCEnergy",eW);
	Tch->SetBranchAddress("ScintPos",scintPos);
	Tch->SetBranchAddress("MWPCPos",mwpcPos);
	Tch->SetBranchAddress("time",time);
	Tch->SetBranchAddress("primTheta",&costheta);
	Tch->SetBranchAddress("primKE",&ePrim);
	if(Tch->GetBranch("primPos"))
		Tch->SetBranchAddress("primPos",primPos);
	else
		primPos[0] = primPos[1] = primPos[2] = primPos[3] = 0;
}

void PenelopeToPMT::setReadpoints() {
	Tch->SetBranchAddress("Epe",&fEdep[EAST]);
	Tch->SetBranchAddress("Epw",&fEdep[WEST]);
	Tch->SetBranchAddress("Egea",&fEW[EAST]);
	Tch->SetBranchAddress("Egwa",&fEW[WEST]);
	
	Tch->SetBranchAddress("eposx",&fMWPCpos[EAST][X_DIRECTION]);
	Tch->SetBranchAddress("eposy",&fMWPCpos[EAST][Y_DIRECTION]);
	Tch->SetBranchAddress("wposx",&fMWPCpos[WEST][X_DIRECTION]);
	Tch->SetBranchAddress("wposy",&fMWPCpos[WEST][Y_DIRECTION]);
	
	Tch->SetBranchAddress("te",&fTime[EAST]);
	Tch->SetBranchAddress("tw",&fTime[WEST]);
	
	Tch->SetBranchAddress("X",&fPrimPos[X_DIRECTION]);
	Tch->SetBranchAddress("Y",&fPrimPos[Y_DIRECTION]);
	Tch->SetBranchAddress("Z",&fPrimPos[Z_DIRECTION]);
	
	Tch->SetBranchAddress("E",&fEprim);
	Tch->SetBranchAddress("W",&fCostheta);
}

void G4toPMT::doUnits() {
	const double wcPosConversion = 10.0*sqrt(0.6);
	for(unsigned int i=0; i<3; i++)
		primPos[i] *= 10.0;
	for(Side s = EAST; s <= WEST; ++s) {
		if(matchPenelope)
			eQ[s] = eDep[s];
		for(AxisDirection d=X_DIRECTION; d<=Y_DIRECTION; ++d) {
			mwpcPos[s][d] *= wcPosConversion;
			scintPos[s][d] *= wcPosConversion;
			if(matchPenelope)
				scintPos[s][d] = wires[s][d].center = mwpcPos[s][d] = primPos[d];
			wires[s][d].center = mwpcPos[s][d];
		}
	}
	costheta=cos(costheta);
}

void PenelopeToPMT::doUnits() {
	const double wcPosConversion = 10.0;
	for(unsigned int i=0; i<3; i++)
		primPos[i] = fPrimPos[i]*wcPosConversion;
	for(Side s = EAST; s <= WEST; ++s) {
		eQ[s] = fEdep[s]*0.001;	// really Edep and not Equenched
		eW[s] = fEW[s]*0.001;
		time[s] = fTime[s];
		for(AxisDirection d=X_DIRECTION; d<=Y_DIRECTION; ++d) {			
			scintPos[s][d] = primPos[d];	// fake scintillator pos from primary
			mwpcPos[s][d]  = primPos[d]; 	// same position in scintillator
			wires[s][d].center = mwpcPos[s][d];
		}
	}
	ePrim = fEprim*0.001;
	costheta = fCostheta;
}


//-------------------------------------------


void MixSim::startScan(unsigned int startRandom) {
	for(std::vector<Sim2PMT*>::iterator it = subSims.begin(); it != subSims.end(); it++)
		(*it)->startScan(startRandom);
}

bool MixSim::nextPoint() {
	assert(cumStrength.size());
	std::vector<double>::iterator sline = std::lower_bound(cumStrength.begin(),cumStrength.end(),mc_rnd_source.Uniform(0,cumStrength.back()));
	unsigned int i = sline-cumStrength.begin();
	assert(i<subSims.size());
	currentSim = subSims[i];
	currentSim->nextPoint();
	
	// copy over simulation info
	for(Side s = EAST; s <= WEST; ++s) {
		eQ[s]=currentSim->eQ[s];
		eDep[s]=currentSim->eDep[s];
		eW[s]=currentSim->eW[s];
		time[s]=currentSim->time[s];
		for(unsigned int a=0; a<3; a++) {
			scintPos[s][a]=currentSim->scintPos[s][a];
			mwpcPos[s][a]=currentSim->mwpcPos[s][a];
		}
	}
	for(unsigned int a=0; a<4; a++)
		primPos[a] = currentSim->primPos[a];
	costheta = currentSim->costheta;
	ePrim = currentSim->ePrim;
	physicsWeight = currentSim->physicsWeight;
	
	// copy over data info
	for(Side s = EAST; s <= WEST; ++s) {
		scints[s] = currentSim->scints[s];
		led_pd[s] = currentSim->led_pd[s];
		mwpcs[s] = currentSim->mwpcs[s];
		mwpcEnergy[s] = currentSim->mwpcEnergy[s];
		for(unsigned int d = X_DIRECTION; d <= Y_DIRECTION; d++)
			wires[s][d] = currentSim->wires[s][d];
	}
	runClock = currentSim->runClock;
	fPID = currentSim->fPID;
	fType = currentSim->fType;
	fSide = currentSim->fSide;
	
	return true;
}

void MixSim::addSim(Sim2PMT* S, double r0, double thalf) {
	subSims.push_back(S);
	initStrength.push_back(r0);
	halflife.push_back(thalf);
	cumStrength.push_back((cumStrength.size()?cumStrength.back():0)+exp((t0-t1)*log(2)/thalf)*r0);
}

void MixSim::setTime(double t) {
	for(unsigned int i=0; i<halflife.size(); i++)
		cumStrength[i] = (i?cumStrength[i-1]:0)+exp((t0-t1)*log(2)/halflife[i])*initStrength[i];
}

void MixSim::setAFP(AFPState a) {
	for(std::vector<Sim2PMT*>::iterator it = subSims.begin(); it != subSims.end(); it++)
		(*it)->setAFP(a);
	Sim2PMT::setAFP(a);
}

void MixSim::setCalibrator(PMTCalibrator& PCal) {
	for(std::vector<Sim2PMT*>::iterator it = subSims.begin(); it != subSims.end(); it++)
		(*it)->setCalibrator(PCal);
	Sim2PMT::setCalibrator(PCal);	
}

