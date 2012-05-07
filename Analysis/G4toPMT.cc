#include "G4toPMT.hh"
#include "PathUtils.hh"
#include "BetaSpectrum.hh"
#include <cassert>
#include <cmath>

TRandom3 mc_rnd_source;	

Sim2PMT::Sim2PMT(const std::string& treeName): ProcessedDataScanner(treeName,false),
reSimulate(true), fakeClip(false), nSimmed(0), nCounted(0), mwpcAccidentalProb(4.3e-4), afp(AFP_OTHER) {
	offPos[X_DIRECTION] = offPos[Y_DIRECTION] = 0.;
	for(Side s = EAST; s <= WEST; ++s) {
		PGen[s].setSide(s);
		mwpcThresh[s] = 0.02;
	}
	totalTime = BlindTime(1.0);
	fPID = PID_BETA;	// only simulating beta events
}

void Sim2PMT::setOffset(double dx, double dy) {
	printf("Set simulation position offset to (%.2f,%.2f)\n",dx,dy);
	offPos[X_DIRECTION]=dx;
	offPos[Y_DIRECTION]=dy;
}

Stringmap Sim2PMT::evtInfo() {
	Stringmap m = ProcessedDataScanner::evtInfo();
	m.insert("Eprim",ePrim);
	m.insert("costh",costheta);
	m.insert("weight",physicsWeight);
	return m;
}

void Sim2PMT::calcReweight() {
	physicsWeight = 1.0;
	if(afp==AFP_ON||afp==AFP_OFF)
		physicsWeight *= (1.0+correctedAsymmetry(ePrim,costheta*(afp==AFP_OFF?1:-1)))*spectrumCorrectionFactor(ePrim);
	if(fakeClip) {
		const double R = 70.*sqrt(0.6); // wirechamber entrance window radius, projected back to decay trap
		// event origin distance from edge
		double l = R-sqrt(primPos[X_DIRECTION]*primPos[X_DIRECTION]+primPos[Y_DIRECTION]*primPos[Y_DIRECTION]);	
		if(l<=0) { physicsWeight=0; return; }
		double pt = sqrt(1-costheta*costheta)*ePrim; // transverse momentum
		double r = pt/(300.0*1.0); // larmor radius, same in decay trap as projected back from window
		if(l>=2*r) return;
		double cosalpha = (l*l+2.*R*(r-l))/(2.*r*(R-l));
		physicsWeight *= acos(cosalpha)/3.1415926535;
		if(!(physicsWeight==physicsWeight)) physicsWeight=0;
	}
}

void Sim2PMT::setCalibrator(PMTCalibrator& PCal) {
	for(Side s = EAST; s <= WEST; ++s)
		PGen[s].setCalibrator(&PCal);
	ActiveCal = &PCal;
	evtRun = ActiveCal->rn;
}

bool Sim2PMT::nextPoint() {
	bool np = ProcessedDataScanner::nextPoint();
	reverseCalibrate();
	calcReweight();
	nSimmed++;
	nCounted+=simEvtCounts();
	return np;
}

void Sim2PMT::reverseCalibrate() {
	
	doUnits();
	
	assert(ActiveCal);
	evtRun = ActiveCal->rn;
	
	// simulated event time stamp
	runClock = mc_rnd_source.Uniform(0.,ActiveCal->totalTime);
	
	// apply position offsets; set wires position
	for(AxisDirection d = X_DIRECTION; d <= Y_DIRECTION; ++d) {
		for(Side s = EAST; s <= WEST; ++s) {
			scintPos[s][d] += offPos[d];
			mwpcPos[s][d] += offPos[d];
			wires[s][d].center = mwpcPos[s][d];
		}
		primPos[d] += offPos[d];
	}
	
	// simulate event on both sides
	for(Side s = EAST; s <= WEST; ++s) {
		// wirechamber accidentals
		if(!eW[s] && mc_rnd_source.Uniform(0.,1.)<mwpcAccidentalProb)
			eW[s] = mwpcThresh[s]+0.1;
		
		mwpcEnergy[s] = eW[s];
		// simulate detector energy response, or use un-smeared original data 
		if(reSimulate) {
			PGen[s].evtm = runClock.t[BOTH];
			PGen[s].setPosition(scintPos[s][X_DIRECTION], scintPos[s][Y_DIRECTION],
								wires[s][X_DIRECTION].center-scintPos[s][X_DIRECTION], wires[s][Y_DIRECTION].center-scintPos[s][Y_DIRECTION]);
			scints[s] = PGen[s].generate(eQ[s]);
			passesScint[s] = PGen[s].triggered();
		} else {
			scints[s].energy = eQ[s];
			passesScint[s] = (eQ[s] > 0);
		}
	}
	
	classifyEvent();
}

void Sim2PMT::classifyEvent() {
	
	bool passesMWPC[2];
	bool is2fold[2];
	
	for(Side s = EAST; s <= WEST; ++s) {
		passesMWPC[s] = (eW[s] > mwpcThresh[s]);		
		is2fold[s] = passesMWPC[s] && passesScint[s];
	}
	
	if(is2fold[EAST] || is2fold[WEST]) fPID = PID_BETA;
	else fPID = PID_SINGLE;
	
	fType = TYPE_IV_EVENT;
	fSide = NOSIDE;
	for(Side s = EAST; s<=WEST; ++s) {
		if(is2fold[s]) {
			if(passesScint[otherSide(s)])
				fType = TYPE_I_EVENT;
			else
				fType = passesMWPC[otherSide(s)]?TYPE_II_EVENT:TYPE_0_EVENT;
		}
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
	// wirechamber position projection plus empirical window-diameter-matching fudge factor
	const double wcPosConversion = 10.0*sqrt(0.6)*(51.96/52.8);
	for(unsigned int i=0; i<3; i++)
		primPos[i] *= 10.0;
	for(Side s = EAST; s <= WEST; ++s) {
		for(AxisDirection d=X_DIRECTION; d<=Y_DIRECTION; ++d) {
			mwpcPos[s][d] *= wcPosConversion;
			scintPos[s][d] *= wcPosConversion;
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

void G4toPMT_SideSwap::doUnits() {
	std::swap(eQ[EAST],eQ[WEST]);
	std::swap(eDep[EAST],eDep[WEST]);
	std::swap(eW[EAST],eW[WEST]);
	for(AxisDirection d=X_DIRECTION; d<=Y_DIRECTION; ++d) {
		std::swap(mwpcPos[EAST][d],mwpcPos[WEST][d]);
		std::swap(scintPos[EAST][d],scintPos[WEST][d]);
	}
	G4toPMT::doUnits();
}

//-------------------------------------------


void MixSim::startScan(bool startRandom) {
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

unsigned int MixSim::getnFiles() const {
	unsigned int nf = 0;
	for(std::vector<Sim2PMT*>::const_iterator it = subSims.begin(); it != subSims.end(); it++)
		nf += (*it)->getnFiles();
	return nf;
}
