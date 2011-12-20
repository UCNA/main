#include "G4toPMT.hh"
#include "PathUtils.hh"
#include "BetaSpectrum.hh"
#include <cassert>
#include <cmath>

Sim2PMT::Sim2PMT(const std::string& treeName): ProcessedDataScanner(treeName,false),
reSimulate(true), afp(AFP_OTHER) {
	for(Side s = EAST; s <= WEST; s = nextSide(s)) {
		PGen[s].setSide(s);
		PGen[s].larmorField = 0;
	}
	inputEnergy[EAST]=inputEnergy[WEST]=NULL;
	totalTime = BlindTime(1.0);
	fPID = PID_BETA;	// only simulating beta events
}

void Sim2PMT::calcReweight() {
	physicsWeight = 1.0; //= spectrumCorrectionFactor(ePrim) for beta spectrum
	if(afp==AFP_ON||afp==AFP_OFF)
		physicsWeight *= 1.0+correctedAsymmetry(ePrim,costheta*(afp==AFP_OFF?1:-1));
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
	for(Side s = EAST; s <= WEST; s = nextSide(s)) {
		if(matchPenelope)
			eQ[s] = eDep[s];
		for(unsigned int i=0; i<2; i++) {
			mwpcPos[s][i] *= wcPosConversion;
			scintPos[s][i] *= wcPosConversion;
			if(matchPenelope)
				scintPos[s][i] = wires[s][i].center = mwpcPos[s][i] = primPos[i];
			wires[s][i].center = mwpcs[s].pos[i]=mwpcPos[s][i];
		}
	}
	costheta=cos(costheta);
}

void PenelopeToPMT::doUnits() {
	const double wcPosConversion = 10.0;
	for(unsigned int i=0; i<3; i++)
		primPos[i] = fPrimPos[i]*wcPosConversion;
	for(Side s = EAST; s <= WEST; s = nextSide(s)) {
		eQ[s] = fEdep[s]*0.001;	// really Edep and not Equenched
		eW[s] = fEW[s]*0.001;
		time[s] = fTime[s];
		for(unsigned int i=0; i<2; i++) {			
			scintPos[s][i] = primPos[i];	// fake scintillator pos from primary
			mwpcPos[s][i]  = primPos[i]; 	// same position in scintillator
			wires[s][i].center = mwpcs[s].pos[i] = mwpcPos[s][i];
		}
	}
	ePrim = fEprim*0.001;
	costheta = fCostheta;
}

void Sim2PMT::setCalibrator(PMTCalibrator& PCal) {
	for(Side s = EAST; s <= WEST; s = nextSide(s))
		PGen[s].setCalibrator(&PCal);
	ActiveCal = &PCal;
}

bool Sim2PMT::nextPoint() {
	bool np = ProcessedDataScanner::nextPoint();
	reverseCalibrate();
	calcReweight();
	return np;
}

void Sim2PMT::reverseCalibrate() {
	
	doUnits();
	
	bool passesMWPC[2];
	bool passesScint[2];
	bool is2fold[2];
	
	// simulate event on both sides
	for(Side s = EAST; s <= WEST; s = nextSide(s)) {
		
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
	for(Side s = EAST; s<=WEST; s=nextSide(s)) {
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

void Sim2PMT::genType0(unsigned int nToSim, double wx, double wy) {
	
	assert(Tch->GetEntries());
	
	printf("Smearing simulated events (target %i)...\n",nToSim);
	if(wx && wy)
		printf("\twith position cut (%g,%g)\n",wx,wy);
	unsigned int nsimmed[2] = {0,0};
	
	startScan(nToSim);
	
	while(1) {
		// run until quota of points generated
		if(!nextPoint())
			if(!nToSim)
				break;
		if(nToSim && nsimmed[EAST] >= nToSim && nsimmed[WEST] >= nToSim)
			break;
		
		recalibrateEnergy();
	
		// save type-0 events
		for(Side s = EAST; s <= WEST; s = nextSide(s)) {
			bool wc_cut = eW[s] > 0.05 && eW[otherSide(s)] < 0.05;
			if(inputEnergy[s] && wc_cut && !eQ[otherSide(s)])
				inputEnergy[s]->Fill(eQ[s]);
			if(wx && wy && !(wires[s][X_DIRECTION].center*wires[s][X_DIRECTION].center/wx/wx
							 + wires[s][Y_DIRECTION].center*wires[s][Y_DIRECTION].center/wy/wy <= 1.0))
				continue;
			if(fType == TYPE_0_EVENT && fSide == s) {
				nsimmed[s]++;
			}
		}
	}
	
	printf(" Done [%i, %i / %i completed].\n",nsimmed[EAST],nsimmed[WEST],(int)Tch->GetEntries());
}
