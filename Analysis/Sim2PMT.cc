#include "Sim2PMT.hh"
#include "PathUtils.hh"
#include "BetaSpectrum.hh"
#include <cassert>
#include <cmath>

TRandom3 mc_rnd_source;	


//--------------------------------------------------------------

void SimPositioner::applyOffset(Sim2PMT& S) {
	calcOffset(S);
	for(AxisDirection d = X_DIRECTION; d <= Y_DIRECTION; ++d) {
		for(Side s = EAST; s <= WEST; ++s) {
			S.scintPos[s][d] += offPos[d];
			S.mwpcPos[s][d] += offPos[d];
		}
		S.primPos[d] += offPos[d];
	}
}

void SourcedropPositioner::calcOffset(const Sim2PMT& S) {
	while(true) {
		offPos[X_DIRECTION] = mc_rnd_source.Uniform(-1.,1.);
		offPos[Y_DIRECTION] = mc_rnd_source.Uniform(-1.,1.);
		if(offPos[X_DIRECTION]*offPos[X_DIRECTION]+offPos[Y_DIRECTION]*offPos[Y_DIRECTION]<=1.) break;
	}
	offPos[X_DIRECTION] *= r0;
	offPos[Y_DIRECTION] *= r0;
	offPos[X_DIRECTION] += x0-S.primPos[X_DIRECTION];
	offPos[Y_DIRECTION] += y0-S.primPos[Y_DIRECTION];
}

//--------------------------------------------------------------


Sim2PMT::Sim2PMT(const std::string& treeName): ProcessedDataScanner(treeName,false),
SP(NULL), reSimulate(true), fakeClip(false), weightAsym(true),
nSimmed(0), nCounted(0), mwpcAccidentalProb(0), afp(AFP_OTHER) {
	for(Side s = EAST; s <= WEST; ++s) {
		PGen[s].setSide(s);
		mwpcThresh[s] = 0; //(s==EAST?0.98:0.10);
		for(AxisDirection d = X_DIRECTION; d <= Y_DIRECTION; ++d)
			for(unsigned int c = 0; c < kMaxCathodes; c++)
				cathodes[s][d][c] = 0;
		fTaggedBack[s] = false;
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
	physicsWeight = 1.0;
	if(weightAsym && (afp==AFP_ON || afp==AFP_OFF))
		physicsWeight *= (1.0+correctedAsymmetry(ePrim,costheta*(afp==AFP_ON?1:-1)))*spectrumCorrectionFactor(ePrim);
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
	if(SP) SP->applyOffset(*this);
	for(AxisDirection d = X_DIRECTION; d <= Y_DIRECTION; ++d)
		for(Side s = EAST; s <= WEST; ++s)
			wires[s][d].rawCenter = wires[s][d].center = mwpcPos[s][d];
	
	// simulate event on both sides
	for(Side s = EAST; s <= WEST; ++s) {
		// wirechamber accidentals
		if(!eW[s] && mc_rnd_source.Uniform(0.,1.)<mwpcAccidentalProb)
			eW[s] = mwpcThresh[s]+0.1;
		
		mwpcEnergy[s] = eW[s];
		// simulate detector energy response, or use un-smeared original data 
		if(reSimulate) {
			PGen[s].evtm = runClock[BOTH];
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
	
	primSide = costheta>0?WEST:EAST;
	
	for(Side s = EAST; s <= WEST; ++s) {
		passesMWPC[s] = (mwpcEnergy[s] > mwpcThresh[s]);		
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
	
	// Type II/III separation
	fProbIII = ((fType==TYPE_II_EVENT)?
				WirechamberCalibrator::sep23Prob(fSide,getEnergy(),mwpcEnergy[fSide])
				: 0.);
	if(fProbIII>0.5) fType=TYPE_III_EVENT;
}

float Sim2PMT::getEtrue() {
	if(fSide>WEST) return 0;
	return PGen[fSide].getCalibrator()->Etrue(fSide,fType,scints[EAST].energy.x,scints[WEST].energy.x);
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
		for(AxisDirection d = X_DIRECTION; d <= Y_DIRECTION; ++d)
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

