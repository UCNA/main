#include "SimNonlinearity.hh"
#include <cassert>
#include <TRandom3.h>

TRandom3 snl_rnd_source; //< random number generator

SimNonlinearizer::SimNonlinearizer(): newrunGainerr(0), afpCorrGain(1.0), range(1500) {
	for(Side s = EAST; s <= WEST; s = nextSide(s))
		for(unsigned int t=0; t<nBetaTubes; t++)
			gainfactor[s][t]=1.0;
	setReserr(1.0);
	
	snl_rnd_source.SetSeed(0);	// randomize random seed set from system clock
	
	// error bounds PRL10
	/*
	errlim = new TGraph(6);
	errlim->SetPoint(0,0.0,5.0);
	errlim->SetPoint(1,100.0,5.0);
	errlim->SetPoint(2,250,250*0.02);
	errlim->SetPoint(3,500,500*0.013);
	errlim->SetPoint(4,900,900*0.025);
	errlim->SetPoint(5,1000,1000*0.025);
	*/
	
	// error bounds 2011
	errlim = new TGraph(4);
	errlim->SetPoint(0,0.0,2.5);
	errlim->SetPoint(1,200,200*0.0125);
	errlim->SetPoint(2,500,500*0.0125);
	errlim->SetPoint(3,1000,500*0.0125);
}

float SimNonlinearizer::delinearize(Side s, unsigned int t, float l0) const {
	assert(s==EAST || s==WEST);
	assert(t<nBetaTubes);
	float lout = l0;
	for(unsigned int n=0; n<abscoeffs[s][t].size(); n++)
		lout += cos(3.1415926535*n*l0/range)*abscoeffs[s][t][n];
	for(unsigned int n=0; n<relcoeffs[s][t].size(); n++)
		lout += cos(3.1415926535*n*l0/range)*relcoeffs[s][t][n]*l0;
	return lout*gainfactor[s][t];
}

void SimNonlinearizer::makeRanderr(bool correlate) {
	do {
		randomAbscurve(5.0,5,false);
		randomRelcurve(0.025,5,false);
		unifyErrors();	
		if(correlate)
			correlateErrors();
	} while(!checkLimits());
}

bool SimNonlinearizer::checkLimits() const {
	for(Side s = EAST; s <= WEST; s = nextSide(s))
		for(unsigned int t=0; t<nBetaTubes; t++)
			for(float x = 0; x <= 1500; x += 1500/50.0)
				if( fabs(delinearize(s, t, x)-x) > fabs(errlim->Eval(x)) )
					return false;
	return true;
}

float SimNonlinearizer::maxError(Side s, unsigned int t, float emin, float emax) const {
	float de = 0;
	for(float x = emin; x <= emax; x += (emax-emin)/20.0) {
		float d = delinearize(s,t,x)-x;
		if(fabs(d)>fabs(de))
			de=d;
	}
	return de;
}

float SimNonlinearizer::maxError(float emin, float emax) const {
	float de = 0;
	for(Side s = EAST; s <= WEST; s = nextSide(s)) {
		for(unsigned int t=0; t<nBetaTubes; t++) {
			float d = maxError(s,t,emin,emax);
			if(fabs(d)>fabs(de))
				de=d;
		}
	}
	return de;
}


void SimNonlinearizer::startNewRun(AFPState afp) {
	for(Side s = EAST; s <= WEST; s = nextSide(s)) {
		for(unsigned int t=0; t<nBetaTubes; t++) {
			if(afpCorrGain != 1.0) {
				if(afpCorrGain > 0) {
					if(afp == AFP_OFF)
						gainfactor[s][t] = 1.0/sqrt(afpCorrGain);
					if(afp == AFP_ON)
						gainfactor[s][t] = sqrt(afpCorrGain);
				} else { 
					if(afp == AFP_OFF) {
						gainfactor[EAST][t] = 1.0/sqrt(-afpCorrGain);
						gainfactor[WEST][t] = sqrt(-afpCorrGain);
					}
					if(afp == AFP_ON) {
						gainfactor[EAST][t] = sqrt(-afpCorrGain);
						gainfactor[WEST][t] = 1.0/sqrt(-afpCorrGain);
					}
				}
			}
			if(newrunGainerr)
				gainfactor[s][t]=snl_rnd_source.Gaus(1.0,newrunGainerr);		
		}
	}
	
}

Stringmap SimNonlinearizer::toStringmap() const {
	Stringmap m;
	for(Side s = EAST; s <= WEST; s = nextSide(s)) {
		for(unsigned int t=0; t<nBetaTubes; t++) {
			std::string st = ctos(sideNames(s))+itos(t);
			m.insert(std::string("gainfactor_")+st,gainfactor[s][t]);
			if(resErr[s][t] != 1.0)
				m.insert(std::string("reserr_")+st,resErr[s][t]);
			if(abscoeffs[s][t].size())
				m.insert(std::string("abscoeffs_")+st,vtos(abscoeffs[s][t]));
			if(relcoeffs[s][t].size())
				m.insert(std::string("relcoeffs_")+st,vtos(relcoeffs[s][t]));
		}
	}
	m.insert("afpCorrGain",afpCorrGain);
	m.insert("range",range);
	return m;
}

void SimNonlinearizer::fixCalPoint(Side s, unsigned int t, float l0, float l1) {
	if(l1==10000)
		l1 = l0;
	gainfactor[s][t]=1.0;
	gainfactor[s][t]=l1/delinearize(s,t,l0);
}

void SimNonlinearizer::setOffset(float dl, Side s0) {
	for(Side s=EAST; s<=WEST; s=nextSide(s)) {
		if(!(s0==BOTH || s==s0)) continue;
		for(unsigned int t=0; t<nBetaTubes; t++){
			abscoeffs[s][t].clear();
			abscoeffs[s][t].push_back(dl);
		}
	}
}

void SimNonlinearizer::randomOffsets(float sigma, float mu, bool gaussStats) {
	for(Side s=EAST; s<=WEST; s=nextSide(s)) {
		for(unsigned int t=0; t<nBetaTubes; t++){
			abscoeffs[s][t].clear();
			if(gaussStats)
				abscoeffs[s][t].push_back(snl_rnd_source.Gaus(mu,sigma));
			else
				abscoeffs[s][t].push_back(snl_rnd_source.Uniform(-sigma,sigma));
		}
	}
}

void SimNonlinearizer::randomGainflucts(float sigma, float mu, bool gaussStats) {
	for(Side s=EAST; s<=WEST; s=nextSide(s)) {
		for(unsigned int t=0; t<nBetaTubes; t++){
			relcoeffs[s][t].clear();
			if(gaussStats)
				relcoeffs[s][t].push_back(snl_rnd_source.Gaus(mu,sigma));
			else
				relcoeffs[s][t].push_back(snl_rnd_source.Uniform(-sigma,sigma));
		}
	}
}

void SimNonlinearizer::randomRelcurve(float sigmaRel, unsigned int nTerms, bool gaussStats) {
	for(Side s=EAST; s<=WEST; s=nextSide(s)) {
		for(unsigned int t=0; t<nBetaTubes; t++){
			relcoeffs[s][t].clear();
			for(unsigned int n=0; n<nTerms; n++) {
				if(gaussStats)
					relcoeffs[s][t].push_back(snl_rnd_source.Gaus(0,sigmaRel));
				else
					relcoeffs[s][t].push_back(snl_rnd_source.Uniform(-sigmaRel,sigmaRel));				
			}
		}
	}
}

void SimNonlinearizer::randomAbscurve(float sigma, unsigned int nTerms, bool gaussStats) {
	for(Side s=EAST; s<=WEST; s=nextSide(s)) {
		for(unsigned int t=0; t<nBetaTubes; t++){
			abscoeffs[s][t].clear();
			for(unsigned int n=0; n<nTerms; n++) {
				if(gaussStats)
					abscoeffs[s][t].push_back(snl_rnd_source.Gaus(0,sigma));
				else
					abscoeffs[s][t].push_back(snl_rnd_source.Uniform(-sigma,sigma));				
			}
		}
	}
}


void SimNonlinearizer::setRelerr(float relerr, Side s0) {
	for(Side s=EAST; s<=WEST; s=nextSide(s)) {
		if(!(s0==BOTH || s==s0)) continue;
		for(unsigned int t=0; t<nBetaTubes; t++){
			relcoeffs[s][t].clear();
			relcoeffs[s][t].push_back(relerr);
		}
	}
}	

void SimNonlinearizer::setReserr(float reserr, Side s0) {
	for(Side s=EAST; s<=WEST; s=nextSide(s)) {
		if(!(s0==BOTH || s==s0)) continue;
		for(unsigned int t=0; t<nBetaTubes; t++)
			resErr[s][t] = reserr;
	}
}	

void SimNonlinearizer::unifyErrors() {
	for(Side s=EAST; s<=WEST; s=nextSide(s)) {
		for(unsigned int t=1; t<nBetaTubes; t++) {
			relcoeffs[s][t] = relcoeffs[s][0];
			abscoeffs[s][t] = abscoeffs[s][0];
			resErr[s][t] = resErr[s][0];
			gainfactor[s][t] = gainfactor[s][0];
		}
	}	
}

void SimNonlinearizer::correlateErrors(bool correlate) {
	for(unsigned int t=0; t<nBetaTubes; t++) {
		relcoeffs[WEST][t] = relcoeffs[EAST][t];
		abscoeffs[WEST][t] = abscoeffs[EAST][t];
		if(correlate) {
			resErr[WEST][t] = resErr[EAST][t];
			gainfactor[WEST][t] = gainfactor[EAST][t];
		} else {
			for(std::vector<float>::iterator it = relcoeffs[WEST][t].begin(); it!=relcoeffs[WEST][t].end(); it++)
				(*it)*=-1.0;
			for(std::vector<float>::iterator it = abscoeffs[WEST][t].begin(); it!=abscoeffs[WEST][t].end(); it++)
				(*it)*=-1.0;
			resErr[WEST][t] = 1.0/resErr[EAST][t];
			gainfactor[WEST][t] = 1.0/gainfactor[EAST][t];			
		}
	}
}

void SimNonlinearizer::addGainNoise(float sigmaRel) {
	for(Side s=EAST; s<=WEST; s=nextSide(s))
		for(unsigned int t=0; t<nBetaTubes; t++)
			gainfactor[s][t] *= snl_rnd_source.Gaus(1.0,sigmaRel);
}

