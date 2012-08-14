#include "PMTGenerator.hh"
#include "QFile.hh"
#include "strutils.hh"
#include "PathUtils.hh"
#include <algorithm>
#include <Math/DistFunc.h>

TRandom3 PMTGenerator::sim_rnd_source;

PMTGenerator::PMTGenerator(Side s, float xx, float yy):
x(xx), y(yy), xw(xx), yw(yy), evtm(0), presmear(0), dgain(16.0), pedcorr(0.2), crosstalk(0.010), xscatter(0.), mySide(s) { }

void PMTGenerator::setCalibrator(PMTCalibrator* P) { 
	assert(P);
	currentCal = P;
	setPosition(x,y,xw-x,yw-y);
}

void PMTGenerator::setPosition(float xx, float yy, float dxw, float dyw) {
	assert(currentCal);
	x = xx;
	y = yy;
	xw = x+dxw;
	yw = y+dyw;
	for(Side s = EAST; s <= WEST; ++s)
		for(unsigned int t=0; t<nBetaTubes; t++)
			pmtRes[s][t] = currentCal->nPE(s, t, 300.0, x, y, evtm)/300.0;
}

void PMTGenerator::setSide(Side s) { mySide = s; }

void preuncorrelate(float* v, double a, unsigned int n = nBetaTubes) {
	float d = (1.+(n-1)*a)*(1.-a);
	float od = -a/d;
	d = ((n-2)*a+1)/d;
	std::vector<float> v2(n);
	for(unsigned int i=0; i<n; i++)
		for(unsigned int j=0; j<n; j++)
			v2[i] += v[j]*(i==j?d:od);
	for(unsigned int i=0; i<n; i++)
		v[i]=v2[i];
}

void recorrelate(float* v, double a, unsigned int n = nBetaTubes) {
	std::vector<float> v2(n);
	for(unsigned int i=0; i<n; i++)
		for(unsigned int j=0; j<n; j++)
			v2[i] += v[j]*(i==j?1:a);
	for(unsigned int i=0; i<n; i++)
		v[i]=v2[i];
}

ScintEvent PMTGenerator::generate(float en) {
	if(en<=0) {
		for(unsigned int t=0; t<nBetaTubes; t++)
			sevt.tuben[t]=sevt.adc[t]=0;
		sevt.energy=0;
		return sevt;
	}
	
	// PE-smeared energy plus crosstalk
	float nPE[nBetaTubes];
	float tubeRes[nBetaTubes];
	for(unsigned int t=0; t<nBetaTubes; t++) {
		if(currentCal->scaleNoiseWithL)
			tubeRes[t] = pmtRes[mySide][t]*en;
		else
			tubeRes[t] = currentCal->nPE(mySide,t,en,x,y,evtm);
		if(presmear) {
			if(en*presmear > tubeRes[t]+0.1)
				tubeRes[t] = 1.0/(1.0/tubeRes[t]-1.0/(en*presmear));
			else
				tubeRes[t] = 1.0/(1.0/tubeRes[t]-1.0/(tubeRes[t]+0.1));
		}
		nPE[t] = tubeRes[t];
	}
	preuncorrelate(nPE, crosstalk);
	for(unsigned int t=0; t<nBetaTubes; t++)
		nPE[t] = sim_rnd_source.PoissonD(dgain*sim_rnd_source.PoissonD(nPE[t]>0?nPE[t]:0))/dgain;
	recorrelate(nPE, crosstalk);
	for(unsigned int t=0; t<nBetaTubes; t++)
		sevt.tuben[t] = nPE[t]/tubeRes[t]*en;
	
	
	// target pedestal noise widths squared
	float pedw[nBetaTubes];
	for(unsigned int t=0; t<nBetaTubes; t++) {
		pedw[t] = currentCal->getPedwidth(currentCal->sensorNames[mySide][t],evtm);
		pedw[t] *= pedw[t];
	}
	preuncorrelate(pedw,pedcorr);
	for(unsigned int t=0; t<nBetaTubes; t++)
		pedw[t] = pedw[t]<0?0:sim_rnd_source.Gaus(0.,sqrt(pedw[t]));
	recorrelate(pedw,pedcorr);
	
	for(unsigned int t=0; t<nBetaTubes; t++) {
		// adc value
		sevt.adc[t] = currentCal->invertCorrections(mySide,t,sevt.tuben[t].x,x,y,evtm); 
		// correlated pedestal noise
		sevt.adc[t] += pedw[t];
		// quantization noise
		float pedx = currentCal->getPedestal(currentCal->sensorNames[mySide][t],evtm);
		sevt.adc[t] = int(sevt.adc[t]+pedx+0.5)-pedx;
	}
	currentCal->calibrateEnergy(mySide, xw, yw, sevt, evtm, xscatter);
	return sevt;
}

unsigned int PMTGenerator::triggers() {
	
	nTrigs = 0;
	unsigned int nZero = 0;
	for(unsigned int t=0; t<nBetaTubes; t++) {
		pmtTriggered[t] = sim_rnd_source.Uniform(0.0,1.0) < currentCal->trigEff(mySide,t,sevt.adc[t]);
		if(sevt.adc[t] == 0) nZero++;
		if(pmtTriggered[t]) nTrigs++;
	}
	// check for all-zero events (probably negative input energy)
	if(nZero == nBetaTubes) {
		for(unsigned int t=0; t<nBetaTubes; t++)
			pmtTriggered[t] = false;
		nTrigs = 0;
	}
	return nTrigs;
}

