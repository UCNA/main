#include "PMTGenerator.hh"
#include "QFile.hh"
#include "strutils.hh"
#include "PathUtils.hh"
#include <algorithm>
#include <Math/DistFunc.h>

TRandom3 PMTGenerator::sim_rnd_source;

PMTGenerator::PMTGenerator(Side s, float xx, float yy):
x(xx), y(yy), xw(xx), yw(yy), evtm(0), presmear(0), pedcorr(0.4), mySide(s) { }

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

ScintEvent PMTGenerator::generate(float en) {
	float nPEtot = 0;
	float restot = 0;
	float tubeRes;
	if(en<=0) {
		for(unsigned int t=0; t<nBetaTubes; t++)
			sevt.tuben[t]=sevt.adc[t]=0;
		sevt.energy=0;
		return sevt;
	}
	// PE-smeared energy
	for(unsigned int t=0; t<nBetaTubes; t++) {
		if(currentCal->scaleNoiseWithL)
			tubeRes = pmtRes[mySide][t]*en;
		else
			tubeRes = currentCal->nPE(mySide,t,en,x,y,evtm);
		if(presmear) {
			if(en*presmear > tubeRes+0.1)
				tubeRes = 1.0/(1.0/tubeRes-1.0/(en*presmear));
			else
				tubeRes = 1.0/(1.0/tubeRes-1.0/(tubeRes+0.1));
		}
		restot += tubeRes;
		float nPE = sim_rnd_source.PoissonD(tubeRes);      // photoelectrons produced
		nPEtot += nPE;
		sevt.tuben[t] = nPE/tubeRes*en;
	}
	
	// Calculations for correlated pedestals
	// target pedestal noise widths
	float pedw[nBetaTubes];
	for(unsigned int t=0; t<nBetaTubes; t++)
		pedw[t] = currentCal->getPedwidth(currentCal->sensorNames[mySide][t],evtm);
	// generate uncorrelated fluctuations
	float pedc[nBetaTubes];
	float d = (1.+(nBetaTubes-1)*pedcorr)*(1.-pedcorr);
	float od = -pedcorr/d;
	d = ((nBetaTubes-2)*pedcorr+1)/d;
	for(unsigned int t=0; t<nBetaTubes; t++) {
		pedc[t] = 0;
		for(unsigned int t2=0; t2<nBetaTubes; t2++)
			pedc[t] += pedw[t2]*pedw[t2]*(t2==t?d:od);
		pedc[t] = pedc[t]<0?0:sim_rnd_source.Gaus(0.,sqrt(pedc[t]));
	}
		   
	for(unsigned int t=0; t<nBetaTubes; t++) {
		// adc value
		sevt.adc[t] = currentCal->invertCorrections(mySide,t,sevt.tuben[t].x,x,y,evtm); 
		// correlated pedestal noise
		for(unsigned int t2=0; t2<nBetaTubes; t2++)
			sevt.adc[t] += pedc[t2]*(t2==t?1.:pedcorr);
		// quantization noise
		float pedx = currentCal->getPedestal(currentCal->sensorNames[mySide][t],evtm);
		sevt.adc[t] = int(sevt.adc[t]+pedx+0.5)-pedx;
	}
	currentCal->calibrateEnergy(mySide, xw, yw, sevt, evtm);
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

