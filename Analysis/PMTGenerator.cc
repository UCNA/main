#include "PMTGenerator.hh"
#include "QFile.hh"
#include "strutils.hh"
#include "PathUtils.hh"
#include <TRandom3.h>
#include <algorithm>
#include <Math/DistFunc.h>

TRandom3 sim_rnd_source;		

PMTGenerator::PMTGenerator(Side s, float xx, float yy):
x(xx), y(yy), xw(xx), yw(yy), presmear(0), threshnoise(0.0), propnoise(0.04), mySide(s) { }

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
			pmtRes[s][t] = currentCal->nPE(s, t, 300.0, x, y, 0)/300.0;
}

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
	//double correrr = sim_rnd_source.Gaus(0.0,1.0);
	for(unsigned int t=0; t<nBetaTubes; t++) {
		if(currentCal->scaleNoiseWithL)
			tubeRes = pmtRes[mySide][t]*en;
		else
			tubeRes = currentCal->nPE(mySide,t,en,x,y,0);
		if(presmear) {
			if(en*presmear > tubeRes+0.1)
				tubeRes = 1.0/(1.0/tubeRes-1.0/(en*presmear));
			else
				tubeRes = 1.0/(1.0/tubeRes-1.0/(tubeRes+0.1));
		}
		restot += tubeRes;
		float nPE = sim_rnd_source.PoissonD(tubeRes); // photoelectrons produced
		nPE += sim_rnd_source.Gaus(0.0,threshnoise); // constant threshold noise
		nPE += sim_rnd_source.Gaus(0.0,tubeRes*propnoise); // proportional noise
		nPEtot += nPE;
		sevt.tuben[t] = nPE/tubeRes*en;
		
		// adc value
		sevt.adc[t] = currentCal->invertCorrections(mySide,t,sevt.tuben[t].x,x,y,0); 
		// pedestal and quantization noise
		float pedw = currentCal->getPedwidth(currentCal->sensorNames[mySide][t],0.0);
		sevt.adc[t] += sim_rnd_source.Gaus(0.,pedw);
		float pedx = currentCal->getPedestal(currentCal->sensorNames[mySide][t],0.0);
		float gmsf = currentCal->gmsFactor(mySide,t,0.0);
		sevt.adc[t] = (int(sevt.adc[t]/gmsf+pedx+0.5)-pedx)*gmsf;
	}
	currentCal->calibrateEnergy(mySide, xw, yw, sevt, 0);
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

