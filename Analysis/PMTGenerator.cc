#include "PMTGenerator.hh"
#include "QFile.hh"
#include "strutils.hh"
#include "PathUtils.hh"
#include <TRandom3.h>
#include <algorithm>
#include <Math/DistFunc.h>

TRandom3 sim_rnd_source;		

PMTGenerator::PMTGenerator(Side s, float xx, float yy): calcADC(true),
x(xx), y(yy), xw(xx), yw(yy), presmear(0), mySide(s) { }

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
		float nPE = sim_rnd_source.PoissonD(tubeRes);
		nPE += sim_rnd_source.Gaus(0.0,0.7); // analog smoothing by PMT
		nPEtot += nPE;
		sevt.tuben[t] = nPE/tubeRes*en;
		
		if(calcADC) { // ADC + pedestal noise
			sevt.adc[t] = currentCal->invertCorrections(mySide,t,sevt.tuben[t].x,x,y,0); 
			float pedw = currentCal->getPedwidth(currentCal->sensorNames[mySide][t],0.0);
			sevt.adc[t] += sim_rnd_source.Gaus(0.0,pedw);
			//sevt.adc[t] += sim_rnd_source.Gaus(0.0,25);
			//sevt.adc[t] += (correrr*0.7+sim_rnd_source.Gaus(0.0,0.7))*20.0;
		}
	}
	if(calcADC)
		currentCal->calibrateEnergy(mySide, xw, yw, sevt, 0);
	else
		sevt.energy = nPEtot/restot*en;
	return sevt;
}

unsigned int PMTGenerator::triggers() {
	if(!calcADC) {
		nTrigs = nBetaTubes;
	} else {
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
	}
	return nTrigs;
}

