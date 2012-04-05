#include "PMTGenerator.hh"
#include "QFile.hh"
#include "strutils.hh"
#include "PathUtils.hh"
#include <TRandom3.h>
#include <algorithm>
#include <Math/DistFunc.h>

TRandom3 sim_rnd_source;		

PMTGenerator::PMTGenerator(Side s, float xx, float yy): calcADC(true),
x(xx), y(yy), dsx(0), dsy(0), dwx(0), dwy(0),
presmear(0), larmorField(0), mySide(s) { }

void PMTGenerator::setCalibrator(PMTCalibrator* P) { 
	assert(P);
	currentCal = P;
}

void PMTGenerator::setOffsets(float xxs, float yys, float xxw, float yyw) {
	dsx = xxs;
	dsy = yys;
	dwx = xxw;
	dwy = yyw;
	if(!(dsx==dsx)) dsx = -100;
	if(!(dsy==dsy)) dsy = -100;
	if(!(dwx==dwx)) dwx = -100;
	if(!(dwy==dwy)) dwy = -100;
	for(Side s = EAST; s <= WEST; ++s) {
		resTot[s] = 0;
		for(unsigned int t=0; t<nBetaTubes; t++) {
			pmtRes[s][t] = currentCal->nPE(s, t, 300.0, x+dsx, y+dsy, 0)/300.0;
			resTot[s] += pmtRes[s][t];
		}
	}
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
	if(larmorField) {
		float pz = sim_rnd_source.Uniform(0.0,en);
		float pxy = sqrt(en*en-pz*pz);
		float rg = pxy/(300.0*larmorField);
		float theta0 = sim_rnd_source.Uniform(0.0,2*PI);
		float theta1 = sim_rnd_source.Uniform(0.0,2*PI);
		float theta2 = sim_rnd_source.Uniform(0.0,2*PI);
		setOffsets(rg*(sin(theta0)+sin(theta1)), rg*(sin(theta0)+sin(theta1)), rg*(cos(theta2)-cos(theta1)), rg*(sin(theta2)-sin(theta1)));
	}
	//double correrr = sim_rnd_source.Gaus(0.0,1.0);
	for(unsigned int t=0; t<nBetaTubes; t++) {
		if(currentCal->scaleNoiseWithL)
			tubeRes = pmtRes[mySide][t]*en;
		else
			tubeRes = currentCal->nPE(mySide,t,en,x+dsx,y+dsy,0);
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
			sevt.adc[t] = currentCal->invertCorrections(mySide,t,sevt.tuben[t].x,x+dsx,y+dsy,0); 
			float pedw = currentCal->getPedwidth(currentCal->sensorNames[mySide][t],0.0);
			sevt.adc[t] += sim_rnd_source.Gaus(0.0,pedw);
			//sevt.adc[t] += sim_rnd_source.Gaus(0.0,25);
			//sevt.adc[t] += (correrr*0.7+sim_rnd_source.Gaus(0.0,0.7))*20.0;
		}
	}
	if(calcADC)
		currentCal->calibrateEnergy(mySide, x+dwx, y+dwy, sevt, 0);
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

void PMTGenerator::setPosition(float xx, float yy) {
	x = xx;
	y = yy;
	for(Side s = EAST; s <= WEST; ++s) {
		resTot[s] = 0;
		for(unsigned int t=0; t<nBetaTubes; t++) {
			pmtRes[s][t] = currentCal->nPE(s, t, 300.0, x, y, 0)/300.0;
			resTot[s] += pmtRes[s][t];
		}
	}
}
