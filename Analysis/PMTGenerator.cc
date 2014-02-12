#include "PMTGenerator.hh"
#include "QFile.hh"
#include "strutils.hh"
#include "PathUtils.hh"
#include <algorithm>
#include <Math/DistFunc.h>

TRandom3 PMTGenerator::sim_rnd_source;


double TriggerProb::calcProb() {
	float p0 = tubeProbs[0];
	float p1 = tubeProbs[1];
	float p2 = tubeProbs[2];
	float p3 = tubeProbs[3];
		
	// prob 2 of 3 p1*(1-(1-p2)*(1-p3))+(1-p1)*p2*p3; :
	float p2of3 = p1*p2 + p2*p3 + p3*p1 -2*p1*p2*p3;
	
	return p0*(1-(1-p1)*(1-p2)*(1-p3)) + (1-p0)*p2of3;
}

double TriggerProbMLP::unfold(double x) {
	if(x<1e-6) x = 1e-6;
	if(x>0.999999) x = 0.999999;
	return atan(tan((x-0.5)*M_PI)/5.);
}

void TriggerProbMLP::condition(Double_t* aIn) {
	unsigned int nn=nBetaTubes;
	for(unsigned int t1=1; t1<nBetaTubes; t1++)
		for(unsigned int t2=0; t2<t1; t2++)
			aIn[nn++] = unfold(aIn[t1]*aIn[t2]);
	for(unsigned int t=4; t<nBetaTubes; t++) aIn[t] = unfold(aIn[t]);
}

double TriggerProbMLP::calcProb() {
	double p0 = TriggerProb::calcProb();
	if(p0<0.005 || p0 > 0.99) return  p0;
	condition(tubeProbs);
	return TMLP->Evaluate(0,tubeProbs)+0.5;
}
	
	
	
// note, crosstalk = 0.010 was used for 2010 data analysis.

PMTGenerator::PMTGenerator(Side s, float xx, float yy):
x(xx), y(yy), xw(xx), yw(yy), evtm(0), presmear(0), dgain(16.0), pedcorr(0.2), crosstalk(0.010), xscatter(0.), trigThreshScale(1.0),
currentCal(NULL), TProb(new TriggerProb()), mySide(s) {
	for(Side s = EAST; s <= WEST; ++s)
		for(unsigned int t=0; t<nBetaTubes; t++)
			lightBal[s][t] = 1.;
}

void PMTGenerator::setCalibrator(PMTCalibrator* P) { 
	assert(P);
	currentCal = P;
	setPosition(x,y,xw-x,yw-y);
}

void PMTGenerator::setTriggerProb(TriggerProb* TP) {
	assert(TP);
	TProb = TP;
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

void PMTGenerator::setLightbal(Side s, float l1, float l2, float l3, float l4) {
	assert(s==EAST || s==WEST);
	lightBal[s][0] = l1;
	lightBal[s][1] = l2;
	lightBal[s][2] = l3;
	lightBal[s][3] = l4;
}

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

// input array of uncorrelated normalized (mu=0, sigma=1) random variables and desired correlation
// re-correlates variables, applies optional final scaling
void makeCorrelated(float* v, unsigned int n, float c, float* sigma = NULL) {
	// calculate recorrelation matrix parameters
	float a = sqrt(1+(n-1)*(n-1)-(n-1)*(n-2)*c + 2*(n-1)*sqrt((1-c)*(1+(n-1)*c)))/n;
	float b = ( 2*c*(n-1) + (n-2)*(sqrt((1-c)*(1+(n-1)*c))-1) ) / (n*n*a);
	
	// apply recorrelation matrix
	std::vector<float> v2(n);
	for(unsigned int i=0; i<n; i++)
		for(unsigned int j=0; j<n; j++)
			v2[i] += v[j]*(i==j?a:b);
	for(unsigned int i=0; i<n; i++)
		v[i]=v2[i]*(sigma?sigma[i]:1);
}

// makes correlation "waves" where some events are 4-way correlated, others uncorrelated
void makeCorrelated4(float* v, unsigned int n, float c, float* sigma = NULL) {
	std::vector<float> v2(n);
	if(PMTGenerator::sim_rnd_source.Uniform(0,1)<c) {
		float r = PMTGenerator::sim_rnd_source.Gaus(0.,1.);
		for(unsigned int i=0; i<n; i++) v2[i] = r;
	} else {
		for(unsigned int i=0; i<n; i++) v2[i] = v[i];
	}
	for(unsigned int i=0; i<n; i++)
		v[i]=v2[i]*(sigma?sigma[i]:1);
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
		nPE[t] = tubeRes[t]*lightBal[mySide][t];
	}
	preuncorrelate(nPE, crosstalk);
	for(unsigned int t=0; t<nBetaTubes; t++) {
		nPE[t] = sim_rnd_source.PoissonD(nPE[t]>0?nPE[t]:0);				//< primary photoelectrons
		nPE[t] = sim_rnd_source.PoissonD(dgain*(nPE[t]>0?nPE[t]:0))/dgain;	//< first gain stage electron multiplication
	}
	recorrelate(nPE, crosstalk);
	for(unsigned int t=0; t<nBetaTubes; t++)
		sevt.tuben[t] = nPE[t]/tubeRes[t]*en;
	
	
	// correlated pedestal noise
	float pedw[nBetaTubes]; // pedestal widths
	float dped[nBetaTubes]; // noise term on each pedestal
	for(unsigned int t=0; t<nBetaTubes; t++) {
		pedw[t] = currentCal->getPedwidth(currentCal->sensorNames[mySide][t], evtm);
		dped[t] = sim_rnd_source.Gaus(0.,1.);
	}
	makeCorrelated(dped,nBetaTubes,pedcorr,pedw);
	
	for(unsigned int t=0; t<nBetaTubes; t++) {
		// adc value
		sevt.adc[t] = currentCal->invertCorrections(mySide, t, sevt.tuben[t].x, x, y, evtm);
		// correlated pedestal noise
		sevt.adc[t] += dped[t];
		// quantization noise
		float pedx = currentCal->getPedestal(currentCal->sensorNames[mySide][t], evtm);
		sevt.adc[t] = int(sevt.adc[t]+pedx+0.5)-pedx;
	}
	currentCal->calibrateEnergy(mySide, xw, yw, sevt, evtm, xscatter);
	return sevt;
}

unsigned int PMTGenerator::triggers() {
	nTrigs = 0;
	unsigned int nZero = 0;
	for(unsigned int t=0; t<nBetaTubes; t++) {
		TProb->tubeProbs[t] = currentCal->trigEff(mySide,t,sevt.adc[t]*trigThreshScale);
		pmtTriggered[t] = sim_rnd_source.Uniform(0.0,1.0) < TProb->tubeProbs[t];
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

bool PMTGenerator::triggered() {
	triggers();
	return sim_rnd_source.Uniform(0.0,1.0) < TProb->calcProb();
}

//-------------------------------------------------------------

void PMTGenerator::calcCathodeSignals(Side s, AxisDirection d, const float* cath_chg, float* cath_adc, wireHit& w) const {
	// generate simulated cathode ADCs
	for(unsigned int c=0; c<currentCal->nWires(s,d); c++) {
		// base ADC conversion
		cath_adc[c] = currentCal->ccloud_eta[s]->eval(s,0,xw,yw,true) * currentCal->cathseg_energy_norm[s][d][c] * cath_chg[c];
		// pedestals and electronics noise
		double ped = currentCal->cathPeds0[s][d][c];
		cath_adc[c] += ped + sim_rnd_source.Gaus(0.,currentCal->cathPedW0[s][d][c]);
		// clamp to ADC range and quantize
		cath_adc[c] = int(cath_adc[c]<0 ? 0 : (cath_adc[c]>4095 ? 4095:cath_adc[c]));
		// pedestal subtract
		cath_adc[c] -= ped;
	}
	
	// calculate wireplane info from simulated cathodes, but keep original center location
	double old_center = w.center;
	double old_rawcenter = w.rawCenter;
	w = currentCal->calcHitPos(s,d,cath_adc,&currentCal->cathPeds0[s][d][0]);
	w.center = old_center;
	w.rawCenter = old_rawcenter;
}

