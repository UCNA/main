#include "EnergyCalibrator.hh"
#include "GainStabilizer.hh"
#include "GraphUtils.hh"
#include "SQL_Utils.hh"
#include "SMExcept.hh"
#include <utility>
#include <TRandom3.h>

TRandom3 pcal_random_src;

PedestalCorrector::~PedestalCorrector() {
	for(std::map<std::string,TGraph*>::iterator it = pedestals.begin(); it != pedestals.end(); it++)
		delete(it->second);
}
bool PedestalCorrector::checkPedestals(const std::string& sensorName) {
	std::map<std::string,TGraph*>::const_iterator it = pedestals.find(sensorName);
	std::map<std::string,TGraph*>::const_iterator itw = pedwidths.find(sensorName);
	if( it != pedestals.end() && itw != pedwidths.end() ) 
		return true;
	TGraph* tg = pCDB->getPedestals(myRun,sensorName);
	TGraph* tgw = pCDB->getPedwidths(myRun,sensorName);
	if(!tg || !tgw)
		return false;
	pedestals.insert(std::make_pair(sensorName,tg));
	pedwidths.insert(std::make_pair(sensorName,tgw));
	return true;	
}
float PedestalCorrector::getPedestal(const std::string& sensorName, float time) {
	std::map<std::string,TGraph*>::const_iterator it = pedestals.find(sensorName);
	if( it != pedestals.end() ) 
		return it->second->Eval(time);
	TGraph* tg = pCDB->getPedestals(myRun,sensorName);
	if(!tg) {
		SMExcept e("missingPed");
		e.insert("sensor",sensorName);
		throw(e);
	}
	pedestals.insert(std::make_pair(sensorName,tg));
	return getPedestal(sensorName,time);
}
float PedestalCorrector::getPedwidth(const std::string& sensorName, float time) {
	std::map<std::string,TGraph*>::const_iterator it = pedwidths.find(sensorName);
	if( it != pedwidths.end() ) 
		return it->second->Eval(time);
	TGraph* tg = pCDB->getPedwidths(myRun,sensorName);
	if(!tg) {
		SMExcept e("missingPed");
		e.insert("sensor",sensorName);
		throw(e);
	}
	pedwidths.insert(std::make_pair(sensorName,tg));
	return getPedwidth(sensorName,time);
}
void PedestalCorrector::insertPedestal(const std::string& sensorName, TGraph* g) {
	if(g->GetN()<2)
		throw(SMExcept("tooFewPoints"));
	pedestals.erase(sensorName);
	pedestals.insert(std::make_pair(sensorName,g));
}






LinearityCorrector::LinearityCorrector(RunNum myRun, CalDB* cdb):
scaleNoiseWithL(true), P(cdb->getPositioningCorrector(myRun)), GS(NULL), rn(myRun), LCRef(NULL), CDB(cdb) {
	
	for(Side s = EAST; s <= WEST; ++s)
		for(unsigned int t=0; t<nBetaTubes; t++)
			sensorNames[s][t] = sideSubst("ADC%c",s)+itos(pmtHardwareNum(s,t))+"Beta";
	
	if(rn < 5000) {
		printf("******* Bogus Linearity Calibration Run! ********\n");
		return;
	}
	
	P->setNormAvg();
	
	rGMS = CDB->getGMSRun(myRun);
	if(isRefRun())
		LCRef = this;
	else
		LCRef = getCachedRun(rGMS,CDB);
	
	for(Side s = EAST; s<=WEST; ++s) {
		
		for(unsigned int t=0; t<nBetaTubes; t++) {
			linearityFunctions[s][t] = CDB->getLinearity(myRun,s,t);
			linearityInverses[s][t] = invertGraph(linearityFunctions[s][t]);
		}
		
		if(isRefRun()) {
			for(unsigned int t=0; t<nBetaTubes; t++) {
				if(linearityInverses[s][t]->GetN())
					expected_adc[s][t] = linearityInverses[s][t]->Eval(CDB->getEcalEvis(myRun,s,t)*eta(s,t,CDB->getEcalX(myRun,s),CDB->getEcalY(myRun,s)));
				gms0[s][t] = expected_adc[s][t]/CDB->getEcalADC(myRun,s,t);
				deltaL[s][t] = ( CDB->getNoiseWidth(myRun,s,t) * dLinearity(s,t,CDB->getNoiseADC(myRun,s,t),0.0) / 
								sqrt(linearityCorrector(s,t,CDB->getNoiseADC(myRun,s,t),0.0)) );
				deltaADC[s][t] = CDB->getNoiseWidth(myRun,s,t)/sqrt(CDB->getNoiseADC(myRun,s,t));
			}
		} else {
			for(unsigned int t=0; t<nBetaTubes; t++) {
				deltaL[s][t] = LCRef->getDeltaL(s,t);
				deltaADC[s][t] = LCRef->getDeltaL(s,t);
				gms0[s][t] = LCRef->gmsFactor(s,t,0.0);
			}
		}
	}
}

LinearityCorrector::~LinearityCorrector() {
	if(rn<5000)
		return;
	for(Side s = EAST; s<=WEST; ++s) {
		for(unsigned int t=0; t<nBetaTubes; t++) {
			//TODO
			//if(linearityFunctions[s][t]) delete(linearityFunctions[s][t]); //TODO why does this crash?
			//if(ledPeaks[s][t]) delete(ledPeaks[s][t]);
		}
	}
}


float LinearityCorrector::linearityCorrector(Side s, unsigned int t, float adc, float time) const {
	assert(t<=nBetaTubes);
	if(t<nBetaTubes) {
		if(linearityFunctions[s][t]->GetN())
			return linearityFunctions[s][t]->Eval(adc*gmsFactor(s,t,time));
		else {
			assert(IGNORE_DEAD_DB);
			return 0;
		}
	}
	return 0; //< TODO ref linearity PMT?
}
float LinearityCorrector::invertLinearityStabilized(Side s, unsigned int t, float l) const {
	assert(s<=WEST && t<nBetaTubes && linearityInverses[s][t]);
	return linearityInverses[s][t]->Eval(l);
}
float LinearityCorrector::invertLinearity(Side s, unsigned int t, float l, float time) const {
	return invertLinearityStabilized(s,t,l)/gmsFactor(s,t,time);
}
float LinearityCorrector::dLinearity(Side s, unsigned int t, float adc, float time) const {
	const float h = 0.1;
	return (linearityCorrector(s,t,adc+h,time)-linearityCorrector(s,t,adc-h,time))/(2.0*h);
}
float LinearityCorrector::dInverse(Side s, unsigned int t, float l, float time) const {
	const float h = 0.1;
	return (invertLinearity(s, t, l+h, time)-invertLinearity(s, t, l-h, time))/(2.0*h);
}
float_err LinearityCorrector::invertLinearity(Side s, unsigned int t, float_err l, float time) const {
	return float_err(invertLinearity(s,t,l.x,time),l.err*dInverse(s, t, l.x,time));
}

float LinearityCorrector::gmsFactor(Side s, unsigned int t, float time) const {
	if(GS)
		return GS->gmsFactor(s,t,time);
	return gms0[s][t];
}
LinearityCorrector* LinearityCorrector::getCachedRun(RunNum r,CalDB* cdb) {
	std::map<RunNum,LinearityCorrector*>::iterator it = cachedRuns.find(r);
	if(it!=cachedRuns.end())
		return it->second;
	LinearityCorrector* C = new LinearityCorrector(r,cdb);
	cachedRuns.insert(std::make_pair(r,C));
	return C;
}

std::map<RunNum,LinearityCorrector*> LinearityCorrector::cachedRuns = std::map<RunNum,LinearityCorrector*>();




//------------------------------------------------------------------------------------




PMTCalibrator::PMTCalibrator(RunNum rn, CalDB* cdb): LinearityCorrector(rn,cdb),
PedestalCorrector(rn,cdb), EvisConverter(rn,cdb), WirechamberCalibrator(rn,cdb) {
	if(myRun < 5000) {
		printf("******* Bogus PMT Calibration Run! ********\n");
		return;
	}
	printf("Creating PMTCalibrator for %i.\n",myRun);
	GS = new TweakedGainStabilizer(new ChrisGainStabilizer(myRun,CDB,this));
	for(Side s = EAST; s <= WEST; ++s) {
		for(unsigned int t=0; t<nBetaTubes; t++) {
			pmtEffic[s][t] = CDB->getTrigeff(myRun,s,t);
			clipThreshold[s][t] = 4000;
			if(checkPedestals(sensorNames[s][t]))
				clipThreshold[s][t] -= 1.1*getPedestal(sensorNames[s][t],0);
			else
				clipThreshold[s][t] -= 500;
		}
	}
	printSummary();
}
PMTCalibrator::~PMTCalibrator() {
	printf("Deleting PMTCalibrator for %i.\n",myRun);
	for(Side s = EAST; s <= WEST; ++s)
		for(unsigned int t=0; t<nBetaTubes; t++)
			if(pmtEffic[s][t])
				delete(pmtEffic[s][t]);
}

float PMTCalibrator::lightResolution(Side s, unsigned int t, float l, float time) const {
	if(scaleNoiseWithL)
		return sqrt(l)*deltaL[s][t];
	float adc = invertLinearity(s, t, l, time);
	if(adc<=0)
		return 0;
	return sqrt(adc)*deltaADC[s][t]*dLinearity(s, t, adc, time);
}
float PMTCalibrator::adcResolution(Side s, unsigned int t, float adc, float time) const {
	float l = linearityCorrector(s, t, adc, time);
	return lightResolution(s,t,l,time)*dInverse(s,t,l,time);
}
float PMTCalibrator::energyResolution(Side s, unsigned int t, float e0, float x, float y,float time) const {
	if(t==nBetaTubes)
		return combinedResolution(s,e0,x,y,time);
	return lightResolution(s,t,e0*eta(s,t,x,y),time)/eta(s,t,x,y);
}
float PMTCalibrator::combinedResolution(Side s, float e0, float x, float y,float time) const {
	float de = 0;
	for(unsigned int t=0; t<nBetaTubes; t++)
		de += 1.0/pow((double)energyResolution(s,t,e0,x,y,time),(double)2.0);
	return sqrt(1.0/de);
}
float PMTCalibrator::nPE(Side s, unsigned int t, float e0, float x, float y,float time) const {
	if(e0<=0)
		return 0;
	if(t==nBetaTubes) {
		float nsum = 0;
		for(unsigned int i=0; i<nBetaTubes; i++)
			nsum += nPE(s,i,e0,x,y,time);
		return nsum;
	}
	return pow(double(e0/energyResolution(s,t,e0,x,y,time)),2.0);
}
float PMTCalibrator::pmtSumPE(Side s, float e0, float x, float y, float time) const {
	if(e0<=0)
		return 0;
	float nsum = 0;
	for(unsigned int i=0; i<nBetaTubes; i++)
		nsum += 1.0/nPE(s,i,e0,x,y,time);
	return nBetaTubes*nBetaTubes/nsum;
}

float PMTCalibrator::nPE(Side s, unsigned int t, float adc, float time) const {
	float l0 = linearityCorrector(s,t,adc,time);
	if(l0<=0)
		return 0;
	return pow(double(l0/lightResolution(s,t,l0,time)),2.0);
}


float_err PMTCalibrator::calibratedEnergy(Side s, unsigned int t, float x, float y, float adc, float time) const {
	float l = linearityCorrector(s,t,adc,time); //< observed light
	float eta0 = eta(s,t,x,y); //< positioning factor
	float_err E;
	E.x = l/eta0;
	E.err = energyResolution(s,t,E.x,x,y,time);
	if(!(E.x==E.x))
		E = 0.;
	return E;
}

float PMTCalibrator::trigEff(Side s, unsigned int t, float adc) const {
	assert(s==EAST || s==WEST);
	assert(t<nBetaTubes);
	assert(pmtEffic[s][t]);
	return pmtEffic[s][t]->effic(adc);
}

void PMTCalibrator::pedSubtract(Side s, float* adc, float time) {
	for(unsigned int t=0; t<nBetaTubes; t++)
		adc[t] -= getPedestal(sensorNames[s][t],time);
}
void PMTCalibrator::calibrateEnergy(Side s, float x, float y, ScintEvent& evt, float time, float poserr) const {
	evt.energy.x = evt.energy.err = 0;
	float weight[nBetaTubes];
	// count clipped PMTs
	unsigned int nclipped = 0;
	for(unsigned int t=0; t<nBetaTubes; t++)
		nclipped += evt.adc[t]>clipThreshold[s][t]-300;
	
	for(unsigned int t=0; t<nBetaTubes; t++) {
		
		float eta0 = eta(s,t,x,y);
		if(poserr) eta0 *= pcal_random_src.Gaus(1.,poserr);
		float l0 = linearityCorrector(s,t,evt.adc[t],time); // tube observed light
		if(l0 != l0)
			l0 = 0;
		float E0 = l0/eta0; // tube observed energy keV
		
		if(evt.adc[t] < 5 || l0 < 50)
			weight[t] = eta0*50/pow(lightResolution(s,t,50,time),2.0); // = nPE/keV @ 50keV to avoid uncertainty at 0
		else
			weight[t] = eta0*l0/pow(lightResolution(s,t,l0,time),2.0); // = nPE/keV for this position
		
		// remove bad weights
		if(!(weight[t]>0.01 && weight[t]<10.0))
			weight[t] = 0;
		
		evt.nPE[t] = E0*weight[t];
		
		// de-weight for clipping, unless all PMTs clipped
		if(nclipped<nBetaTubes) {
			if(evt.adc[t]>clipThreshold[s][t]-300)
				weight[t] *= (clipThreshold[s][t]-evt.adc[t])/300.0;
			if(evt.adc[t]>clipThreshold[s][t])
				weight[t] = 0;
		}
		
		evt.energy.x += E0*weight[t];
		evt.energy.err += weight[t];
		evt.tuben[t].x = E0;
	}
	evt.energy.x /= evt.energy.err;
	evt.energy.err = sqrt(evt.energy.x/evt.energy.err);
	for(unsigned int t=0; t<nBetaTubes; t++)
		evt.tuben[t].err = sqrt(evt.energy.x/weight[t]);
	
	if(!(evt.energy.x==evt.energy.x)) evt.energy = 0.;	// NaN test
}

void PMTCalibrator::summedEnergy(Side s, float x, float y, ScintEvent& evt, float time) const {
	evt.energy.x = evt.energy.err = 0;
	for(unsigned int t=0; t<nBetaTubes; t++) {
		float weight = eta(s,t,x,y);
		float E0 = linearityCorrector(s,t,evt.adc[t],time)/weight;
		evt.energy.x += E0*weight;
		evt.energy.err += weight;
		evt.tuben[t].x = E0;
	}
	evt.energy.x /= evt.energy.err;
	evt.energy.err = 0;
	for(unsigned int t=0; t<nBetaTubes; t++)
		evt.tuben[t].err = 0;
}

float PMTCalibrator::invertCorrections(Side s, unsigned int t, float e0, float x, float y, float time) const {
	return invertLinearity(s,t,e0*eta(s,t,x,y),time);
}
float_err PMTCalibrator::invertCorrections(Side s, unsigned int t, float_err e0, float x, float y, float time) const {
	float_err adc;
	adc.x = invertCorrections(s,t,e0.x,x,y,time);
	adc.err = e0.err * dInverse(s, t, e0.x*eta(s,t,x,y), time) * eta(s,t,x,y);
	return adc;
}

void PMTCalibrator::printSummary() {
	if(myRun < 5000) {
		printf("-* Energy Calibrations %i BOGUS!!!! *-\n\n",myRun);
		return;
	}
	printf("-- Energy Calibrations %i (GMS to %i) -- [%s]\n",myRun,rGMS,CDB->getName().c_str());
	if(GS)
		GS->printSummary();
	for(Side s = EAST; s <= WEST; ++s) {
		printf("%c Res500:",sideNames(s));
		for(unsigned int t=0; t<nBetaTubes; t++)
			printf("\t%.1f",energyResolution(s,t,500.0,0,0,0));
		printf("\t  energy resolution [keV] at 500keV at center\n");
		printf("%c Adc100:",sideNames(s));
		for(unsigned int t=0; t<nBetaTubes; t++)
			printf("\t%.1f",calibratedEnergy(s,t,0,0,100.0,0).x);
		printf("\t  energy [keV] for ADC=100 at center\n");
		printf("%c Trig50:",sideNames(s));
		for(unsigned int t=0; t<nBetaTubes; t++) {
			if(pmtEffic[s][t])
				printf("\t%.1f",calibratedEnergy(s,t,0,0,pmtEffic[s][t]->getThreshold(),0).x);
			else
				printf("\t???");
		}
		printf("\t  keV at 50%% PMT trigger threshold\n");
		printf("%c TrigPE:",sideNames(s));
		for(unsigned int t=0; t<nBetaTubes; t++) {
			if(pmtEffic[s][t])
				printf("\t%.1f",nPE(s,t,pmtEffic[s][t]->getThreshold(),0));
			else
				printf("\t???");
		}		
		printf("\t  photoelectrons at 50%% PMT trigger threshold\n");
	}
	WirechamberCalibrator::printSummary();
	printf("----------------------------------------------\n\n");
}

Stringmap PMTCalibrator::calSummary() const { 
	Stringmap m = wirecalSummary();
	if(GS) m += GS->gmsSummary();
	m.insert("run",itos(myRun));
	m.insert("tstart",itos(CDB->startTime(myRun)));
	m.insert("tend",itos(CDB->endTime(myRun)));
	
	for(Side s = EAST; s <= WEST; ++s) {
		for(unsigned int t=0; t<nBetaTubes; t++) {
			std::string tname = ctos(sideNames(s))+itos(t);
			m.insert(tname+"_deltaL",getDeltaL(s,t));
			m.insert(tname+"_adc100",calibratedEnergy(s,t,0,0,100.0,0).x);
			m.insert(tname+"_res500",energyResolution(s,t,500.0,0,0,0));
			if(pmtEffic[s][t])
				m.insert(tname+"_trig50",calibratedEnergy(s,t,0,0,pmtEffic[s][t]->getThreshold(),0).x);
		}
	}
	
	return m;
}	
