#include "GainStabilizer.hh"
#include "EnergyCalibrator.hh"

float GainStabilizer::gmsFactor(Side s, unsigned int t, float time) const {
	return LCor->getGMS0(s,t);
}

Stringmap GainStabilizer::gmsSummary() const {
	Stringmap m;
	for(Side s = EAST; s <= WEST; ++s) {
		for(unsigned int t=0; t<nBetaTubes; t++) {
			std::string tname = sideSubst("%c",s)+itos(t);
			m.insert(tname+"_gms0",LCor->getGMS0(s, t));
			m.insert(tname+"_gms",gmsFactor(s,t,0));
		}
	}
	return m;
}

ChrisGainStabilizer::ChrisGainStabilizer(RunNum myRun, CalDB* cdb, LinearityCorrector* myCorrecter):
GainStabilizer(myRun, cdb, myCorrecter) {
	for(Side s = EAST; s <= WEST; ++s) {
		for(unsigned int t=0; t<nBetaTubes; t++) {
			pulserPeak[s][t] = CDB->getRunMonitor(rn,LCor->sensorNames[s][t],"Chris_peak");
			pulser0[s][t] = CDB->getRunMonitorStart(LCor->rGMS,LCor->sensorNames[s][t],"Chris_peak");
			if(!pulser0[s][t] || !pulserPeak[s][t])
				printf("*** Missing Chris Pulser data to calibrate %i%c%i! ***\n",rn,sideNames(s),t);
		}
	}
}

float ChrisGainStabilizer::gmsFactor(Side s, unsigned int t, float time) const {
	if(pulser0[s][t] && pulserPeak[s][t]
	   && pulser0[s][t]>800 && pulserPeak[s][t]->Eval(time)>800)
		return LCor->getGMS0(s,t)*pulser0[s][t]/pulserPeak[s][t]->Eval(time);
	else
		return LCor->getGMS0(s,t);
}

void ChrisGainStabilizer::printSummary() {
	for(Side s = EAST; s <= WEST; ++s) {
		
		printf("%c t0 GMS0:",sideNames(s));
		for(unsigned int t=0; t<nBetaTubes; t++)
			printf("\t%.3f",LCor->getGMS0(s,t));
		printf("\t  GMS correction at reference run\n");
		
		printf("%c pulser0:",sideNames(s));
		for(unsigned int t=0; t<nBetaTubes; t++) {
			if(pulser0[s][t])
				printf("\t%.1f",pulser0[s][t]);
			else
				printf("\t???");
		}
		printf("\t  Chris Pulser at reference run\n");
		printf("%c pulser1:",sideNames(s));
		for(unsigned int t=0; t<nBetaTubes; t++) {
			if(pulserPeak[s][t])
				printf("\t%.1f",pulserPeak[s][t]->Eval(0));
			else
				printf("\t???");
		}
		printf("\t  Chris Pulser at current run\n");
		printf("%cChrisGMS:",sideNames(s));
		for(unsigned int t=0; t<nBetaTubes; t++)
			printf("\t%.3f",gmsFactor(s,t,0)/LCor->getGMS0(s,t));
		printf("\t  Chris Pulser correction since reference run\n");
	}
}

TweakedGainStabilizer::TweakedGainStabilizer(GainStabilizer* BG): GainStabilizer(BG->rn,BG->CDB,BG->LCor), baseGain(BG) {
	for(Side s=EAST; s<=WEST; ++s)
		for(unsigned int t=0; t<=nBetaTubes; t++)
			CDB->getGainTweak(rn,s,t,eOrig[s][t],eFinal[s][t]);
}
float TweakedGainStabilizer::gmsFactor(Side s, unsigned int t, float time) const {
	assert(s<=WEST && t<=nBetaTubes);
	if(t==nBetaTubes) return baseGain->gmsFactor(s,t,time)*eFinal[s][t]/eOrig[s][t];
	return baseGain->gmsFactor(s,t,time)*LCor->invertLinearityStabilized(s,t,eFinal[s][t])/LCor->invertLinearityStabilized(s,t,eOrig[s][t]);
}
Stringmap TweakedGainStabilizer::gmsSummary() const {
	Stringmap m = baseGain->gmsSummary();
	for(Side s = EAST; s <= WEST; ++s) {
		for(unsigned int t=0; t<=nBetaTubes; t++) {
			std::string tname = sideSubst("%c",s)+itos(t);
			m.insert(tname+"_eOrig",eOrig[s][t]);
			m.insert(tname+"_eFinal",eFinal[s][t]);
		}
	}
	return m;
}
void TweakedGainStabilizer::printSummary() {
	baseGain->printSummary();
	for(Side s = EAST; s <= WEST; ++s) {
		printf("%c tweak e0:",sideNames(s));
		for(unsigned int t=0; t<=nBetaTubes; t++)
			printf("\t%.1f",eOrig[s][t]);
		printf("\t  GMS tweaking original energy\n");
		
		printf("%c tweak e1:",sideNames(s));
		for(unsigned int t=0; t<=nBetaTubes; t++)
			printf("\t%.1f",eFinal[s][t]);
		printf("\t  GMS tweaking final energy\n");	
	}
}

/*
 LEDCorrector::LEDCorrector(RunNum myRun, CalDB* cdb): CDB(cdb), rn(myRun) {
 if(myRun < 5000) {
 printf("******* Bogus LED Calibration Run! ********\n");
 return;
 }
 for(Side s = EAST; s<=WEST; ++s) {
 refLinearity[s] = CDB->getLinearity(myRun,s,4);
 assert(refLinearity[s]);
 refInverses[s] = invertGraph(refLinearity[s]);
 for(unsigned int p=0; p<2; p++)
 co60Peaks[s][p] = CDB->getCo60(myRun, s, p);
 }
 refLED[EAST]=refLED[WEST]=NULL;
 refLED[EAST] = CDB->getLED(myRun, "ADCRefELED");
 refLED[WEST] = CDB->getLED(myRun, "ADCRefWLED");
 if(!(refLED[EAST] && refLED[WEST]))
 printf("*** WARNING: missing reference LED signal, will fake it...\n");
 }
 LEDCorrector::~LEDCorrector() {
 if(rn<5000)
 return;
 for(Side s = EAST; s <= WEST; ++s) {
 //if(refLED[s]) delete(refLED[s]);
 //if(co60Peaks[s][0]) delete(co60Peaks[s][0]);
 //if(co60Peaks[s][0]) delete(co60Peaks[s][1]);
 //if(refLinearity[s]) delete(refLinearity[s]);
 }
 }
 float LEDCorrector::ledBrightness(Side s, float time) {
 if(rn<5000)
 return 1.0;
 if(!refLED[s])
 return 1.0;
 return refLinearity[s]->Eval(refLED[s]->Eval(time) * 2500/adc_Co(s,time))/2500;
 }
 float LEDCorrector::adc_Co(Side s, float time) { 
 if(rn<5000)
 return 1.0;
 if(!(co60Peaks[s][0] && co60Peaks[s][1])) {
 printf("*** WARNING: Co60 Peaks MISSING!\n");
 return 2500.0;
 }
 return 0.5*co60Peaks[s][0]->Eval(time)+0.5*co60Peaks[s][1]->Eval(time);
 }
 */

/*
 bool LinearityCorrector::setMuonGMS() {
 if(muonGMS) return true;
 ledGMS = false;
 kurieGMS = false;
 muonGMS = true;
 LinearityCorrector* LCRef = NULL;
 if(!isRefRun()) {
 LCRef= getCachedRun(rGMS,CDB);
 assert(LCRef->setMuonGMS());
 }
 for(Side s = EAST; s <= WEST; ++s) {
 for(unsigned int t=0; t<nBetaTubes; t++) {
 if(!muonPeaks[s][t])
 muonPeaks[s][t] = CDB->getMuonPeak(rn, sensorNames[s][t]);
 if(!muonPeaks[s][t]) {
 printf("*** Warning: missing Muon peaks %c%i; GMS disabled!\n",sideNames(s),t);
 muonGMS = false;
 break;
 }
 if(isRefRun())
 muon0[s][t] = getMuonPeak(s,t,0.0);
 else
 muon0[s][t] = LCRef->getMuonPeak(s,t,0.0);
 }		
 }
 return muonGMS;
 }
 bool LinearityCorrector::setLEDGMS() {
 if(ledGMS) return true;
 ledGMS = true;
 kurieGMS = false;
 muonGMS = false;
 LinearityCorrector* LCRef = NULL;
 if(!isRefRun()) {
 LCRef= getCachedRun(rGMS,CDB);
 assert(LCRef->setLEDGMS());
 }	
 for(Side s = EAST; s<=WEST; ++s) {
 for(unsigned int t=0; t<nBetaTubes; t++) {
 if(!ledPeaks[s][t])
 ledPeaks[s][t] = CDB->getLED(rn, sensorNames[s][t]);
 if(!ledPeaks[s][t]) {
 printf("*** Warning: missing LED peaks %c%i; GMS disabled!\n",sideNames(s),t);
 ledGMS = false;
 break;
 }
 if(isRefRun()) {
 if(linearityFunctions[s][t]->GetN())
 energy_peg[s][t] = linearityFunctions[s][t]->Eval(expected_adc[s][t]/CDB->getEcalADC(rn,s,t)*getLEDadc(s,t,0)) / ledBrightness(s,0);
 } else
 energy_peg[s][t] = LCRef->getEnergyPeg(s,t);
 }
 }
 return ledGMS;
 }
 bool LinearityCorrector::setKurieGMS() {
 if(kurieGMS) return true;
 kurieGMS = true;
 ledGMS = false;
 muonGMS = false;
 for(Side s = EAST; s <= WEST; ++s) {
 for(unsigned int t=0; t<nBetaTubes; t++) {
 float ken = CDB->getKurieEnergy(rn,s,t);
 float kadc = CDB->getKurieADC(rn,s,t);
 if(ken*kadc == 0) {
 printf("*** Warning: missing Kurie GMS data for Run %i %c%i; defaulting to Muon GMS!\n",rn,sideNames(s),t);
 kurieGMS = false;
 break;
 }
 if(linearityInverses[s][t]->GetN())
 kgms[s][t] = linearityInverses[s][t]->Eval(eta(s,t,0,0)*ken)/kadc;
 else
 assert(IGNORE_DEAD_DB);
 }
 }
 if(!kurieGMS)
 setMuonGMS();
 return kurieGMS;
 }
 bool LinearityCorrector::disableGMS() {
 muonGMS = ledGMS = kurieGMS = false;
 return true;
 }
 */

/*
 float LinearityCorrector::getLEDadc(Side s, unsigned int t, float time) {
 if(t<nBetaTubes) {
 if(ledPeaks[s][t])
 return ledPeaks[s][t]->Eval(time);
 return 1000;
 }
 if(refLED[s])
 return refLED[s]->Eval(time);
 return 1000;
 }
 */

/*
 if(muonGMS)
 return gms0[s][t]*muon0[s][t]/getMuonPeak(s, t, time);
 if(ledGMS)
 return linearityInverses[s][t]->Eval(ledBrightness(s,time)*energy_peg[s][t])/getLEDadc(s, t, time);
 if(kurieGMS)
 return kgms[s][t];
 */
