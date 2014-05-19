#include "GainStabilizer.hh"
#include "EnergyCalibrator.hh"

float GainStabilizer::gmsFactor(Side, unsigned int, float) const { return 1.0; }

Stringmap GainStabilizer::gmsSummary() const {
	Stringmap m;
	for(Side s = EAST; s <= WEST; ++s) {
		for(unsigned int t=0; t<nBetaTubes; t++) {
			std::string tname = sideSubst("%c",s)+itos(t);
			m.insert(tname+"_gms",gmsFactor(s,t,0));
		}
	}
	return m;
}

ChrisGainStabilizer::ChrisGainStabilizer(RunNum myRun, CalDB* cdb, LinearityCorrector* myCorrecter):
GainStabilizer(myRun, cdb, myCorrecter) {
	for(Side s = EAST; s <= WEST; ++s) {
		for(unsigned int t=0; t<nBetaTubes; t++) {
			try {
				pulserPeak[s][t] = CDB->getRunMonitor(rn,LCor->sensorNames[s][t],"Chris_peak");
				pulser0[s][t] = CDB->getRunMonitorStart(LCor->rGMS,LCor->sensorNames[s][t],"Chris_peak");
			} catch (const Stringmap& e) {
				pulserPeak[s][t] = NULL;
				pulser0[s][t] = 0;
				e.display();
			}
			if(!pulser0[s][t] || !pulserPeak[s][t])
				printf("*** Missing Chris Pulser data to calibrate %i%c%i! ***\n",rn,sideNames(s),t);
		}
	}
}

float ChrisGainStabilizer::gmsFactor(Side s, unsigned int t, float time) const {
	if(pulser0[s][t] && pulserPeak[s][t]
	   && pulser0[s][t]>800 && pulserPeak[s][t]->Eval(time)>800)
		return pulser0[s][t]/pulserPeak[s][t]->Eval(time);
	else
		return 1.0;
}

Stringmap ChrisGainStabilizer::gmsSummary() const {
	Stringmap m = GainStabilizer::gmsSummary();
	for(Side s = EAST; s <= WEST; ++s) {
		for(unsigned int t=0; t<nBetaTubes; t++) {
			std::string tname = sideSubst("%c",s)+itos(t);
			if(!pulserPeak[s][t]) continue;
			m.insert(tname+"_nBiPts",itos(pulserPeak[s][t]->GetN()));
			m.insert(tname+"_BiMean",pulserPeak[s][t]->GetMean(2));
			m.insert(tname+"_BiRMS",pulserPeak[s][t]->GetRMS(2));
			m.insert(tname+"_Pulser0",pulser0[s][t]);
		}
	}
	return m;
}

void ChrisGainStabilizer::printSummary() {
	for(Side s = EAST; s <= WEST; ++s) {
			
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
			printf("\t%.3f",gmsFactor(s,t,0));
		printf("\t  Chris Pulser correction since reference run\n");
	}
}

TweakedGainStabilizer::TweakedGainStabilizer(GainStabilizer* BG): GainStabilizer(BG->rn,BG->CDB,BG->LCor), baseGain(BG) {
	for(Side s=EAST; s<=WEST; ++s)
		for(unsigned int t=0; t<=nBetaTubes; t++)
			CDB->getGainTweak(rn,s,t,eOrig[s][t],eFinal[s][t]);
}
float TweakedGainStabilizer::gmsFactor(Side s, unsigned int t, float time) const {
	smassert(s<=WEST && t<=nBetaTubes);
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
