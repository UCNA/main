#include "AsymmetryErrors.hh"
#include "PathUtils.hh"
#include <TRandom3.h>

void AsymData::write(QFile& qf, std::string pfx) const {
	qf.insert(pfx+"_SR_Asym",AsymSR->toStringmap());
	qf.insert(pfx+"_Bonehead_Asym",AsymBoner->toStringmap());
	Stringmap m;
	for(Side s = EAST; s<=WEST; ++s) {
		for(unsigned int t=0; t<=nBetaTubes; t++) {
			for(int afp = 0; afp < 2; afp++) {
				m.insert(std::string("ep_")+ctos(sideNames(s))+itos(t)+"_"+itos(afp),kurieEp[s][afp][t].x);
				m.insert(std::string("dep_")+ctos(sideNames(s))+itos(t)+"_"+itos(afp),kurieEp[s][afp][t].err);
			}
			m.insert(std::string("ep_")+ctos(sideNames(s))+itos(t),kurieEpAvg[s][t].x);
			m.insert(std::string("dep_")+ctos(sideNames(s))+itos(t),kurieEpAvg[s][t].err);
		}
	}
	qf.insert(pfx+"_KurieFits",m);
}

AsymErrorExplorer::AsymErrorExplorer(RunGeometry g, RunNum rn): CDB("mpm_debug") {
	
	std::string basePath = "../PostPlots/SimNonlinearity/";
	makePath(basePath);
	OM = new OutputManager("SimNonlinearities","../PostPlots/SimNonlinearity/");
	
	nBetas = 16000000;
	refMultiple = 4.0;
	flipon_effic = 0.66;
	nBinsE = 40;
	eMax = 800;
	fitMin = 225;
	fitMax = 675;
	
	// get a calibrator for run
	PMTCalibrator* PCal = new PMTCalibrator(rn,&CDB);
	
	// set up MC generator
	PGen.setCalibrator(PCal);
	PGen.calcADC = false;
	
	// load Jianglai's montecarlo beta spectra for each side/flipper state
	/*
	hInput[EAST][0] = getSimBetaSpectrum(EAST, g, AFP_OFF);
	hInput[EAST][1] = getSimBetaSpectrum(EAST, g, AFP_ON);
	hInput[WEST][0] = getSimBetaSpectrum(WEST, g, AFP_OFF);
	hInput[WEST][1] = getSimBetaSpectrum(WEST, g, AFP_ON);
	 */
	assert(false);
	n0 = hInput[EAST][0]->Integral();
	
	// simulated spectra for perfect linearity
	simulateSpinstates(OM,"Lin",NULL,int(nBetas*refMultiple),refSpectra);
	
	OM->write();
}


void AsymErrorExplorer::simulateSpinstates(OutputManager* OMz, std::string pfx, SimNonlinearizer* SNL, unsigned int nToSim, AsymData& ad) {
	
	for(Side s = EAST; s <= WEST; ++s)
		for(unsigned int t=0; t<=nBetaTubes; t++)
			ad.hSpinavg[s][t] = OMz->registeredTH1F(pfx+"_SpinAvg_"+ctos(sideNames(s))+itos(t),"Spin-Averaged Spectra",nBinsE,0,eMax);
	
	
	// simulate spectra
	AFPState afpStates[] = {AFP_OFF,AFP_ON};
	PGen.setNonlinearity(SNL);
	for(unsigned int afp = 0; afp < 2; afp++) {
		if(SNL)
			SNL->startNewRun(afpStates[afp]);
		for(Side s = EAST; s <= WEST; ++s) {
			PGen.setSide(s);
			for(unsigned int t=0; t<=nBetaTubes; t++) {
				ad.hSimulated[s][afp][t] = OMz->registeredTH1F(pfx+"_"+ctos(sideNames(s))+itos(afp)+itos(t),"Asymmetry Spectra",nBinsE,0,eMax);
				//PGen.autofillTuben[t]=ad.hSimulated[s][afp][t];
			}
			SectorCutter SC(7,50);
			//PGen.simSpectrumIntegrated(hInput[s][afp],SC,int(nToSim*(afp?flipon_effic:1.0)*hInput[s][afp]->Integral()/n0));
			//PGen.setPosition(0,0);
			//PGen.simSpectrum(hInput[s][afp],int(nToSim*(afp?flipon_effic:1.0)*hInput[s][afp]->Integral()/n0));
			
			for(unsigned int t=0; t<=nBetaTubes; t++) {
				ad.hSpinavg[s][t]->Add(ad.hSimulated[s][afp][t]);
				ad.kurieEp[s][afp][t] = kurieIterator(ad.hSimulated[s][afp][t],neutronBetaEp);
			}
		}
	}
	for(Side s = EAST; s <= WEST; ++s)
		for(unsigned int t=0; t<=nBetaTubes; t++)
			ad.kurieEpAvg[s][t] = kurieIterator(ad.hSpinavg[s][t],neutronBetaEp);
	
	// SR asymmetries
	ad.AsymSR = new SRAsym(ad.hSimulated[EAST][0][nBetaTubes],ad.hSimulated[WEST][0][nBetaTubes],
						   ad.hSimulated[EAST][1][nBetaTubes],ad.hSimulated[WEST][1][nBetaTubes],false);
	ad.AsymSR->fitAsym(fitMin,fitMax);
	ad.AsymSR->toStringmap().display("\t");
	// bonehead asymmetries
	ad.AsymBoner = new SRAsym(ad.hSimulated[EAST][0][nBetaTubes],ad.hSimulated[WEST][0][nBetaTubes],
							  ad.hSimulated[EAST][1][nBetaTubes],ad.hSimulated[WEST][1][nBetaTubes],true);
	ad.AsymBoner->fitAsym(fitMin,fitMax);
	ad.AsymBoner->toStringmap().display("\t");
	
	ad.write(OMz->qOut,pfx);
	if(SNL)
		OMz->qOut.insert(pfx+"_Nonlinearizer",SNL->toStringmap());
	if(refSpectra.AsymSR)
		refSpectra.AsymSR->drawSpectra(1);
	ad.AsymSR->drawSpectra(0,true);
	OMz->printCanvas(std::string("BetaSpectra_")+pfx);
}

void AsymErrorExplorer::gainmatchSNL(SimNonlinearizer& SNL, const AsymData& obs) {
	for(Side s = EAST; s<=WEST; ++s)
		for(unsigned int t=0; t<nBetaTubes; t++)
			SNL.gainfactor[s][t] *= (refSpectra.kurieEp[s][0][t].x+refSpectra.kurieEp[s][1][t].x)/(obs.kurieEp[s][0][t].x+obs.kurieEp[s][1][t].x);
}

AsymData AsymErrorExplorer::processTrial(OutputManager& OMz, SimNonlinearizer& SNL, std::string pfx) {
	
	// do simulation
	AsymData ad;
	simulateSpinstates(&OMz,pfx,&SNL,nBetas,ad);
	refSpectra.write(OMz.qOut,"Ref");
	
	// make plots
	refSpectra.AsymSR->hAnorm->SetMinimum(0);
	refSpectra.AsymSR->hAnorm->SetMaximum(0.07);
	refSpectra.AsymSR->hAnorm->SetLineColor(3);
	refSpectra.AsymSR->hAnorm->Draw();
	ad.AsymSR->hAnorm->SetLineColor(5);
	ad.AsymSR->hAnorm->Draw("Same");
	refSpectra.AsymSR->hAsym->SetLineColor(2);
	refSpectra.AsymSR->hAsym->Draw("Same");
	ad.AsymSR->hAsym->SetLineColor(4);
	ad.AsymSR->hAsym->Draw("Same");
	OMz.printCanvas(std::string("AsymmetrySR_")+pfx);
	
	refSpectra.AsymBoner->hAsym->SetMinimum(0);
	refSpectra.AsymBoner->hAsym->SetMaximum(0.07);
	refSpectra.AsymBoner->hAsym->SetLineColor(2);
	refSpectra.AsymBoner->hAsym->Draw();
	ad.AsymBoner->hAsym->SetLineColor(4);
	ad.AsymBoner->hAsym->Draw("Same");
	OMz.printCanvas(std::string("AsymmetryBoner_")+pfx);
	
	// ratio plot, fit
	TH1F* hRatio = OMz.registeredTH1F(pfx+"_Ratio","Asymmetry Comparison",nBinsE,0,eMax);
	hRatio->Sumw2();
	hRatio->Add(ad.AsymSR->hAsym);
	hRatio->Divide(refSpectra.AsymSR->hAsym);
	TF1 rfit("AsymFit","pol0",fitMin,fitMax);
	hRatio->Fit(&rfit,"QR");
	Stringmap m;
	m.insert("RatioFit",rfit.GetParameter(0));
	m.insert("dRatioFit",rfit.GetParError(0));
	m.insert("fitMin",rfit.GetXmin());
	m.insert("fitMax",rfit.GetXmax());
	OMz.qOut.insert(pfx+"_AsymRatioFit",m);
	hRatio->SetMinimum(0.5);
	hRatio->SetMaximum(1.5);
	hRatio->Draw();
	OMz.printCanvas(std::string("AsymmetryRatio_")+pfx);
	
	OMz.write();
	
	return ad;
}


void AsymmetryErrors() {
	
	int n;
	AsymData ad;
	RunGeometry geoms[] = {GEOMETRY_C,GEOMETRY_B,GEOMETRY_A,GEOMETRY_OTHER};
	
	for(RunGeometry* g = &geoms[0]; g < &geoms[1]; g++) {
		
		AsymErrorExplorer AEE(*g);
	
		// sequential scans
		if(0) {
			
			OutputManager OM_ShiftSequenceCorr(std::string("ShiftSequenceCorr_")+geomName(*g),AEE.OM);
			OutputManager OM_ShiftSequenceAnticorr(std::string("ShiftSequenceAnticorr_")+geomName(*g),AEE.OM);
			OutputManager OM_GainSequenceCorr(std::string("GainSequenceCorr_")+geomName(*g),AEE.OM);
			OutputManager OM_GainSequenceAnticorr(std::string("GainSequenceAnticorr_")+geomName(*g),AEE.OM);
			OutputManager OM_AFPGainSequenceCorr(std::string("AFPGainSequenceCorr_")+geomName(*g),AEE.OM);
			OutputManager OM_AFPGainSequenceAnticorr(std::string("AFPGainSequenceAnticorr_")+geomName(*g),AEE.OM);
			OutputManager OM_ResSequenceCorr(std::string("ResSequenceCorr_")+geomName(*g),AEE.OM);
			OutputManager OM_ResSequenceAnticorr(std::string("ResSequenceAnticorr_")+geomName(*g),AEE.OM);
			int nscan = 20;

			for(n=0; n<2*nscan+1; n++) {
				
				float a = float(n-nscan)/float(nscan);
				
				// energy shifts
				if(0) {
					OutputManager* OM1 = new OutputManager(std::string("Trial_")+itos(n),&OM_ShiftSequenceCorr);
					SimNonlinearizer SNL1;
					SNL1.setOffset(20*a,EAST);
					SNL1.setOffset(20*a,WEST);
					ad = AEE.processTrial(*OM1,SNL1,"Nonlin");
					delete(OM1);		
				}
				if(0) {
					OutputManager* OM1 = new OutputManager(std::string("Trial_")+itos(n),&OM_ShiftSequenceAnticorr);
					SimNonlinearizer SNL1;
					SNL1.setOffset(20*a,EAST);
					SNL1.setOffset(-20*a,WEST);
					ad = AEE.processTrial(*OM1,SNL1,"Nonlin");
					delete(OM1);		
				}
				
				// gain shifts
				if(0) {
					OutputManager* OM1 = new OutputManager(std::string("Trial_")+itos(n),&OM_GainSequenceCorr);
					SimNonlinearizer SNL1;
					SNL1.setRelerr(0.2*a,EAST);
					SNL1.setRelerr(0.2*a,WEST);
					ad = AEE.processTrial(*OM1,SNL1,"Nonlin");
					delete(OM1);		
				}
				if(0) {
					OutputManager* OM1 = new OutputManager(std::string("Trial_")+itos(n),&OM_GainSequenceAnticorr);
					SimNonlinearizer SNL1;
					SNL1.setRelerr(0.2*a,EAST);
					SNL1.setRelerr(-0.2*a,WEST);
					ad = AEE.processTrial(*OM1,SNL1,"Nonlin");
					delete(OM1);		
				}				
							
				// AFP-correlated gain shifts
				if(0) {
					OutputManager* OM1 = new OutputManager(std::string("Trial_")+itos(n),&OM_AFPGainSequenceCorr);
					SimNonlinearizer SNL1;
					SNL1.afpCorrGain = 1.0+0.05*a;
					ad = AEE.processTrial(*OM1,SNL1,"Nonlin");
					delete(OM1);		
				}
				if(0) {
					OutputManager* OM1 = new OutputManager(std::string("Trial_")+itos(n),&OM_AFPGainSequenceAnticorr);
					SimNonlinearizer SNL1;
					SNL1.afpCorrGain = -(1.0+0.05*a);
					ad = AEE.processTrial(*OM1,SNL1,"Nonlin");
					delete(OM1);
				}				
				
				// Energy resolution errors
				if(1) {
					OutputManager* OM1 = new OutputManager(std::string("Trial_")+itos(n),&OM_ResSequenceCorr);
					SimNonlinearizer SNL1;
					SNL1.setReserr(1+0.5*a,EAST);
					SNL1.setReserr(1+0.5*a,WEST);
					ad = AEE.processTrial(*OM1,SNL1,"Nonlin");
					delete(OM1);		
				}
				if(1) {
					OutputManager* OM1 = new OutputManager(std::string("Trial_")+itos(n),&OM_ResSequenceAnticorr);
					SimNonlinearizer SNL1;
					SNL1.setReserr(1+0.5*a,EAST);
					SNL1.setReserr(1-0.5*a,WEST);
					ad = AEE.processTrial(*OM1,SNL1,"Nonlin");
					delete(OM1);		
				}				
				

				
			}
		}
			
		// uglycurves
		if(0) {
			OutputManager OM_UglyAll(std::string("UglyAll_")+geomName(*g),AEE.OM);
			for(n=0; n<200; n++) {
				SimNonlinearizer SNL1;
				SNL1.makeRanderr();
				OutputManager* OM1 = new OutputManager(std::string("Trial_")+itos(n),&OM_UglyAll);
				ad = AEE.processTrial(*OM1,SNL1,"Nonlin");
				delete(OM1);		
			}
		}
		
		// E-W correlated uglycurves
		if(1) {
			OutputManager OM_UglyAll(std::string("UglyAll2011Corr_")+geomName(*g),AEE.OM);
			TRandom3 rnd_trial_source;
			rnd_trial_source.SetSeed(0);
			for(n=0; n<20; n++) {
				unsigned int nrand = rnd_trial_source.Integer(1<<16);
				SimNonlinearizer SNL1;
				SNL1.makeRanderr(true);
				OutputManager* OM1 = new OutputManager(std::string("Trial_")+itos(nrand),&OM_UglyAll);
				ad = AEE.processTrial(*OM1,SNL1,"Nonlin");
				delete(OM1);	
			}
		}		
	}
	
}
