#include "ucnaDataAnalyzer11b.hh"
#include "strutils.hh"
#include "ManualInfo.hh"
#include "GraphicsUtils.hh"
#include "MultiGaus.hh"
#include "SMExcept.hh"
#include <stdio.h>
#include <unistd.h>
#include <TStyle.h>
#include <TDatime.h>

RangeCut::RangeCut(const Stringmap& m): start(m.getDefault("start",0.0)), end(m.getDefault("end",0.0)) {}

ucnaDataAnalyzer11b::ucnaDataAnalyzer11b(RunNum R, std::string bp, CalDB* CDB):
TChainScanner("h1"), OutputManager("spec_"+itos(R),bp+"/hists/"), analyzeLED(false), needsPeds(false),
rn(R), PCal(R,CDB), CDBout(NULL), fAbsTimeEnd(0), deltaT(0), totalTime(0), nLiveTrigs(0), ignore_beam_out(false),
nFailedEvnb(0), nFailedBkhf(0), gvMonChecker(5,5.0), prevPassedCuts(true), prevPassedGVRate(true) {
	plotPath = bp+"/figures/run_"+itos(R)+"/";
	dataPath = bp+"/data/";
	for(Side s = EAST; s <= WEST; ++s)
		for(AxisDirection d = X_DIRECTION; d <= Y_DIRECTION; ++d)
			fMWPC_caths[s][d].resize(kMaxCathodes);
}

void ucnaDataAnalyzer11b::analyze() {
	loadCuts();
	setupOutputTree();
	
	pedestalPrePass();
	printf("\nRun wall time is %.1fs\n\n",wallTime);
	setupHistograms();
	
	printf("Scanning input data...\n");
	startScan();
	nextPoint();	// load data for first event
	while(processEvent()) continue;
	printf("Done.\n");
	
	processBiPulser();
	muonVetoAccidentals();
	calcTrigEffic();
	tallyRunTime();
	locateSourcePositions();
	plotHistos();
	replaySummary();
	
	quickAnalyzerSummary();
}

void ucnaDataAnalyzer11b::loadCut(CutVariable& c, const std::string& cutName) {
	std::vector<Stringmap> v = ManualInfo::MI.getInRange(cutName,rn);
	if(!v.size()) {
		SMExcept e("missingCut");
		e.insert("cutName",cutName);
		e.insert("runNum",rn);
		throw(e);
	}
	c.R = RangeCut(v[0]);
	printf("Loaded cut %s/%i = (%g,%g)\n",cutName.c_str(),rn,c.R.start,c.R.end);
}

void ucnaDataAnalyzer11b::loadCuts() {
	for(Side s = EAST; s <= WEST; ++s) {
		loadCut(fMWPC_anode[s], sideSubst("Cut_MWPC_%c_Anode",s));
		loadCut(fCathMax[s],sideSubst("Cut_MWPC_%c_CathMax",s));
		loadCut(fCathSum[s],sideSubst("Cut_MWPC_%c_CathSum",s));
		loadCut(fCathMaxSum[s],sideSubst("Cut_MWPC_%c_CathMaxSum",s));
		loadCut(fBacking_tdc[s], sideSubst("Cut_TDC_Back_%c",s));
		loadCut(fDrift_tac[s], sideSubst("Cut_ADC_Drift_%c",s));
		loadCut(fScint_tdc[s][nBetaTubes], sideSubst("Cut_TDC_Scint_%c_Selftrig",s));
		ScintSelftrig[s] = fScint_tdc[s][nBetaTubes].R;
		loadCut(fScint_tdc[s][nBetaTubes], sideSubst("Cut_TDC_Scint_%c",s));
		for(unsigned int t=0; t<nBetaTubes; t++)
			loadCut(fScint_tdc[s][t], sideSubst("Cut_TDC_Scint_%c_",s)+itos(t));
	}
	loadCut(fTop_tdc[EAST], "Cut_TDC_Top_E");
	loadCut(fBeamclock,"Cut_BeamBurst");
	loadCut(fWindow,"Cut_ClusterEvt");
	if(ignore_beam_out)
		fBeamclock.R.end = FLT_MAX;
	manualCuts = ManualInfo::MI.getRanges(itos(rn)+"_timecut");
	if(manualCuts.size())
		printf("Manually cutting %i time ranges...\n",(int)manualCuts.size());
}

void ucnaDataAnalyzer11b::checkHeaderQuality() {
	fEvnbGood = fBkhfGood = true;
	for(size_t i=0; i<kNumModules; i++) {
		if(int(r_Evnb[i]-iTriggerNumber)) fEvnbGood = false;
		if(int(r_Bkhf[i])!=17) fBkhfGood = false;
	}
	nFailedEvnb += !fEvnbGood;
	nFailedBkhf += !fBkhfGood;
}

void ucnaDataAnalyzer11b::convertReadin() {
	iSis00 = int(r_Sis00);
	iTriggerNumber = int(r_TriggerNumber);
	for(Side s = EAST; s <= WEST; ++s) {
		for(unsigned int t=0; t<=nBetaTubes; t++)
			fScint_tdc[s][t].val = r_PMTTDC[s][t];
		for(AxisDirection d = X_DIRECTION; d <= Y_DIRECTION; ++d) 
			for(unsigned int c=0; c<kMaxCathodes; c++)
				fMWPC_caths[s][d][c] = r_MWPC_caths[s][d][c];
		fMWPC_anode[s].val = r_MWPC_anode[s];
		for(unsigned int t=0; t<nBetaTubes; t++)
			sevt[s].adc[t] = r_PMTADC[s][t];
	}
	for(unsigned int i=0; i<kNumUCNMons; i++)
		fMonADC[i].val = r_MonADC[i];
}

void ucnaDataAnalyzer11b::calibrateTimes() {
	
	// start/end times
	if(currentEvent==0) {
		fAbsTimeStart = r_AbsTime;
		prevPassedCuts = prevPassedGVRate = true;
		totalTime = deltaT = 0;
		nLiveTrigs = 0;
	}
	if(r_AbsTime>fAbsTimeEnd) fAbsTimeEnd = r_AbsTime;
	
	// convert microseconds to seconds
	fTimeScaler = 1.e-6 * r_Clk;
	fBeamclock.val = 1.e-6 * r_BClk;	
	fDelt0 = 1.e-6 * r_Delt0;
	
	// check for overflow condition
	if(fTimeScaler[BOTH] < totalTime[BOTH]-deltaT[BOTH]-1000.0) {
		printf("\tFixing timing scaler overflow... ");
		deltaT[BOTH] += pow(2,32)*1.e-6;
		for(Side s = EAST; s<=WEST; ++s)
			deltaT[s] = totalTime[s];
	}
	// add overflow wraparound time
	fTimeScaler += deltaT;
	
	if(isUCNMon(UCN_MON_GV))
		gvMonChecker.addCount(fTimeScaler[BOTH]);
	else
		gvMonChecker.moveTimeLimit(fTimeScaler[BOTH]);
	
	// check global time cuts, add new blips as necessary
	bool passedGVRate = (gvMonChecker.getCount() == gvMonChecker.nMax) || (prevPassedGVRate && gvMonChecker.getCount() > gvMonChecker.nMax/2.0);
	prevPassedGVRate = passedGVRate;
	fPassedGlobal = passesBeamCuts() && (ignore_beam_out || fTimeScaler[BOTH]<gvMonChecker.lMax || passedGVRate);
	if(fPassedGlobal != prevPassedCuts) {
		if(!fPassedGlobal) {
			Blip b;
			b.start = 0.5*(fTimeScaler+totalTime);
			cutBlips.push_back(b);
		} else {
			cutBlips.back().end = 0.5*(fTimeScaler+totalTime);
		}
	}
	
	prevPassedCuts = fPassedGlobal;
	totalTime = fTimeScaler;
	nLiveTrigs += fPassedGlobal;
}

unsigned int ucnaDataAnalyzer11b::nFiring(Side s) const {
	unsigned int nf = 0;
	for(unsigned int t=0; t<nBetaTubes; t++)
		nf += pmtFired(s,t);
	return nf;
}

bool ucnaDataAnalyzer11b::isPulserTrigger() {
	if(iSis00 & (1<<5) && !(trig2of4(EAST)||trig2of4(WEST)))
		return true;
	for(Side s = EAST; s <= WEST; ++s) {
		unsigned int nthresh = 0;
		unsigned int nhigh = 0;
		for(unsigned int t=0; t<nBetaTubes; t++) {
			nthresh += (sevt[s].adc[t] > 200);
			nhigh += (sevt[s].adc[t] > 1500);
		}
		if(nhigh == 1 && nthresh == 1)
			return true;
	}
	return false;
}


bool ucnaDataAnalyzer11b::passesBeamCuts() {
	// basic time-since-beam cut
	if(!fBeamclock.inRange())
		return false;	
	// remove manually tagged segments
	for(std::vector< std::pair<double,double> >::const_iterator it = manualCuts.begin(); it != manualCuts.end(); it++)
		if (it->first <= fTimeScaler[BOTH] && fTimeScaler[BOTH] <= it->second)
			return false;
	return true;
}

void ucnaDataAnalyzer11b::reconstructPosition() {
	for(Side s = EAST; s <= WEST; ++s) {
		for(AxisDirection d = X_DIRECTION; d <= Y_DIRECTION; ++d) {
			std::vector<float> cathPeds;
			for(unsigned int c=0; c<cathNames[s][d].size(); c++)
				cathPeds.push_back(PCal.getPedestal(cathNames[s][d][c],fTimeScaler[BOTH]));
			wirePos[s][d] = PCal.calcHitPos(s,d,fMWPC_caths[s][d],cathPeds);
		}
		fMWPC_anode[s].val -= PCal.getPedestal(sideSubst("MWPC%cAnode",s),fTimeScaler[BOTH]);
		fCathSum[s].val = wirePos[s][X_DIRECTION].cathodeSum + wirePos[s][Y_DIRECTION].cathodeSum;
		fCathMax[s].val = wirePos[s][X_DIRECTION].maxValue<wirePos[s][Y_DIRECTION].maxValue?wirePos[s][X_DIRECTION].maxValue:wirePos[s][Y_DIRECTION].maxValue;
		fCathMaxSum[s].val = wirePos[s][X_DIRECTION].maxValue+wirePos[s][Y_DIRECTION].maxValue;
		fPassedAnode[s] = fMWPC_anode[s].inRange();
		fPassedCath[s] = fCathSum[s].inRange();
		fPassedCathMax[s] = fCathMax[s].inRange();
		fPassedCathMaxSum[s] = fCathMaxSum[s].inRange();
	}
}

void ucnaDataAnalyzer11b::reconstructVisibleEnergy() {
	
	for(Side s = EAST; s <= WEST; ++s) {
		// get calibrated energy from the 4 tubes combined; also, wirechamber energy deposition estimate
		if(passedMWPC(s)) {
			PCal.calibrateEnergy(s,wirePos[s][X_DIRECTION].center,wirePos[s][Y_DIRECTION].center,sevt[s],fTimeScaler[BOTH]);
			// second pass with tweaked positions
			for(AxisDirection d = X_DIRECTION; d <= Y_DIRECTION; ++d)
				PCal.tweakPosition(s,d,wirePos[s][d],sevt[s].energy.x);
			PCal.calibrateEnergy(s,wirePos[s][X_DIRECTION].center,wirePos[s][Y_DIRECTION].center,sevt[s],fTimeScaler[BOTH]);
			fEMWPC[s] = PCal.calibrateAnode(fMWPC_anode[s].val,s,wirePos[s][X_DIRECTION].center,
											wirePos[s][Y_DIRECTION].center,fTimeScaler[BOTH]);
		} else {
			PCal.calibrateEnergy(s,0,0,sevt[s],fTimeScaler[BOTH]);
			fEMWPC[s] = PCal.calibrateAnode(fMWPC_anode[s].val,s,0,0,fTimeScaler[BOTH]);
		}
	}
}

void ucnaDataAnalyzer11b::checkMuonVetos() {
	fTaggedTop[WEST] = false;
	fTop_tdc[EAST].val = r_Top_TDC[EAST];
	fTop_adc[EAST] = r_Top_ADC[EAST];
	fTaggedTop[EAST] = fTop_tdc[EAST].inRange();
	for(Side s = EAST; s<=WEST; ++s) {
		fBacking_adc[s] = r_Backing_ADC[s];
		fBacking_tdc[s].val = r_Backing_TDC[s];
		fTaggedBack[s] = fBacking_tdc[s].inRange();
		fDrift_tac[s].val = r_Drift_TAC[s];
		fTaggedDrift[s] = fDrift_tac[s].inRange();
	}
}

void ucnaDataAnalyzer11b::classifyEventType() {
	// PID
	if(isLED()) fPID = PID_LED;	// LED event identified by Sis00
	else if(isPulserTrigger()) fPID = PID_PULSER;
	else if(Is2fold(EAST) || Is2fold(WEST)) {	// passes wirechamber and scintillator cuts on either side
		if(taggedMuon())
			fPID = PID_MUON; //at least one side muon
		else fPID = PID_BETA; //beta-like
	} else fPID = PID_SINGLE; //gamma
	
	// type, side
	fType = TYPE_IV_EVENT;
	fSide = NOSIDE;
	for(Side s = EAST; s<=WEST; ++s) {
		if(Is2fold(s)) {
			if(passedCutTDC(otherSide(s)))
				fType = TYPE_I_EVENT;
			else
				fType = passedMWPC(otherSide(s))?TYPE_II_EVENT:TYPE_0_EVENT;
		}
		if(passedCutTDC(s)&&!passedCutTDC(otherSide(s))) fSide = s;
	}
	// if side is ambiguous, TDCW has a cleaner TDC separation; make an 1-D cut (JL)
	if(passedCutTDC(WEST)&&passedCutTDC(EAST))
		fSide = (fScint_tdc[WEST][nBetaTubes].val < ScintSelftrig[WEST].start)?EAST:WEST;
	
	// Type II/III separation
	fProbIII = ((fType==TYPE_II_EVENT)?
				WirechamberCalibrator::sep23Prob(fSide, sevt[EAST].energy.x + sevt[WEST].energy.x, fEMWPC[fSide])
				: 0.);
	if(fProbIII>0.5) fType=TYPE_III_EVENT;
}

void ucnaDataAnalyzer11b::reconstructTrueEnergy() {
	if((fSide==EAST || fSide==WEST) && fType <= TYPE_III_EVENT)
		fEtrue = PCal.Etrue(fSide,fType,sevt[EAST].energy.x,sevt[WEST].energy.x);
	else
		fEtrue = sevt[EAST].energy.x + sevt[WEST].energy.x;
}

void ucnaDataAnalyzer11b::calcTrigEffic() {
	
	printf("\nCalculating trigger efficiency...\n");
	for(Side s = EAST; s <= WEST; ++s) {
		for(unsigned int t=0; t<nBetaTubes; t++) {
			// efficiency graph
			TGraphAsymmErrors gEffic(hTrigEffic[s][t][0]->GetNbinsX());
			gEffic.BayesDivide(hTrigEffic[s][t][1],hTrigEffic[s][t][0],"w");
			
			// scan for 50% point
			int b = gEffic.GetN();
			double midx,y;
			while(b > 1) {
				gEffic.GetPoint(--b,midx,y);
				if(y < 0.5)
					break;
			}
			
			// fit
			printf("Pre-fit threshold guess: %.1f\n",midx);
			TF1 efficfit("efficfit",&fancyfish,-5,150,4);
			
			efficfit.SetParameter(0,midx);
			efficfit.SetParLimits(0,0,100.0);
			efficfit.FixParameter(1,10.0);
			efficfit.FixParameter(2,10.0);
			efficfit.FixParameter(3,0.999);
			efficfit.SetLineColor(38);
			gEffic.Fit(&efficfit,"QR");
			
			efficfit.SetParameter(1,10.0);
			efficfit.SetParLimits(1,2,200.0);			
			efficfit.SetLineColor(7);
			gEffic.Fit(&efficfit,"QR+");
			
			efficfit.SetParameter(2,10.0);
			efficfit.SetParameter(3,0.999);
			efficfit.SetParLimits(2,0.1,1000.0);
			efficfit.SetParLimits(3,0.75,1.0);
			efficfit.SetLineColor(4);
			gEffic.Fit(&efficfit,"QR+");
			
			
			float_err trigef(efficfit.GetParameter(3),efficfit.GetParError(3));
			float_err trigc(efficfit.GetParameter(0),efficfit.GetParError(0));
			float_err trigw(efficfit.GetParameter(1),efficfit.GetParError(1));
			float_err trign_adj(efficfit.GetParameter(2),efficfit.GetParError(2));
			float trign = trigc.x/trigw.x*trign_adj.x;
			
			// save results
			printf("Poisson CDF Fit: h = %.4f(%.4f), x0 = %.1f(%.1f), dx = %.1f(%.1f), n = %.2f [adjust %.2f(%.2f)]\n",
				   trigef.x, trigef.err, trigc.x, trigc.err, trigw.x, trigw.err, trign, trign_adj.x, trign_adj.err);
			Stringmap m;
			m.insert("effic_params",vtos(efficfit.GetParameters(),efficfit.GetParameters()+4));
			m.insert("effic_params_err",vtos(efficfit.GetParErrors(),efficfit.GetParErrors()+4));
			m.insert("side",ctos(sideNames(s)));
			m.insert("tube",t);
			qOut.insert("trig_effic",m);
			// upload to analysis DB
			if(CDBout) {
				printf("Uploading trigger efficiency...\n");
				std::vector<double> tparams;
				std::vector<double> terrs;
				for(unsigned int i=0; i<4; i++) {
					tparams.push_back(efficfit.GetParameter(i));
					terrs.push_back(efficfit.GetParError(i));
				}
				CDBout->deleteTrigeff(rn,s,t);
				CDBout->uploadTrigeff(rn,s,t,tparams,terrs);
			}
			
			// plot
			gEffic.SetMinimum(-0.10);
			gEffic.SetMaximum(1.10);
			gEffic.Draw("AP");
			gEffic.SetTitle((sideSubst("%c",s)+itos(t)+" PMT Trigger Efficiency").c_str());
			gEffic.GetXaxis()->SetTitle("ADC channels above pedestal");
			gEffic.GetXaxis()->SetLimits(-50,150);
			gEffic.GetYaxis()->SetTitle("Efficiency");
			gEffic.Draw("AP");
			printCanvas(sideSubst("PMTs/TrigEffic_%c",s)+itos(t));
		}
	}
}

bool ucnaDataAnalyzer11b::processEvent() {
	// Stage I: uses read in data, moves to appropriate locations
	convertReadin();
	checkHeaderQuality();
	calibrateTimes();
	checkMuonVetos();
	
	// load data for next point for look-ahead capability; overwrites r_* variables after this point
	bool np = nextPoint();
	fWindow.val = fDelt0 + 1.e-6*r_Delt0;
	
	// Stage II: processing and histograms for all events
	for(Side s = EAST; s <= WEST; ++s)
		PCal.pedSubtract(s, sevt[s].adc, fTimeScaler[BOTH]);
	fillEarlyHistograms();
	
	// LED tree events
	if(isLED() && TLED) {
		reconstructPosition();
		for(Side s = EAST; s <= WEST; ++s)
			PCal.calibrateEnergy(s,0.,0.,sevt[s],fTimeScaler[BOTH]);
		TLED->Fill();
	}
	
	if(!isScintTrigger() || isLED())
		return np;
	
	// Stage III: processing and histograms for beta trigger events
	reconstructPosition();
	reconstructVisibleEnergy();
	classifyEventType();
	reconstructTrueEnergy();
	fillHistograms();
	
	if(fPassedGlobal)
		TPhys->Fill();
	
	return np;
}

void ucnaDataAnalyzer11b::processBiPulser() {
	
	printf("\nFitting Bi Pulser...\n");
	TF1 gausFit("gasufit","gaus",1000,4000);
	QFile pulseLocation(dataPath+"/Monitors/Run_"+itos(rn)+"/ChrisPulser.txt",false);
	for(Side s = EAST; s <= WEST; ++s) {
		for(unsigned int t=0; t<nBetaTubes; t++) {
			std::vector<double> times;
			std::vector<double> centers;
			std::vector<double> dcenters;
			std::vector<double> widths;
			std::vector<double> dwidths;
			for(unsigned int i=0; i<hBiPulser[s][t].size(); i++) {
				Stringmap m;
				m.insert("side",ctos(sideNames(s)));
				m.insert("tube",t);
				m.insert("counts",hBiPulser[s][t][i]->GetEntries());
				
				// initial estimate of peak location: first high isolated peak scanning from right
				int bmax = 0;
				float pmax = 0;
				unsigned int nskip = 0;
				for(unsigned int n = hBiPulser[s][t][i]->GetNbinsX(); n > 0; n--) {
					if(hBiPulser[s][t][i]->GetBinContent(n)>pmax) {
						bmax = n;
						pmax = hBiPulser[s][t][i]->GetBinContent(n);
						nskip = 0;
					} else {
						nskip++;
					}
					if(pmax>hBiPulser[s][t][i]->GetMaximum()*0.3 && nskip>20)
						break;
				}
				float bcenter = hBiPulser[s][t][i]->GetBinCenter(bmax);
				
				// refined fit
				if(!iterGaus(hBiPulser[s][t][i],&gausFit,3,bcenter,200,1.5)) {
					times.push_back((i+0.5)*wallTime/hBiPulser[s][t].size());
					m.insert("time",times.back());
					m.insert("height",gausFit.GetParameter(0));
					m.insert("dheight",gausFit.GetParError(0));
					centers.push_back(gausFit.GetParameter(1));
					m.insert("center",centers.back());
					dcenters.push_back(gausFit.GetParError(1));
					m.insert("dcenter",dcenters.back());
					widths.push_back(gausFit.GetParameter(2));
					m.insert("width",widths.back());
					dwidths.push_back(gausFit.GetParError(2));
					m.insert("dwidth",dwidths.back());
				}
				
				pulseLocation.insert("pulserpeak",m);
			}
			if(CDBout) {
				// upload to DB
				std::string mon_name = PCal.sensorNames[s][t];
				printf("Uploading pulser data '%s'...\n",mon_name.c_str());
				unsigned int cgid = CDBout->uploadGraph(itos(rn)+" "+mon_name+" Pulser Centers",times,centers,std::vector<double>(),dcenters);
				unsigned int wgid = CDBout->uploadGraph(itos(rn)+" "+mon_name+" Pulser Widths",times,widths,std::vector<double>(),dwidths);
				CDBout->deleteRunMonitor(rn,mon_name,"Chris_peak");
				CDBout->addRunMonitor(rn,mon_name,"Chris_peak",cgid,wgid);				
			}
			drawSimulHistos(hBiPulser[s][t]);
			printCanvas(sideSubst("PMTs/BiPulser_%c",s)+itos(t));
		}
	}
	pulseLocation.commit();
}

void ucnaDataAnalyzer11b::tallyRunTime() {
	// complete potential un-closed blip
	if(cutBlips.size() && cutBlips.back().end[BOTH] == 0)
		cutBlips.back().end = totalTime;
	// sum up lost time	
	BlindTime lostTime;
	for(std::vector<Blip>::iterator it = cutBlips.begin(); it != cutBlips.end(); it++)
		lostTime += it->length();
	wallTime = totalTime[BOTH]; // now an informed guess :) 
	totalTime -= lostTime;
	printf("\nFiducial time tally:\n");
	printf("Lost %.1fs run time to %i blips, leaving %.1fs. (%g,%g failed Evnb,Bkhf)\n",
		   lostTime[BOTH],(int)cutBlips.size(),totalTime[BOTH],nFailedEvnb,nFailedBkhf);
}

void ucnaDataAnalyzer11b::replaySummary() {
	if(!CDBout) return;
	sprintf(CDBout->query,"DELETE FROM analysis WHERE run_number = %i",rn);
	CDBout->execute();
	TDatime tNow;
	sprintf(CDBout->query,
			"INSERT INTO analysis(run_number,analysis_time,live_time_e,live_time_w,live_time,total_time,n_trigs,total_trigs,misaligned,tdc_corrupted) \
			VALUES (%i,'%s',%f,%f,%f,%f,%u,%u,%i,%i)",
			int(rn),tNow.AsSQLString(),totalTime[EAST],totalTime[WEST],totalTime[BOTH],
			wallTime,nLiveTrigs,nEvents,int(nFailedEvnb),int(nFailedBkhf));
	CDBout->execute();
	TDatime tStart(fAbsTimeStart);
	std::string sts(tStart.AsSQLString());
	TDatime tEnd(fAbsTimeEnd);
	std::string ste(tEnd.AsSQLString());
	sprintf(CDBout->query,"UPDATE run SET start_time='%s', end_time='%s' WHERE run_number=%i",sts.c_str(),ste.c_str(),int(rn));
	printf("%s\n",CDBout->query);
	CDBout->execute();
}

void ucnaDataAnalyzer11b::quickAnalyzerSummary() const {
	printf("\n--------------- Quick Summary %i ---------------\n",rn);
	float gvcounts = hMonADC[UCN_MON_GV]->Integral();
	printf("GV Mon: %i = %.2f +/- %.2f Hz\n",(int)gvcounts,gvcounts/wallTime,sqrt(10*gvcounts)/10/wallTime);
	float fecounts = hMonADC[UCN_MON_FE]->Integral();
	printf("Fe Mon: %i = %.2f Hz; Fe/GV = %.4f\n",(int)fecounts,fecounts/wallTime,fecounts/gvcounts);
	float scscounts = hMonADC[UCN_MON_SCS]->Integral();
	printf("SCS Mon: %i = %.2f Hz; SCS/GV = %.4f\n",(int)scscounts,scscounts/wallTime,scscounts/gvcounts);
	float scounts[2];
	for(Side s = EAST; s <= WEST; ++s) {
		scounts[s] = hEtrue[s][TYPE_0_EVENT]->Integral()+hEtrue[s][TYPE_I_EVENT]->Integral()+hEtrue[s][TYPE_II_EVENT]->Integral();
		printf("%s Beta Triggers: %i = %.2f +/- %.2f Hz\n",sideWords(s),(int)scounts[s],scounts[s]/wallTime,sqrt(scounts[s])/wallTime);
	}
	printf("Beta/GV = %.2f\n",(scounts[EAST]+scounts[WEST])/gvcounts);
	printf("Bonehead (E-W)/(E+W) = %.2f%%\n",100.0*(scounts[EAST]-scounts[WEST])/(scounts[EAST]+scounts[WEST]));
	TDatime stime;
	stime.Set(fAbsTimeEnd);
	printf("----------------------------------------------------\n");
	printf("%s\t%i\t%i\t%i\t%.1f\n",stime.AsSQLString(),rn,(int)scounts[EAST],(int)scounts[WEST],wallTime);
	printf("----------------------------------------------------\n\n");
}

void printHelp(const char* argname) {
	printf("Syntax: %s <run number(s)> [options...]\n",argname);
	printf("\t<run number(s)>: <number> || <number>-<number>\n");
	printf("\toptions:\n");
	printf("\t\tcutbeam: cut out beam pulses and when GV rate too low (use this for beta/background runs!)\n");
	printf("\t\tnodbout: do not access Calibrations DB for writing\n");
	printf("\t\tnoroot: skip saving output .root file (plots/summary only)\n");
	printf("\t\tledtree: produce separate TTree with only LED events\n");
	printf("\t\tforceped: force recalculation of pedestals\n");
}

int main(int argc, char** argv) {
	
	// check correct arguments
	if(argc<2) {
		printHelp(argv[0]);
		exit(1);
	}
	
	// get run(s)
	std::vector<int> rlist = sToInts(argv[1],"-");
	if(!rlist.size() || !rlist[0] || rlist.size()>2) {
		printf("*** '%s' is not a valid run number! Exiting!\n",argv[1]);
		exit(1);
	}
	if(rlist.size()==1)
		rlist.push_back(rlist[0]);
	
	// other options
	bool cutBeam = false;
	bool nodbout = false;
	bool noroot = false;
	bool ledtree = false;
	bool forceped = false;
	for(int i=2; i<argc; i++) {
		std::string arg(argv[i]);
		if(arg=="cutbeam")
			cutBeam = true;
		else if(arg=="nodbout")
			nodbout = true;
		else if(arg=="noroot")
			noroot = true;
		else if(arg=="ledtree")
			ledtree = true;
		else if(arg=="forceped")
			forceped = true;
		else {
			printHelp(argv[0]);
			exit(1);
		}
	}
	
	gStyle->SetPalette(1);
	gStyle->SetNumberContours(255);
	gStyle->SetOptStat("e");
	
	std::string outDir = getEnvSafe("UCNAOUTPUTDIR");
	
	for(RunNum r = (unsigned int)rlist[0]; r<=(unsigned int)rlist[1]; r++) {
		
		std::string inDir = getEnvSafe("UCNADATADIR");
		if(!fileExists(inDir+"/full"+itos(r)+".root") && r > 16300)
			inDir = "/data/ucnadata/2011/rootfiles/";
		
		ucnaDataAnalyzer11b A(r,outDir,CalDBSQL::getCDB(true));
		A.setIgnoreBeamOut(!cutBeam);
		A.analyzeLED = ledtree;
		A.needsPeds = forceped;
		if(!nodbout) {
			printf("Connecting to output DB...\n");
			A.setOutputDB(CalDBSQL::getCDB(false));
		}
		A.addFile(inDir+"/full"+itos(r)+".root");
		A.analyze();
		A.setWriteRoot(!noroot);
		A.write();
	}
	
	return 0;
}
