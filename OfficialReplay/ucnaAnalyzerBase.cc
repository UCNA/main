#include "ucnaAnalyzerBase.hh"

ucnaAnalyzerBase::ucnaAnalyzerBase(RunNum R, const std::string& bp, const std::string& nm, CalDB* CDB):
TChainScanner("h1"), OutputManager(nm+"_"+itos(R), bp+"/hists/"), rn(R), PCal(R,CDB),
fAbsTimeEnd(0), totalTime(0), deltaT(0), ignore_beam_out(false) {
  printf("I'm in here!\n");
	// beta scintillator TDC timing cuts
	for(Side s = EAST; s <= WEST; ++s) {
		loadRangeCut(rn,fScint_tdc[s][nBetaTubes], sideSubst("Cut_TDC_Scint_%c_Selftrig",s));
		ScintSelftrig[s] = fScint_tdc[s][nBetaTubes].R;
		loadRangeCut(rn,fScint_tdc[s][nBetaTubes], sideSubst("Cut_TDC_Scint_%c",s));
		for(unsigned int t=0; t<nBetaTubes; t++)
			loadRangeCut(rn,fScint_tdc[s][t], sideSubst("Cut_TDC_Scint_%c_",s)+itos(t));
	}
			
	// beam burst time cut
	loadRangeCut(rn,fBeamclock,"Cut_BeamBurst");
	
	// wirechamber and veto cuts
	for(Side s = EAST; s <= WEST; ++s) {
		loadRangeCut(rn,fMWPC_anode[s], sideSubst("Cut_MWPC_%c_Anode",s));
		loadRangeCut(rn,fCathMax[s],sideSubst("Cut_MWPC_%c_CathMax",s));
		loadRangeCut(rn,fCathSum[s],sideSubst("Cut_MWPC_%c_CathSum",s));
		loadRangeCut(rn,fCathMaxSum[s],sideSubst("Cut_MWPC_%c_CathMaxSum",s));
		loadRangeCut(rn,fBacking_tdc[s], sideSubst("Cut_TDC_Back_%c",s));
		loadRangeCut(rn,fDrift_tac[s], sideSubst("Cut_ADC_Drift_%c",s));
	}
	loadRangeCut(rn,fTop_tdc[EAST], "Cut_TDC_Top_E");
	loadRangeCut(rn,fWindow,"Cut_ClusterEvt");
	
	// manually excluded times
	manualCuts = ManualInfo::MI.getRanges(itos(rn)+"_timecut");
	if(manualCuts.size())
		printf("Manually cutting %i time ranges...\n",(int)manualCuts.size());

}

void ucnaAnalyzerBase::setIgnoreBeamOut(bool ibo) {
	ignore_beam_out = ibo;
	if(ignore_beam_out)
		fBeamclock.R.end = FLT_MAX;
}

void ucnaAnalyzerBase::readInTiming() {
	SetBranchAddress("Sis00",&r_Sis00);				// trigger type flags Sis00
	SetBranchAddress("Number",&r_TriggerNumber);	// trigger number
	SetBranchAddress("Clk0",&r_Clk[EAST]);			// E clock
	SetBranchAddress("Clk1",&r_Clk[WEST]);			// W clock
	SetBranchAddress("S83028",&r_Clk[BOTH]);		// unblinded runclock
	SetBranchAddress("S8200",&r_BClk);				// protonClock
	SetBranchAddress("Delt0",&r_Delt0);				// high resolution delta-t since previous event
	SetBranchAddress("Time",&r_AbsTime);			// absolute time during run
	// PMTs and TDCs, in quadrant order (+x towards SCS, +y up, +z from East to West)
	const int pmt_tdc_nums[] = {2,3,0,1, 16, 8,9,14,11, 17}; // including 2of4 TDCs 16,17
	for(Side s = EAST; s <= WEST; ++s)
		for(unsigned int t=0; t<=nBetaTubes; t++)
			SetBranchAddress("Tdc0"+itos(pmt_tdc_nums[t+(nBetaTubes+1)*s]),&r_PMTTDC[s][t]);
}

void ucnaAnalyzerBase::readInWirechambers() {
	const int anode_pdc_nums[] = {30,34};
	for(Side s = EAST; s<=WEST; ++s) {
		for(AxisDirection d = X_DIRECTION; d <= Y_DIRECTION; ++d) {
			cathNames[s][d].clear();
			std::vector<std::string> cchans = PCal.getCathChans(s,d);
			for(unsigned int i=0; i<cchans.size(); i++) {
				SetBranchAddress(cchans[i],&r_MWPC_caths[s][d][i]);
				cathNames[s][d].push_back(sideSubst("MWPC%c",s)+(d==X_DIRECTION?"x":"y")+itos(i+1));
			}
		}
		SetBranchAddress("Pdc"+itos(anode_pdc_nums[s]),&r_MWPC_anode[s]);
	}
}

void ucnaAnalyzerBase::readInUCNMon() {
	const int ucn_mon_adc_nums[kNumUCNMons] = {38,39,310,311};
	for(unsigned int n=0; n<kNumUCNMons; n++)
		SetBranchAddress("Pdc"+itos(ucn_mon_adc_nums[n]),&r_MonADC[n]);
}

void ucnaAnalyzerBase::readInPMTADC() {
	// PMTs and TDCs, in quadrant order (+x towards SCS, +y up, +z from East to West)
	const int pmt_adc_nums[] = {2,3,0,1,4,5,6,7};
	for(Side s = EAST; s <= WEST; ++s)
		for(unsigned int t=0; t<nBetaTubes; t++)
			SetBranchAddress("Qadc"+itos(pmt_adc_nums[t+4*s]),&r_PMTADC[s][t]);
}

void ucnaAnalyzerBase::readInMuonVetos() {
	const int back_tdc_nums[] = {18,20};
	const int back_adc_nums[] = {8,10};
	const int drift_tac_nums[] = {313,315};
	for(Side s = EAST; s<=WEST; ++s) {
		SetBranchAddress("Tdc0"+itos(back_tdc_nums[s]),&r_Backing_TDC[s]);
		SetBranchAddress("Pdc"+itos(drift_tac_nums[s]),&r_Drift_TAC[s]);
		SetBranchAddress("Qadc"+itos(back_adc_nums[s]),&r_Backing_ADC[s]);
	}
	SetBranchAddress("Tdc019",&r_Top_TDC[EAST]);
	SetBranchAddress("Qadc9",&r_Top_ADC[EAST]);
}

void ucnaAnalyzerBase::readInHeaderChecks() {
	for(size_t i=0; i<kNumModules; i++) {
		SetBranchAddress("Evnb"+itos(i),&r_Evnb[i]);
		SetBranchAddress("Bkhf"+itos(i),&r_Bkhf[i]);
	}
}

unsigned int ucnaAnalyzerBase::nFiring(Side s) const {
	unsigned int nf = 0;
	for(unsigned int t=0; t<nBetaTubes; t++) nf += pmtFired(s,t);
	return nf;
}

void ucnaAnalyzerBase::calibrateTimes() {
	
	// start/end times
	if(currentEvent==0) {
		fAbsTimeStart = r_AbsTime;
		deltaT = 0;
		totalTime = 0;
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
	totalTime = fTimeScaler;
	
	// scintillator TDC read-in
	for(Side s = EAST; s <= WEST; ++s)
		for(unsigned int t=0; t<=nBetaTubes; t++)
			fScint_tdc[s][t].val = r_PMTTDC[s][t];
}

bool ucnaAnalyzerBase::passesBeamCuts() const {
	// basic time-since-beam cut
	if(!fBeamclock.inRange())
		return false;	
	// remove manually tagged segments
	for(std::vector< std::pair<double,double> >::const_iterator it = manualCuts.begin(); it != manualCuts.end(); it++)
		if (it->first <= fTimeScaler[BOTH] && fTimeScaler[BOTH] <= it->second)
			return false;
	return true;
}
