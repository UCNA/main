#include "ucnaAnalyzerBase.hh"

ucnaAnalyzerBase::ucnaAnalyzerBase(RunNum R, const std::string& bp, const std::string& nm, CalDB* CDB):
TChainScanner("h1"), OutputManager(nm+"_"+itos(R), bp+"/hists/"), rn(R), PCal(R,CDB),
fAbsTimeEnd(0), totalTime(0), deltaT(0), ignore_beam_out(false) {

	// beta scintillator TDC timing cuts
	for(Side s = EAST; s <= WEST; ++s) {
		loadRangeCut(rn,fScint_tdc[s][nBetaTubes], sideSubst("Cut_TDC_Scint_%c_Selftrig",s));
		ScintSelftrig[s] = fScint_tdc[s][nBetaTubes].R;
		loadRangeCut(rn,fScint_tdc[s][nBetaTubes], sideSubst("Cut_TDC_Scint_%c",s));
		for(unsigned int t=0; t<nBetaTubes; t++)
			loadRangeCut(rn,fScint_tdc[s][t], sideSubst("Cut_TDC_Scint_%c_",s)+itos(t));
	}
	
	// input tree branches
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
	
	// beam burst time cut
	loadRangeCut(rn,fBeamclock,"Cut_BeamBurst");
	if(ignore_beam_out)
		fBeamclock.R.end = FLT_MAX;
	// manually excluded times
	manualCuts = ManualInfo::MI.getRanges(itos(rn)+"_timecut");
	if(manualCuts.size())
		printf("Manually cutting %i time ranges...\n",(int)manualCuts.size());

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
