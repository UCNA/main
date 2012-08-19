#include "ucnaDataAnalyzer11b.hh"

void ucnaDataAnalyzer11b::setReadpoints() {
	
	// clock, trigger
	SetBranchAddress("Sis00",&r_Sis00);				// trigger type flags Sis00
	SetBranchAddress("Number",&r_TriggerNumber);	// trigger number
	SetBranchAddress("Clk0",&r_Clk[EAST]);		// E clock
	SetBranchAddress("Clk1",&r_Clk[WEST]);		// W clock
	SetBranchAddress("S83028",&r_Clk[BOTH]);		// unblinded runclock
	SetBranchAddress("S8200",&r_BClk);				// protonClock
	SetBranchAddress("Delt0",&r_Delt0);				// high resolution delta-t since previous event
	SetBranchAddress("Time",&r_AbsTime);			// absolute time during run
	
	//TCh->SetBranchAddress("S8300","");					// E Beta Counter
	//TCh->SetBranchAddress("S8301","");					// W Beta Counter
	//TCh->SetBranchAddress("S8302","");					// UCN Mon 1 trigger counter
	//TCh->SetBranchAddress("S8303","");					// UCN Mon 2 trigger counter
	//TCh->SetBranchAddress("S8304","");					// UCN Mon 3 trigger counter
	//TCh->SetBranchAddress("S8308","");					// E Gated trigger counter
	//TCh->SetBranchAddress("S8309","");					// W Gated trigger counter
	
	// ucn monitors
	const int ucn_mon_adc_nums[kNumUCNMons] = {38,39,310,311};
	for(unsigned int n=0; n<kNumUCNMons; n++)
		SetBranchAddress("Pdc"+itos(ucn_mon_adc_nums[n]),&r_MonADC[n]);
	
	// PMTs and TDCs, in quadrant order (+x towards SCS, +y up, +z from East to West)
	const int pmt_adc_nums[] = {2,3,0,1,4,5,6,7};
	const int pmt_tdc_nums[] = {2,3,0,1,8,9,14,11};
	for(Side s = EAST; s <= WEST; ++s) {
		for(unsigned int t=0; t<nBetaTubes; t++) {
			SetBranchAddress("Qadc"+itos(pmt_adc_nums[t+4*s]),&r_PMTADC[s][t]);
			SetBranchAddress("Tdc0"+itos(pmt_tdc_nums[t+4*s]),&r_PMTTDC[s][t]);
		}
	}
	SetBranchAddress("Tdc016",&r_PMTTDC[EAST][nBetaTubes]);		// E 2of4 TDC
	SetBranchAddress("Tdc017",&r_PMTTDC[WEST][nBetaTubes]);		// W 2of4 TDC
	
	// Wirechambers	
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
	
	// muon vetos
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
	
	// header checks
	for(size_t i=0; i<kNumModules; i++) {
		SetBranchAddress("Evnb"+itos(i),&r_Evnb[i]);
		SetBranchAddress("Bkhf"+itos(i),&r_Bkhf[i]);
	}
	
	
}

void ucnaDataAnalyzer11b::setupOutputTree() {
	openOutfile();
	TPhys = (TTree*)addObject(new TTree("phys","physics quantities"));
	TPhys->SetMaxVirtualSize(1000000);
	
	TPhys->Branch("TriggerNum",&iTriggerNumber,"TriggerNum/I");
	TPhys->Branch("Sis00",&iSis00,"Sis00/I");
	TPhys->Branch("DeltaT",&fDelt0,"DeltaT/F");
	TPhys->Branch("EvtN",&currentEvent,"EvtN/I");
	
	for(Side s = EAST; s <= WEST; ++s) {
		TPhys->Branch(sideSubst("Time%c",s).c_str(),&fTimeScaler[s],sideSubst("Time%c/F",s).c_str());
		TPhys->Branch(sideSubst("Tof%c",s).c_str(),&fBeamclock.val,sideSubst("Tof%c/F",s).c_str());
		
		TPhys->Branch(sideSubst("TDC%c",s).c_str(),&fScint_tdc[s][nBetaTubes].val,sideSubst("TDC%c/F",s).c_str());
		for(unsigned int t=0; t<nBetaTubes; t++)
			TPhys->Branch((sideSubst("TDC%c",s)+itos(t+1)).c_str(),&fScint_tdc[s][t].val,(sideSubst("TDC%c",s)+itos(t+1)+"/F").c_str());
		
		for(AxisDirection d = X_DIRECTION; d <= Y_DIRECTION; ++d) {
			TPhys->Branch((d==X_DIRECTION?sideSubst("x%cmpm",s):sideSubst("y%cmpm",s)).c_str(),
						  &wirePos[s][d],"center/F:width/F:maxValue/F:cathSum/F:maxWire/i:nClipped/i:mult/i:err/I:rawCenter/F");
			std::string cathname = sideSubst("Cathodes_%c",s)+(d==X_DIRECTION?"x":"y");
			TPhys->Branch(cathname.c_str(),&fMWPC_caths[s][d][0],(cathname+"["+itos(kMaxCathodes)+"]/F").c_str());
		}
		
		TPhys->Branch(sideSubst("Scint%c",s).c_str(),&sevt[s],
					  "q1/F:q2/F:q3/F:q4/F:e1/F:de1/F:e2/F:de2/F:e3/F:de3/F:e4/F:de4/F:energy/F:denergy/F:nPE1/F:nPE2/F:nPE3/F:nPE4/F");
		
		TPhys->Branch(sideSubst("EMWPC_%c",s).c_str(),&fEMWPC[s],sideSubst("EMWPC_%c/F",s).c_str());
		
		TPhys->Branch(sideSubst("Anode%c",s).c_str(),&fMWPC_anode[s].val,sideSubst("ano%c/F",s).c_str());
		TPhys->Branch(sideSubst("Evis%c",s).c_str(),&sevt[s].energy.x,sideSubst("Evis%c/F",s).c_str());
		TPhys->Branch(sideSubst("CathSum%c",s).c_str(),&fCathSum[s].val,sideSubst("CathSum%c/F",s).c_str());
		TPhys->Branch(sideSubst("CathMax%c",s).c_str(),&fCathMax[s].val,sideSubst("CathMax%c/F",s).c_str());
		TPhys->Branch(sideSubst("PassedAno%c",s).c_str(),&fPassedAnode[s],sideSubst("PassedAno%c/I",s).c_str());
		TPhys->Branch(sideSubst("PassedCath%c",s).c_str(),&fPassedCath[s],sideSubst("PassedCath%c/I",s).c_str());
		
		TPhys->Branch(sideSubst("TaggedBack%c",s).c_str(),&fTaggedBack[s],sideSubst("TaggedBack%c/I",s).c_str());
		TPhys->Branch(sideSubst("TaggedTop%c",s).c_str(),&fTaggedTop[s],sideSubst("TaggedTop%c/I",s).c_str());
		TPhys->Branch(sideSubst("TaggedDrift%c",s).c_str(),&fTaggedDrift[s],sideSubst("TaggedDrift%c/I",s).c_str());
		TPhys->Branch(sideSubst("%sBackADC",s).c_str(),&fBacking_adc[s],sideSubst("%sBackADC/F",s).c_str());
		TPhys->Branch(sideSubst("%sBackTDC",s).c_str(),&fBacking_tdc[s],sideSubst("%sBackTDC/F",s).c_str());
		TPhys->Branch(sideSubst("%sDriftVetoADC",s).c_str(),&fDrift_tac[s].val,sideSubst("%sDriftVetoADC/F",s).c_str());
	}
	TPhys->Branch("EastTopVetoADC",&fTop_adc[EAST],"EastTopVetoADC/F");
	TPhys->Branch("EastTopVetoTDC",&fTop_tdc[EAST].val,"EastTopVetoTDC/F");
	
	TPhys->Branch("EvnbGood",&fEvnbGood,"EvnbGood/I");
	TPhys->Branch("BkhfGood",&fBkhfGood,"BkhfGood/I");
	
	TPhys->Branch("PID",&fPID,"PID/I");
	TPhys->Branch("Type",&fType,"Type/I");
	TPhys->Branch("Side",&fSide,"Side/I");
	TPhys->Branch("ProbIII",&fProbIII,"ProbIII/F");
	TPhys->Branch("Etrue",&fEtrue,"Etrue/F");
	
	// LED events tree
	if(analyzeLED) {
		TLED = (TTree*)addObject(new TTree("LED","LED events"));
		TLED->SetMaxVirtualSize(1000000);
		TLED->Branch("EvtN",&currentEvent,"EvtN/I");
		TLED->Branch("Time",&fTimeScaler[BOTH],"Time/F");
		for(Side s = EAST; s <= WEST; ++s) {
			TLED->Branch(sideSubst("Scint%c",s).c_str(),&sevt[s],
						 "q1/F:q2/F:q3/F:q4/F:e1/F:de1/F:e2/F:de2/F:e3/F:de3/F:e4/F:de4/F:energy/F:denergy/F:nPE1/F:nPE2/F:nPE3/F:nPE4/F");
			TLED->Branch(sideSubst("Anode%c",s).c_str(),&fMWPC_anode[s].val,sideSubst("ano%c/F",s).c_str());
			TLED->Branch(sideSubst("CathMax%c",s).c_str(),&fCathMax[s].val,sideSubst("CathMax%c/F",s).c_str());
		}
	} else {
		TLED = NULL;
	}
}
