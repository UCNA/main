#include "ucnaDataAnalyzer11b.hh"

void ucnaDataAnalyzer11b::setReadpoints() {
	readInTiming();
	readInWirechambers();
	readInPMTADC();
	readInUCNMon();
	readInMuonVetos();
	readInHeaderChecks();
}

void ucnaDataAnalyzer11b::setupOutputTree() {
	openOutfile();
	TPhys = (TTree*)addObject(new TTree("phys","physics quantities"));
	TPhys->SetMaxVirtualSize(1000000);
	
	TPhys->Branch("TriggerNum",&iTriggerNumber,"TriggerNum/I");
	TPhys->Branch("Sis00",&SIS00,"Sis00/I");
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
						  &wirePos[s][d],"center/F:width/F:maxValue/F:cathSum/F:maxWire/i:nClipped/i:mult/i:err/I:rawCenter/F:height/F");
			std::string cathname = sideSubst("Cathodes_%c",s)+(d==X_DIRECTION?"x":"y");
			TPhys->Branch(cathname.c_str(),r_MWPC_caths[s][d],(cathname+"["+itos(kMaxCathodes)+"]/F").c_str());
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
	TPhys->Branch("Erecon",&fEtrue,"Erecon/F");
	
	// LED events tree
	if(analyzeLED) {
		TLED = (TTree*)addObject(new TTree("LED","LED events"));
		TLED->SetMaxVirtualSize(1000000);
		TLED->Branch("EvtN",&currentEvent,"EvtN/I");
		TLED->Branch("Time",&fTimeScaler[BOTH],"Time/F");
		TLED->Branch("Sis00",&SIS00,"Sis00/I");
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
