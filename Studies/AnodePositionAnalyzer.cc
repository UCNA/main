#include "AnodePositionAnalyzer.hh"

AnodePositionAnalyzer::AnodePositionAnalyzer(RunAccumulator* RA, unsigned int nr): PositionBinnedAnalyzer(RA,"AnodePos",nr) {
	TH1F hTemplate("hAnode","Anode Signal",175,0,3500);
	for(Side s = EAST; s <= WEST; ++s)
		sectAnode[s] = allocateSegmentHistograms(hTemplate,AFP_OTHER,s);
}

void AnodePositionAnalyzer::fillCoreHists(ProcessedDataScanner& PDS, double weight) {
	Side s = PDS.fSide;
	if(!(PDS.fType == TYPE_0_EVENT && PDS.fPID == PID_BETA && PDS.getEnergy() > 250)) return;
	unsigned int m = sects.sector(PDS.wires[s][X_DIRECTION].center,PDS.wires[s][Y_DIRECTION].center);
	if(m>=sects.nSectors()) return;
	sectAnode[s][m]->h[currentGV]->Fill(PDS.mwpcs[s].anode,weight);
}

void AnodePositionAnalyzer::genAnodePosmap() {
	std::string pmapname = "Anode_"+itos(myA->runCounts.counts.begin()->first)+"-"+itos(myA->runCounts.counts.rbegin()->first)+"/"+itos(time(NULL));
	CalDBSQL* CDBout = CalDBSQL::getCDB(false);
	unsigned int pmid = CDBout->newPosmap(pmapname,sects.n,sects.r);
	float x,y;
	TF1 fLandau("landauFit","landau",0,3500);
	fLandau.SetParLimits(0, 0, 2e5);
	fLandau.SetParameter(2,50.);
	fLandau.SetParLimits(2,20.,100.);
	for(Side s=EAST; s<=WEST; ++s) {
		for(unsigned int m=0; m<sects.nSectors(); m++) {
			printf("\n------------- %s %i ------------------\n",sideWords(s),m);
			sects.sectorCenter(m,x,y);
			TH1* hAnode = sectAnode[s][m]->h[GV_OPEN];
			double c = hAnode->GetBinCenter(hAnode->GetMaximumBin());
			double mx = hAnode->GetMaximum();
			printf("Max = %.3g at %.1f\n",mx,c);
			fLandau.SetRange(c-300>0?0:c-300, c+500);
			fLandau.SetParameter(0,mx*c/fLandau.GetParameter(2));
			fLandau.SetParameter(1,c);
			hAnode->Fit(&fLandau,"ERMB");
			CDBout->addPosmapPoint(pmid,s,0,m,fLandau.GetParameter(1),1.0,x,y);
		}
	}
	Stringmap pmsm;
	pmsm.insert("pmid",pmid);
	pmsm.insert("name",pmapname);
	myA->qOut.insert("anode_posmap",pmsm);
}
