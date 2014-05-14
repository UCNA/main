#include "HighEnergyExcessPlugin.hh"
#include "GraphUtils.hh"
#include "GraphicsUtils.hh"

void fitHighEnergyExcess(RunAccumulator* RA, quadHists* qh, double e0, double e1) {
	for(AFPState afp = AFP_OFF; afp <= AFP_ON; ++afp) {
		TH1* hEn = qh->fgbg[afp]->h[GV_OPEN];
		TH1* hEnBG = qh->fgbg[afp]->h[GV_CLOSED];
		
		int b0 = hEnBG->FindBin(e0);
		int b1 = hEnBG->FindBin(e1)-1;
		if(b1 == hEnBG->GetNbinsX()) ++b1;
		Double_t d_nBG;
		double nBG = hEnBG->IntegralAndError(b0,b1,d_nBG);
		Double_t d_xs;
		double xs = hEn->IntegralAndError(b0,b1,d_xs);
		
		AnaNumber AN("bg_xs_"+qh->name+"_"+itos(e0)+"-"+itos(e1));
		AN.s = qh->mySide;
		AN.value = xs;
		AN.err = d_xs;
		RA->uploadAnaNumber(AN, GV_CLOSED, afp);
		
		AN.name = "bg_cts_"+qh->name+"_"+itos(e0)+"-"+itos(e1);
		AN.value = nBG;
		AN.err = d_nBG;
		RA->uploadAnaNumber(AN, GV_CLOSED, afp);
	}
}

HighEnergyExcessPlugin::HighEnergyExcessPlugin(OctetAnalyzer* OA): OctetAnalyzerPlugin(OA,"highenergy") {
	for(Side s = EAST; s <= WEST; ++s) {
		qExcessr2[s] = registerCoreHist("Excessr2","Excess events radial distribution",10,0,80*80,s);
		qExcessr2[s]->setAxisTitle(X_DIRECTION,"event radius squared [mm^{2}]");
		qExcessTheta[s] = registerCoreHist("ExcessTheta","Excess events angular distribution",10,-M_PI,M_PI,s);
		qExcessTheta[s]->setAxisTitle(X_DIRECTION,"event position angle [radians]");
		qExcessSpectra[s] = registerCoreHist("ExcessE","Excess high energy event spectrum",72,800,8000,s);
		qExcessSpectra[s]->setAxisTitle(X_DIRECTION,"Energy [keV]");
		qExcessGamma[s] = registerCoreHist("ExcessGamma","Excess high energy gamma event spectrum",72,800,8000,s);
		qExcessGamma[s]->setAxisTitle(X_DIRECTION,"Energy [keV]");
	}
}

void HighEnergyExcessPlugin::fillCoreHists(ProcessedDataScanner& PDS, double weight) {
	Side s = PDS.fSide;
	if(!(s==EAST || s==WEST)) return;
	if(PDS.fPID == PID_SINGLE && PDS.fType == TYPE_IV_EVENT)
		qExcessGamma[s]->fillPoint->Fill(PDS.getErecon(),weight);
	if(PDS.fPID != PID_BETA) return;
	if(PDS.fType <= TYPE_I_EVENT)
		qExcessSpectra[s]->fillPoint->Fill(PDS.getErecon(),weight);
	if(PDS.fType <= TYPE_I_EVENT && PDS.getErecon()>1000) {
		qExcessr2[s]->fillPoint->Fill(PDS.radius2(s),weight);
		qExcessTheta[s]->fillPoint->Fill(atan2(PDS.wires[s][Y_DIRECTION].center,PDS.wires[s][X_DIRECTION].center),weight);
	}
}

void HighEnergyExcessPlugin::calculateResults() {
	if(myA->grouping < GROUP_OCTET) return;
	for(Side s = EAST; s <= WEST; ++s) {
		fitHighEnergyExcess(myA,qExcessSpectra[s],1000,2200);
		fitHighEnergyExcess(myA,qExcessSpectra[s],2200,7000);
		fitHighEnergyExcess(myA,qExcessGamma[s],200,1000);
		fitHighEnergyExcess(myA,qExcessGamma[s],1000,2200);
		fitHighEnergyExcess(myA,qExcessGamma[s],2200,7000);
	}
}

void HighEnergyExcessPlugin::makePlots() {
	if(myA->grouping == GROUP_RANGE) {
		drawQuadSides(qExcessGamma[EAST],qExcessGamma[WEST],true,"Energy");
		drawQuadSides(qExcessSpectra[EAST],qExcessSpectra[WEST],true,"Energy");
		drawQuadSides(qExcessr2[EAST],qExcessr2[WEST],true,"Positions");
		drawQuadSides(qExcessTheta[EAST],qExcessTheta[WEST],true,"Positions");
		
		std::vector<TH1*> hToPlot;
		for(Side s = EAST; s <= WEST; ++s) {
			TH1* hXS = myA->flipperSummedRate(qExcessSpectra[s], GV_OPEN);
			hXS->GetYaxis()->SetTitleOffset(1.45);
			hXS->Scale(1000000);
			hXS->GetYaxis()->SetTitle("event rate [uHz/keV]");
			hXS->SetTitle("Excess high energy beta events");
			hXS->SetMinimum(0);
			hXS->SetMaximum(30);
			hXS->SetLineStyle(1+s);
			hXS->SetLineColor(1);
			hToPlot.push_back(hXS);
		}
		drawSimulHistos(hToPlot);
		printCanvas("Energy/HighEnergyBetas");
	}
}

void HighEnergyExcessPlugin::compareMCtoData(AnalyzerPlugin* AP) {
	if(myA->grouping < GROUP_RANGE) return;
	
	// re-cast to correct type
	HighEnergyExcessPlugin& dat = *(HighEnergyExcessPlugin*)AP;
	
	printf("Comparing data/MC high energy excess events...\n");
	
	std::vector<TH1*> hToPlot;
	std::vector<TH1*> hToPlotR;
	for(Side s = EAST; s <= WEST; ++s) {
		TH1* hXS = myA->flipperSummedRate(dat.qExcessSpectra[s], GV_OPEN);
		hXS->SetLineStyle(1+s);
		hXS->GetYaxis()->SetTitleOffset(1.45);
		hXS->GetXaxis()->SetRangeUser(800,5000);
		hXS->Scale(1000000);
		hXS->GetYaxis()->SetTitle("event rate [uHz/keV]");
		hXS->SetTitle("Excess high energy beta events");
		hXS->SetLineColor(1);
		hXS->SetMinimum(-5);
		hXS->SetMaximum(25);
		hToPlot.push_back(hXS);
		
		TH1* hXSR = myA->flipperSummedRate(dat.qExcessr2[s], GV_OPEN);
		hXSR->SetLineStyle(1+s);
		hXSR->SetLineWidth(2);
		hXSR->GetYaxis()->SetTitleOffset(1.45);
		hXSR->Scale(1000000);
		hXSR->GetYaxis()->SetTitle("event rate [uHz/mm^{2}]");
		hXSR->SetTitle("Excess events radial distribution");
		hXSR->SetLineColor(1);
		hToPlotR.push_back(hXSR);

	}
	drawSimulHistos(hToPlot);
	TH1* hXSim = myA->flipperSummedRate(qExcessSpectra[EAST], GV_OPEN);
	hXSim->Scale(1000000);
	hXSim->SetMarkerStyle(33);
	hXSim->SetMarkerSize(0.75);
	hXSim->Draw("P Same");
	printCanvas("DataComparison/HighEnergyBetas");
	drawSimulHistos(hToPlotR);
	TH1* hXSimR = myA->flipperSummedRate(qExcessr2[EAST], GV_OPEN);
	hXSimR->Scale(1000000);
	hXSimR->SetMarkerStyle(33);
	hXSimR->SetMarkerSize(1.25);
	hXSimR->SetLineColor(0);
	hXSimR->SetLineWidth(0);
	hXSimR->Draw("P Same");
	for(int r = 10; r <= 70; r += 10) drawVLine(r*r, myA->defaultCanvas, 1, 3);
	printCanvas("DataComparison/HighEnergyRadial");
}


//-------------------------------------------------------------------------------------------------------------

NGBGAnalyzer::NGBGAnalyzer(OutputManager* pnt, const std::string& nm, const std::string& inflName):
OctetAnalyzer(pnt,nm,inflName) {
	addPlugin(myHEE = new HighEnergyExcessPlugin(this));
	addPlugin(myAsym = new AsymmetryPlugin(this));
}

//-------------------------------------------------------------//

void NGBGSpectra(std::string datname) {
	OutputManager OM("NGBG",getEnvSafe("UCNA_ANA_PLOTS")+"/NGBG/");
	G4toPMT g2p;
	g2p.addFile(getEnvSafe("G4OUTDIR")+"/"+datname+"/analyzed_*.root");
	g2p.setAFP(AFP_OFF);
	g2p.weightAsym = false;
	PMTCalibrator PCal(15668);
	g2p.setCalibrator(PCal);
	
	NGBGAnalyzer AH(&OM,datname);
	AH.loadSimData(g2p);
	
	AH.calculateResults();
	AH.makePlots();
	AH.write();
	AH.setWriteRoot(true);
}

