#include "PositionsPlugin.hh"
#include "GraphicsUtils.hh"

PositionsPlugin::PositionsPlugin(OctetAnalyzer* OA): OctetAnalyzerPlugin(OA,"position"), offSects(5,45.0) {
	for(unsigned int m=0; m<offSects.nSectors(); m++) {
		for(AxisDirection d = X_DIRECTION; d <= Y_DIRECTION; ++d) {
			poff[d].push_back(registerFGBGPair((d==X_DIRECTION?"pOff_X_":"pOff_Y_")+itos(m), "Type I Position Offsets",50,-25, 25));
			poff[d].back()->setAxisTitle(X_DIRECTION,"E-W offset [mm]");
		}
	}
	
	for(Side s = EAST; s <= WEST; ++s) {
		for(EventType t=TYPE_0_EVENT; t<=TYPE_III_EVENT; ++t) {
			TH2F hPositionsTemplate(("hPos_Type_"+itos(t)).c_str(),
									("Type "+itosRN(t)+" Positions").c_str(),
									200,-65,65,200,-65,65);
			qPositions[s][t] = registerCoreHist(hPositionsTemplate,s);
			qPositions[s][t]->setAxisTitle(X_DIRECTION,"x Position [mm]");
			qPositions[s][t]->setAxisTitle(Y_DIRECTION,"y Position [mm]");
			myA->ignoreMissingHistos = true;
			qRadius2[s][t] = registerCoreHist("hR2_Type_"+itos(t), "Type "+itosRN(t)+" radius^{2}", 60, 0, 60*60, s);
			qRadius2[s][t]->setAxisTitle(X_DIRECTION,"radius^{2} [mm^{2}]");
			myA->ignoreMissingHistos = false;
		}
	}
}

void PositionsPlugin::fillCoreHists(ProcessedDataScanner& PDS, double weight) {
	Side s = PDS.fSide;
	if(PDS.fPID!=PID_BETA || !(s==EAST||s==WEST)) return;
	if(PDS.fType <= TYPE_III_EVENT) {
		((TH2F*)qPositions[s][PDS.fType]->fillPoint)->Fill(PDS.wires[s][X_DIRECTION].center,PDS.wires[s][Y_DIRECTION].center,weight);
		qRadius2[s][PDS.fType]->fillPoint->Fill(PDS.radius2(s),weight);
	}
	if(PDS.fType != TYPE_I_EVENT) return;
	double x = 0.5*(PDS.wires[EAST][X_DIRECTION].center+PDS.wires[WEST][X_DIRECTION].center);
	double y = 0.5*(PDS.wires[EAST][Y_DIRECTION].center+PDS.wires[WEST][Y_DIRECTION].center);
	unsigned int m = offSects.sector(x,y);
	if(m>=offSects.nSectors()) return;
	for(AxisDirection d = X_DIRECTION; d <= Y_DIRECTION; ++d)
		poff[d][m]->h[currentGV]->Fill(PDS.wires[EAST][d].center-PDS.wires[WEST][d].center,weight);
}

void PositionsPlugin::calculateResults() {
	float x,y;
	TF1 gausFit("gasufit","gaus",-25,25);
	gausFit.SetLineColor(2);
	for(unsigned int m=0; m<offSects.nSectors(); m++) {
		Stringmap odat;
		offSects.sectorCenter(m,x,y);
		odat.insert("m",itos(m));
		odat.insert("x",x);
		odat.insert("y",y);
		for(AxisDirection d = X_DIRECTION; d <= Y_DIRECTION; ++d) {
			std::string axname = (d==X_DIRECTION?"x":"y");
			poff[d][m]->h[GV_OPEN]->Fit(&gausFit,"QR");
			odat.insert("d"+axname,gausFit.GetParameter(1));
			odat.insert("d_d"+axname,gausFit.GetParError(1));
			odat.insert("w"+axname,gausFit.GetParameter(2));
			odat.insert("d_w"+axname,gausFit.GetParError(2));
		}
		myA->qOut.insert("posOffset",odat);
	}
	for(EventType tp = TYPE_0_EVENT; tp <= TYPE_III_EVENT; ++tp) {
		hSuperSumPos[tp] = calculateSuperSum("Event_Radius_SuperSum_Tp"+itos(tp),qRadius2[EAST][tp],qRadius2[WEST][tp],GV_OPEN,true);
		hSuperSumPos[tp]->Scale(1000);
		hSuperSumPos[tp]->SetTitle(("Type "+itosRN(tp)+" events radius distribution").c_str());
		hSuperSumPos[tp]->GetYaxis()->SetTitle("event rate [mHz/mm^{2}]");
		hSuperSumPos[tp]->GetYaxis()->SetTitleOffset(tp?1.7:1.5);
	}
}

void PositionsPlugin::makePlots() {
	if(myA->grouping < GROUP_RANGE) return;
	
	myA->defaultCanvas->cd();
	
	for(EventType tp=TYPE_0_EVENT; tp<=TYPE_III_EVENT; ++tp) {
		std::vector<TH1*> hToPlot;
		for(Side s = EAST; s <= WEST; ++s) {
			myA->defaultCanvas->SetRightMargin(0.08);
			myA->defaultCanvas->SetLeftMargin(0.10);
			
			qPositions[s][tp]->setDrawRange(0,false);
			drawQuad(qPositions[s][tp],"Positions","COL");
			drawQuad(qRadius2[s][tp],"Positions","E L");
			
			TH1* hBgPos = myA->flipperSummedRate(qPositions[s][tp], GV_CLOSED, false);
			hBgPos->SetTitle(("Type "+itosRN(tp)+" background positions").c_str());
			hBgPos->Draw(tp==TYPE_0_EVENT?"scat=0.5":"");
			printCanvas(sideSubst("Positions/BGPos_%c_Type_",s)+itos(tp));
			delete hBgPos;
			
			TH1* hBgRad = myA->flipperSummedRate(qRadius2[s][tp], GV_CLOSED);
			hBgRad->SetLineStyle(1+s);
			hBgRad->SetLineColor(1);
			hBgRad->GetYaxis()->SetTitleOffset(1.5);
			hBgRad->SetTitle(("Type "+itosRN(tp)+" background radius^{2}").c_str());
			hBgRad->Scale(1000000);
			hBgRad->GetYaxis()->SetTitle("event rate [uHz/mm^{2}]");
			hToPlot.push_back(hBgRad);
		}
		myA->defaultCanvas->SetRightMargin(0.04);
		myA->defaultCanvas->SetLeftMargin(0.14);
		drawSimulHistos(hToPlot,"HIST");
		for(int r = 10; r <= 50; r += 10) drawVLine(r*r, myA->defaultCanvas, 1, 3);
		printCanvas("Positions/BGRad_Type_"+itos(tp));
	}
}

void PositionsPlugin::compareMCtoData(AnalyzerPlugin* AP) {
	if(myA->grouping < GROUP_RANGE) return;
	
	// re-cast to correct type
	PositionsPlugin& dat = *(PositionsPlugin*)AP;
	
	myA->defaultCanvas->SetRightMargin(0.04);
	myA->defaultCanvas->SetLeftMargin(0.14);
		
	/*
	for(unsigned int t=TYPE_0_EVENT; t<=TYPE_III_EVENT; t++) {
		std::vector<TH1*> hToPlot;
		for(Side s = EAST; s <= WEST; ++s) {
			for(AFPState afp = AFP_OFF; afp <= AFP_ON; ++afp) {
				TH1* hR2sim = qRadius2[s][t]->fgbg[afp]->h[GV_OPEN];
				hR2sim->SetMarkerColor(2+2*s);
				hR2sim->SetMarkerStyle(22+4*afp);
				hR2sim->SetMarkerSize(0.25);
				hToPlot.push_back(hR2sim);
				TH1* hR2dat = dat.qRadius2[s][t]->fgbg[afp]->h[GV_OPEN];
				hR2dat->SetMarkerColor(2+2*s);
				hR2dat->SetMarkerStyle(20+4*afp);
				hR2dat->SetMarkerSize(0.25);
				hR2dat->Scale(hR2sim->Integral()/hR2dat->Integral());
				hToPlot.push_back(hR2dat);
			}
		}
		drawSimulHistos(hToPlot,"HIST P E");
		printCanvas("DataComparison/Radius2_Type_"+itos(t));
	}
	*/
	
	for(EventType tp = TYPE_0_EVENT; tp <= TYPE_III_EVENT; ++tp) {
		hSuperSumPos[tp]->SetMarkerSize(0.5);
		hSuperSumPos[tp]->SetMarkerStyle(33);
		hSuperSumPos[tp]->Scale(dat.hSuperSumPos[tp]->Integral()/hSuperSumPos[tp]->Integral());
		
		dat.hSuperSumPos[tp]->Draw();
		hSuperSumPos[tp]->Draw("P SAME");
		for(int r = 10; r <= 70; r += 10) drawVLine(r*r, myA->defaultCanvas, 1, 3);
		printCanvas("DataComparison/Radius2_Type_"+itos(tp));
	}
}
