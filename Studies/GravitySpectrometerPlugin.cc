#include "GravitySpectrometerPlugin.hh"
#include "GraphUtils.hh"

GravitySpectrometerPlugin::GravitySpectrometerPlugin(OctetAnalyzer* OA): OctetAnalyzerPlugin(OA,"gravSpec") {
	myA->ignoreMissingHistos = true;
	qHeight = registerCoreHist("qHeight", "event height", 20, -60, 60);
	qHeight->setAxisTitle(X_DIRECTION,"y Position [mm]");
	TH2F hPositionsTemplate("hGrav","Event Height v Time",150,0,300,20,-60,60);
	qTime = registerCoreHist(hPositionsTemplate);
	qTime->setTimeScaling(false);
	qTime->setAxisTitle(X_DIRECTION,"time [s]");
	qTime->setAxisTitle(Y_DIRECTION,"y Position [mm]");
	myA->ignoreMissingHistos = false;
}

void GravitySpectrometerPlugin::fillCoreHists(ProcessedDataScanner& PDS, double weight) {
	Side s = PDS.fSide;
	if(PDS.fPID!=PID_BETA || !(s==EAST||s==WEST) || PDS.fType > TYPE_III_EVENT) return;
	qHeight->fillPoint->Fill(PDS.wires[s][Y_DIRECTION].center,PDS.physicsWeight);
	((TH2F*)qTime->fillPoint)->Fill(PDS.runClock[s],PDS.wires[s][Y_DIRECTION].center,PDS.physicsWeight);
}

void GravitySpectrometerPlugin::calculateResults() {
	
	TF1 lineFit("lineFit","[0]*(1+[1]*x)",-50,50);
	
	// time filling profile
	for(AFPState afp = AFP_OFF; afp <= AFP_ON; ++afp)
		timeProf[afp] = ((TH2F*)qTime->fgbg[afp]->h[GV_OPEN])->ProjectionX();
	// on vs off filling ratio
	hTimeRat = (TH1*)timeProf[AFP_OFF]->Clone();
	hTimeRat->Divide(timeProf[AFP_ON]);

	// on vs off height profile
	hAFPRat = (TH1*)qHeight->fgbg[AFP_OFF]->h[GV_OPEN]->Clone();
	hAFPRat->Divide(qHeight->fgbg[AFP_ON]->h[GV_OPEN]);
	hAFPRat->Fit(&lineFit,"QR");
	Stringmap m;
	m.insert("p0",lineFit.GetParameter(0));
	m.insert("d_p0",lineFit.GetParError(0));
	m.insert("p1",lineFit.GetParameter(1));
	m.insert("d_p1",lineFit.GetParError(1));
	myA->qOut.insert("height_Off_v_On", m);
	hAFPRat->GetYaxis()->SetRangeUser(lineFit.GetParameter(0)*0.97,lineFit.GetParameter(0)*1.03);
	
	// time evolution of height profile
	for(AFPState afp = AFP_OFF; afp <= AFP_ON; ++afp) {
		((TH2F*)qTime->fgbg[afp]->h[GV_OPEN])->RebinX(2);
		timeSlices[afp] = sliceTH2(*(TH2F*)qTime->fgbg[afp]->h[GV_OPEN], X_DIRECTION);
		timeEvol[afp] = new TGraphErrors((int)timeSlices[afp].size());
		for(unsigned int i=0; i<timeSlices[afp].size(); i++) {
			timeSlices[afp][i]->Divide(qHeight->fgbg[afp]->h[GV_OPEN]);
			timeSlices[afp][i]->Fit(&lineFit,"QR");
			double tm = ((TH2F*)qTime->fgbg[afp]->h[GV_OPEN])->GetXaxis()->GetBinCenter(i+1);
			timeEvol[afp]->SetPoint(i,tm,lineFit.GetParameter(1));
			timeEvol[afp]->SetPointError(i,0,lineFit.GetParError(1));
			Stringmap m2;
			m2.insert("afp",afpWords(afp));
			m2.insert("p0",lineFit.GetParameter(0));
			m2.insert("d_p0",lineFit.GetParError(0));
			m2.insert("p1",lineFit.GetParameter(1));
			m2.insert("d_p1",lineFit.GetParError(1));
			m2.insert("time",tm);
			myA->qOut.insert("height_timeslice", m2);
			m2.display();			
		}
	}
}

void GravitySpectrometerPlugin::makePlots() {
	hAFPRat->Draw();
	myA->printCanvas("Gravity/Height_Off_v_On");
	hTimeRat->GetYaxis()->SetRangeUser(1.2,1.6);
	hTimeRat->Draw();
	myA->printCanvas("Gravity/Time_Off_v_On");
	for(AFPState afp = AFP_OFF; afp <= AFP_ON; ++afp) {
		qHeight->fgbg[afp]->h[GV_OPEN]->Draw();
		myA->printCanvas("Gravity/HeightProf_"+afpWords(afp));
		qTime->fgbg[afp]->h[GV_OPEN]->Draw("COL");
		myA->printCanvas("Gravity/Pos_v_Time_"+afpWords(afp));
		timeProf[afp]->Draw();
		myA->printCanvas("Gravity/TimeProf_"+afpWords(afp));
		timeEvol[afp]->Draw("AP");
		myA->printCanvas("Gravity/TimeEvol_"+afpWords(afp));
	}
}
