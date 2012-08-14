#include "WirechamberGainAnalyzer.hh"
#include "GraphUtils.hh"
#include <TH2F.h>
#include <TH3F.h>
#include <TProfile.h>

Stringmap cathseg2sm(const CathodeSeg& c) {
	Stringmap m;
	m.insert("side",sideSubst("%c",c.s));
	m.insert("plane",c.d==X_DIRECTION?"x":"y");
	m.insert("i",c.i);
	m.insert("max",c.max.x);
	m.insert("d_max",c.max.err);
	m.insert("height",c.height.x);
	m.insert("d_height",c.height.err);
	m.insert("width",c.width.x);
	m.insert("d_width",c.width.err);
	m.insert("center",c.center.x);
	m.insert("d_center",c.center.x);
	m.insert("position",c.pos);
	m.insert("fill_frac",c.fill_frac);
	return m;
}

CathodeSeg sm2cathseg(const Stringmap& m) {
	CathodeSeg c;
	c.s=(m.getDefault("side","E")=="E")?EAST:WEST;
	c.d=m.getDefault("plane","x")=="x"?X_DIRECTION:Y_DIRECTION;
	c.i=m.getDefault("i",0);
	c.pos=m.getDefault("position",0);
	c.height=float_err(m.getDefault("height",0),m.getDefault("d_height",0));
	c.width=float_err(m.getDefault("width",0),m.getDefault("d_width",0));
	c.center=float_err(m.getDefault("center",0),m.getDefault("d_center",0));
	return c;
}

//---------------------------------------------------

CathodeGainAnalyzer::CathodeGainAnalyzer(RunAccumulator* RA): AnalyzerPlugin(RA,"cathodegain") {
	for(Side s = EAST; s <= WEST; ++s) {
		for(AxisDirection d = X_DIRECTION; d <= Y_DIRECTION; ++d) {
			for(unsigned int c=0; c < kMaxCathodes; c++) {
				TH2F hCathNormTemplate((std::string("hCathNorm_")+(d==X_DIRECTION?"x_":"y_")+itos(c)).c_str(),
									   "Normalized Cathode",51,-0.5,0.5,100,0,5);
				cathNorm[s][d][c] = registerFGBGPair(hCathNormTemplate,AFP_OTHER,s);
				cathNorm[s][d][c]->setAxisTitle(X_DIRECTION,"Normalized Position");
				cathNorm[s][d][c]->setAxisTitle(Y_DIRECTION,"Normalized Cathode");
			}
		}
	}
}

void CathodeGainAnalyzer::fillCoreHists(ProcessedDataScanner& PDS, double weight) {
	const Side s = PDS.fSide;
	if(!(PDS.fPID == PID_BETA && (s==EAST||s==WEST))) return;
	if(PDS.fType == TYPE_0_EVENT && PDS.mwpcs[s].anode < 1000) { // avoid most massively clipped events
		unsigned int n;
		float c;
		for(AxisDirection d = X_DIRECTION; d <= Y_DIRECTION; ++d) {
			PDS.ActiveCal->toLocal(s,d,PDS.wires[s][d].center,n,c);
			assert(n<kMaxCathodes);
			((TH2F*)cathNorm[s][d][n]->h[currentGV])->Fill(c,PDS.cathodes[s][d][n]/PDS.mwpcs[s].anode,weight);
		}
	}
}

//---------------------------------------------------

AnodeGainAnalyzer::AnodeGainAnalyzer(RunAccumulator* RA): AnalyzerPlugin(RA,"wirechambergain") {
	TProfile anodeGaincorrTemplate("pAnodeGaincorr","Anode Gain Correction",2,-0.5,1.5);
	anodeGaincorr = registerFGBGPair(anodeGaincorrTemplate, AFP_OTHER, BOTH);
	anodeGaincorr->doSubtraction = false;
	anodeGaincorr->setAxisTitle(X_DIRECTION, "Side");
	anodeGaincorr->setAxisTitle(Y_DIRECTION, "Anode Gain Correction Factor");
	for(EventType tp = TYPE_0_EVENT; tp <= TYPE_III_EVENT; ++tp) {
		TH2F hAnodeTemplate(("hAnodeEnergy_Tp"+itos(tp)).c_str(),"MWPC Energy",10,0,1000,50,0,10);
		for(Side s = EAST; s <= WEST; ++s) {
			gAnode[s][tp] = NULL;
			anodeCal[s][tp] = registerFGBGPair(hAnodeTemplate, AFP_OTHER, s);
			anodeCal[s][tp]->setAxisTitle(X_DIRECTION, "Scintillator Energy [keV]");
			anodeCal[s][tp]->setAxisTitle(Y_DIRECTION, "Primary Side MWPC Energy [keV]");
		}
	}
}

AnodeGainAnalyzer::~AnodeGainAnalyzer() {
	for(Side s = EAST; s <= WEST; ++s) {
		for(EventType tp = TYPE_0_EVENT; tp <= TYPE_III_EVENT; ++tp) {
			if(gAnode[s][tp]) {
				delete gAnode[s][tp];
				gAnode[s][tp] = NULL;
			}
		}
	}
}

void AnodeGainAnalyzer::fillCoreHists(ProcessedDataScanner& PDS, double weight) {
	const Side s = PDS.fSide;
	if(!(PDS.fPID == PID_BETA && (s==EAST||s==WEST))) return;
	((TProfile*)anodeGaincorr->h[currentGV])->Fill(s,PDS.ActiveCal->wirechamberGainCorr(s,PDS.runClock[BOTH]),weight);
	if(PDS.passesPositionCut(s) && PDS.fType <= TYPE_III_EVENT)
		((TH2F*)(anodeCal[s][PDS.fType]->h[currentGV]))->Fill(PDS.getEnergy(),PDS.mwpcEnergy[s],weight);
}

void AnodeGainAnalyzer::calculateResults() {
	// fit anode spectrum vs. scintillator energy
	TF1 fLandau("landauFit","landau",0,10);
	TF1 fAvg("avgFit","pol0",400,800);
	for(Side s = EAST; s <= WEST; ++s) {
		
		// average gain correction
		Stringmap mGC;
		mGC.insert("side",ctos(sideNames(s)));
		mGC.insert("avg",anodeGaincorr->h[GV_OPEN]->GetBinContent(s+1));
		mGC.insert("d_avg",anodeGaincorr->h[GV_OPEN]->GetBinError(s+1));
		myA->qOut.insert("anodeGaincorr",mGC);
		
		// fit each scintillator energy slice and record results
		fLandau.SetLineColor(2+2*s);
		for(EventType tp = TYPE_0_EVENT; tp <= TYPE_III_EVENT; ++tp) {
			TH2F* hAnode = (TH2F*)(anodeCal[s][tp]->h[GV_OPEN]);
			std::vector<TH1D*> fslices = replaceFitSlicesY(hAnode, &fLandau);
			int nEnergyBins = hAnode->GetNbinsX();
			assert(fslices.size()==4);
			for(int i = 0; i < 3; i++)
				assert(fslices[i]->GetNbinsX()==nEnergyBins);
			gAnode[s][tp] = new TGraphErrors(nEnergyBins);
			for(int i=0; i<nEnergyBins; i++) {
				double e0 = hAnode->GetBinCenter(i+1);
				Stringmap m;
				m.insert("side",ctos(sideNames(s)));
				m.insert("type",itos(tp));
				m.insert("energy",e0);
				m.insert("height",fslices[0]->GetBinContent(i+1));
				m.insert("d_height",fslices[0]->GetBinError(i+1));
				m.insert("mpv",fslices[1]->GetBinContent(i+1));
				m.insert("d_mpv",fslices[1]->GetBinError(i+1));
				m.insert("sigma",fslices[2]->GetBinContent(i+1));
				m.insert("d_sigma",fslices[2]->GetBinError(i+1));
				myA->qOut.insert("anodeCalFit",m);
				gAnode[s][tp]->SetPoint(i,e0,fslices[1]->GetBinContent(i+1));
				gAnode[s][tp]->SetPointError(i,hAnode->GetBinWidth(i+1)*0.5,fslices[2]->GetBinContent(i+1));
			}
			
			// average anode MPV energy
			fslices[1]->Fit(&fAvg,"QR");
			Stringmap m;
			m.insert("side",ctos(sideNames(s)));
			m.insert("type",itos(tp));
			m.insert("avg",fAvg.GetParameter(0));
			m.insert("d_avg",fAvg.GetParError(0));
			myA->qOut.insert("anodeCalAvg",m);

			for(unsigned int i=0; i<fslices.size(); ++i)
				delete fslices[i];
		}
	}
}

void AnodeGainAnalyzer::makePlots() {
	for(Side s = EAST; s <= WEST; ++s) {
		for(EventType tp = TYPE_0_EVENT; tp <= TYPE_III_EVENT; ++tp) {
			if(!gAnode[s][tp]) continue;
			gAnode[s][tp]->Draw("AP");
			gAnode[s][tp]->SetTitle((sideSubst("%s Type ",s)+itos(tp)+" Events").c_str());
			gAnode[s][tp]->GetXaxis()->SetTitle("Scintillator Energy [keV]");
			gAnode[s][tp]->GetYaxis()->SetTitle("Primary Side MWPC Energy [keV]");
			gAnode[s][tp]->GetYaxis()->SetRangeUser(0.0,10.0);
			gAnode[s][tp]->Draw("AP");
			printCanvas(sideSubst("AnodeCal/EMWPC_%c_",s)+itos(tp));
		}
	}
}

void AnodeGainAnalyzer::compareMCtoData(AnalyzerPlugin* AP) {
	// re-cast to correct type
	AnodeGainAnalyzer& dat = *(AnodeGainAnalyzer*)AP;
	
	for(Side s = EAST; s <= WEST; ++s) {
		for(EventType tp = TYPE_0_EVENT; tp <= TYPE_III_EVENT; ++tp) {
			if(!gAnode[s][tp] || !dat.gAnode[s][tp]) continue;
			dat.gAnode[s][tp]->SetLineColor(2);
			dat.gAnode[s][tp]->Draw("AP");
			dat.gAnode[s][tp]->GetXaxis()->SetTitle("Scintillator Energy [keV]");
			dat.gAnode[s][tp]->GetYaxis()->SetTitle("Primary Side MWPC Energy [keV]");
			dat.gAnode[s][tp]->GetYaxis()->SetRangeUser(0.0,10.0);
			dat.gAnode[s][tp]->Draw("AP");
			gAnode[s][tp]->SetLineColor(4);
			gAnode[s][tp]->Draw("P");
			printCanvas(sideSubst("DataComparison/EMWPC_%c_",s)+itos(tp));
		}
	}
}


//---------------------------------------------------


WirechamberSimTypeID::WirechamberSimTypeID(RunAccumulator* RA): AnalyzerPlugin(RA,"wirechamberTypeID") {
	for(EventType tp = TYPE_0_EVENT; tp <= TYPE_III_EVENT; ++tp) {
		TH3F hAnodeTemplate(("hAnodeTypeID_Tp"+itos(tp)).c_str(),"MWPC Energy",10,0,1000,50,0,10,50,0,10);
		for(Side s = EAST; s <= WEST; ++s) {
			if(tp<TYPE_II_EVENT) { anodeTypeID[s][tp] = NULL; continue; }
			anodeTypeID[s][tp] = registerFGBGPair(hAnodeTemplate, AFP_OTHER, s);
			anodeTypeID[s][tp]->setAxisTitle(X_DIRECTION, "Scintillator Energy [keV]");
			anodeTypeID[s][tp]->setAxisTitle(Y_DIRECTION, "Primary Side MWPC energy [keV]");
			anodeTypeID[s][tp]->setAxisTitle(Z_DIRECTION, "Secondary Side MWPC energy [keV]");
		}
	}
}

void WirechamberSimTypeID::fillCoreHists(ProcessedDataScanner& PDS, double weight) {
	assert(PDS.isSimulated());
	Sim2PMT& SDS = *(Sim2PMT*)&PDS;
	const Side s = SDS.fSide; // should be primary side
	EventType tp = SDS.fType;
	EventType rtp = (s!=SDS.primSide)?TYPE_II_EVENT:TYPE_III_EVENT;	// "real" event type
	if(SDS.fPID == PID_BETA && (s==EAST||s==WEST) && PDS.passesPositionCut(s) && (tp==TYPE_II_EVENT||tp==TYPE_III_EVENT))
		((TH3F*)(anodeTypeID[s][rtp]->h[currentGV]))->Fill(SDS.getEnergy(),SDS.mwpcEnergy[s],SDS.mwpcEnergy[otherSide(s)],weight);
}
