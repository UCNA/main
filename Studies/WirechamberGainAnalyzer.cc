#include "WirechamberGainAnalyzer.hh"

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

WirechamberGainAnalyzer::WirechamberGainAnalyzer(RunAccumulator* RA): AnalyzerPlugin(RA,"wirechambergain") {
	for(Side s = EAST; s <= WEST; ++s) {
		anodeCal[s] = registerFGBGPair("AnodeCal","Anode Calibration Events",50, 0, 8, AFP_OTHER, s);
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

void WirechamberGainAnalyzer::fillCoreHists(ProcessedDataScanner& PDS, double weight) {
	const Side s = PDS.fSide;
	if(!(PDS.fType == TYPE_0_EVENT && PDS.fPID == PID_BETA && (s==EAST||s==WEST))) return;
	if(PDS.mwpcs[s].anode < 1000) { // avoid most massively clipped events
		unsigned int n;
		float c;
		for(AxisDirection d = X_DIRECTION; d <= Y_DIRECTION; ++d) {
			PDS.ActiveCal->toLocal(s,d,PDS.wires[s][d].center,n,c);
			assert(n<kMaxCathodes);
			((TH2F*)cathNorm[s][d][n]->h[currentGV])->Fill(c,PDS.cathodes[s][d][n]/PDS.mwpcs[s].anode,weight);
		}
	}
	if(PDS.passesPositionCut(s) && PDS.getEtrue()>225)
		anodeCal[s]->h[currentGV]->Fill(PDS.mwpcEnergy[s]/PDS.ActiveCal->wirechamberGainCorr(s,PDS.runClock[s]),weight);
}

void WirechamberGainAnalyzer::calculateResults() {
	for(Side s = EAST; s <= WEST; ++s) {
		TF1 fLandau("landauFit","landau",0,15);
		fLandau.SetLineColor(2+2*s);
		int fiterr = anodeCal[s]->h[GV_OPEN]->Fit(&fLandau,"Q");
		Stringmap m;
		m.insert("side",ctos(sideNames(s)));
		m.insert("fiterr",itos(fiterr));
		m.insert("height",fLandau.GetParameter(0));
		m.insert("d_height",fLandau.GetParError(0));
		m.insert("mpv",fLandau.GetParameter(1));
		m.insert("d_mpv",fLandau.GetParError(1));
		m.insert("sigma",fLandau.GetParameter(2));
		m.insert("d_sigma",fLandau.GetParError(2));
		myA->qOut.insert("anodeCalFit",m);
	}
}

void WirechamberGainAnalyzer::makePlots() {
	//drawQuadSides(anodeCal[EAST], anodeCal[WEST], true, "AnodeCal");
	//TODO
}

void WirechamberGainAnalyzer::compareMCtoData(AnalyzerPlugin* AP) {
	// re-cast to correct type
	WirechamberGainAnalyzer& dat = *(WirechamberGainAnalyzer*)AP;
}


