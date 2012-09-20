#include "WirechamberGainAnalyzer.hh"
#include "GraphUtils.hh"
#include "GraphicsUtils.hh"
#include <TH2F.h>
#include <TH3F.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TGraphAsymmErrors.h>

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
	for(Side s = EAST; s <= WEST; ++s) {
		norm23[s] = registerFGBGPair("norm23", "Type II/III Normalized MWPC", 50, -1, 4, AFP_OTHER, s);
		norm23[s]->setAxisTitle(X_DIRECTION, "Normalized MWPC signal");
	}
	for(EventType tp = TYPE_0_EVENT; tp <= TYPE_III_EVENT; ++tp) {
		TH2F hAnodeTemplate(("hAnodeEnergy_Tp"+itos(tp)).c_str(),"MWPC Energy",10,0,1000,100,0,10);
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
			for(unsigned int i=0; i<hSlices[s][tp].size(); i++)
				delete hSlices[s][tp][i];
		}
	}
}

void AnodeGainAnalyzer::fillCoreHists(ProcessedDataScanner& PDS, double weight) {
	const Side s = PDS.fSide;
	if(!(PDS.fPID == PID_BETA && (s==EAST||s==WEST))) return;
	((TProfile*)anodeGaincorr->h[currentGV])->Fill(s,PDS.ActiveCal->wirechamberGainCorr(s,PDS.runClock[BOTH]),weight);
	if(PDS.passesPositionCut(s) && PDS.fType <= TYPE_III_EVENT) {
		((TH2F*)(anodeCal[s][PDS.fType]->h[currentGV]))->Fill(PDS.getEnergy(),PDS.mwpcEnergy[s],weight);
		if(PDS.fType >= TYPE_II_EVENT)
			norm23[s]->h[currentGV]->Fill(WirechamberCalibrator::normMWPC(s,PDS.getEnergy(),PDS.mwpcEnergy[s]),weight);
	}
}

void AnodeGainAnalyzer::calculateResults() {
	// fit anode spectrum vs. scintillator energy
	TF1 fLandau("landauFit","landau",0,10);
	fLandau.SetParLimits(1,0.1,10.);
	fLandau.SetParLimits(2,0.20,5.0);
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
			hSlices[s][tp] = sliceTH2(*hAnode,X_DIRECTION);
			gAnode[s][tp] = new TGraphErrors(hSlices[s][tp].size());
			TGraphErrors gAvg(hSlices[s][tp].size());
			for(unsigned int i=0; i<hSlices[s][tp].size(); i++) {
				double e0 = hAnode->GetBinCenter(i+1);
				double c = hSlices[s][tp][i]->GetBinCenter(hSlices[s][tp][i]->GetMaximumBin());
				double mx = hSlices[s][tp][i]->GetMaximum();
				fLandau.SetParameter(2,0.23*c);
				fLandau.SetParameter(0,c*mx/fLandau.GetParameter(2));
				fLandau.SetParameter(1,c);
				fLandau.SetRange(0,2*c);
				printf("\n---- Anode %s Type %i E=%i : estimated h=%.2g at MPV=%.2f (sigma=%.2f) ----\n",
					   sideWords(s),tp,(int)e0,fLandau.GetParameter(0),c,fLandau.GetParameter(2));
				hSlices[s][tp][i]->Fit(&fLandau,"RBNME");
				Stringmap m;
				m.insert("side",ctos(sideNames(s)));
				m.insert("type",itos(tp));
				m.insert("energy",e0);
				m.insert("height",fLandau.GetParameter(0));
				m.insert("d_height",fLandau.GetParError(0));
				m.insert("mpv",fLandau.GetParameter(1));
				m.insert("d_mpv",fLandau.GetParError(1));
				m.insert("sigma",fLandau.GetParameter(2));
				m.insert("d_sigma",fLandau.GetParError(2));
				myA->qOut.insert("anodeCalFit",m);
				gAnode[s][tp]->SetPoint(i,e0,fLandau.GetParameter(1));
				gAnode[s][tp]->SetPointError(i,hAnode->GetBinWidth(i+1)*0.5,fLandau.GetParameter(2));
				gAvg.SetPoint(i,e0,fLandau.GetParameter(1));
				gAvg.SetPointError(i,0,fLandau.GetParError(1));
			}
			
			// average anode MPV energy
			double emin = tp==TYPE_0_EVENT?400:200;
			double emax = tp==TYPE_0_EVENT?800:500;
			fAvg.SetRange(emin,emax);
			gAvg.Fit(&fAvg,"QR");
			Stringmap m;
			m.insert("side",ctos(sideNames(s)));
			m.insert("type",itos(tp));
			m.insert("emin",emin);
			m.insert("emax",emax);
			m.insert("avg",fAvg.GetParameter(0));
			m.insert("d_avg",fAvg.GetParError(0));
			myA->qOut.insert("anodeCalAvg",m);
		}
	}
}

void AnodeGainAnalyzer::makePlots() {
	for(Side s = EAST; s <= WEST; ++s) {
		for(EventType tp = TYPE_0_EVENT; tp <= TYPE_III_EVENT; ++tp) {
			if(hSlices[s][tp].size()) {
				std::vector<TH1*> hToPlot;
				for(unsigned int i=0; i<hSlices[s][tp].size(); i++)
					hToPlot.push_back(hSlices[s][tp][i]);
				drawSimulHistos(hToPlot);
				printCanvas(sideSubst("AnodeCal/hAnode_%c_",s)+itos(tp));
			}
			if(!gAnode[s][tp]) continue;
			gAnode[s][tp]->Draw("AP");
			gAnode[s][tp]->SetTitle((sideSubst("%s Type ",s)+itos(tp)+" Events").c_str());
			gAnode[s][tp]->GetXaxis()->SetTitle("Scintillator Energy [keV]");
			gAnode[s][tp]->GetYaxis()->SetTitle("Primary Side MWPC Energy [keV]");
			gAnode[s][tp]->GetYaxis()->SetRangeUser(0.0,10.0);
			gAnode[s][tp]->GetXaxis()->SetRangeUser(0.0,800.0);
			gAnode[s][tp]->Draw("AP");
			printCanvas(sideSubst("AnodeCal/EMWPC_%c_",s)+itos(tp));
		}
	}
}

void AnodeGainAnalyzer::compareMCtoData(AnalyzerPlugin* AP) {
	// re-cast to correct type
	AnodeGainAnalyzer& dat = *(AnodeGainAnalyzer*)AP;
	
	for(Side s = EAST; s <= WEST; ++s) {
		std::vector<TH1*> hToPlot;
		norm23[s]->h[GV_OPEN]->SetLineColor(4);
		dat.norm23[s]->h[GV_OPEN]->SetLineColor(2);
		norm23[s]->h[GV_OPEN]->Scale(dat.norm23[s]->h[GV_OPEN]->Integral()/norm23[s]->h[GV_OPEN]->Integral());
		hToPlot.push_back(norm23[s]->h[GV_OPEN]);
		hToPlot.push_back(dat.norm23[s]->h[GV_OPEN]);
		drawSimulHistos(hToPlot);
		drawVLine(0,myA->defaultCanvas);
		printCanvas(sideSubst("DataComparison/NormMWPC_%c",s));
		
		for(EventType tp = TYPE_0_EVENT; tp <= TYPE_III_EVENT; ++tp) {
			if(!gAnode[s][tp] || !dat.gAnode[s][tp]) continue;
			dat.gAnode[s][tp]->SetLineColor(2);
			dat.gAnode[s][tp]->Draw("AP");
			dat.gAnode[s][tp]->SetTitle((sideSubst("%s Type ",s)+itos(tp)+" Events").c_str());
			dat.gAnode[s][tp]->GetXaxis()->SetTitle("Scintillator Energy [keV]");
			dat.gAnode[s][tp]->GetYaxis()->SetTitle("Primary Side MWPC Energy [keV]");
			dat.gAnode[s][tp]->GetYaxis()->SetRangeUser(0.0,10.0);
			dat.gAnode[s][tp]->GetXaxis()->SetRangeUser(0.0,800.0);
			dat.gAnode[s][tp]->Draw("AP");
			gAnode[s][tp]->SetLineColor(4);
			gAnode[s][tp]->Draw("P");
			printCanvas(sideSubst("DataComparison/EMWPC_%c_",s)+itos(tp));
		}
		
		// Type 0 upper energy anode slices
		unsigned int i0 = 1;
		for(unsigned int i=i0+1; i<hSlices[s][TYPE_0_EVENT].size(); i++) {
			hSlices[s][TYPE_0_EVENT][i0]->Add(hSlices[s][TYPE_0_EVENT][i]);
			assert(i<dat.hSlices[s][TYPE_0_EVENT].size());
			dat.hSlices[s][TYPE_0_EVENT][i0]->Add(dat.hSlices[s][TYPE_0_EVENT][i]);
		}
		drawHistoPair(dat.hSlices[s][TYPE_0_EVENT][i0], hSlices[s][TYPE_0_EVENT][i0]);
		printCanvas(sideSubst("DataComparison/AnodeSpec_%c",s));
	}
}


//---------------------------------------------------


WirechamberSimTypeID::WirechamberSimTypeID(RunAccumulator* RA): AnalyzerPlugin(RA,"wirechamberTypeID") {
	for(EventType tp = TYPE_0_EVENT; tp <= TYPE_III_EVENT; ++tp) {
		TH3F hAnodeTemplate(("hAnodeTypeID_Tp"+itos(tp)).c_str(),"MWPC Energy",10,0,1000,50,0,10,50,0,10);
		TH2F hAnodeNormTemplate(("hAnodeNormCoords_Tp"+itos(tp)).c_str(),"MWPC Normalized",10,0,1000,50,-1.0,4.0);
		for(Side s = EAST; s <= WEST; ++s) {
			if(tp<TYPE_II_EVENT) {
				anodeTypeID[s][tp] = anodeNormCoords[s][tp] = NULL;
				continue;
			}
			anodeNormCoords[s][tp] = registerFGBGPair(hAnodeNormTemplate, AFP_OTHER, s);
			anodeNormCoords[s][tp]->setAxisTitle(X_DIRECTION, "Scintillator Energy [keV]");
			anodeNormCoords[s][tp]->setAxisTitle(Y_DIRECTION, "MWPC Normalized to II/III Cut");
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
	if(SDS.fPID == PID_BETA && (s==EAST||s==WEST) && PDS.passesPositionCut(s) && (tp==TYPE_II_EVENT||tp==TYPE_III_EVENT)) {
		float Escint = SDS.getEnergy();
		float Emwpc = SDS.mwpcEnergy[s];
		float nMWPC = WirechamberCalibrator::sep23Cut(s,Escint);
		((TH3F*)(anodeTypeID[s][rtp]->h[currentGV]))->Fill(Escint,Emwpc,SDS.mwpcEnergy[otherSide(s)],weight);
		((TH2F*)(anodeNormCoords[s][rtp]->h[currentGV]))->Fill(Escint,(Emwpc-nMWPC)/nMWPC,weight);
	}
}

Double_t crudeEffic(const Double_t* x, const Double_t* params) {
	if(*x<0)
		return params[0]+(0.5-params[0])*exp(*x/params[1]);
	return params[2]+(0.5-params[2])*exp(-*x/params[3]);
}

void WirechamberSimTypeID::make23SepInfo(OutputManager& OM) {
	OM.defaultCanvas->cd();
	OM.defaultCanvas->SetLeftMargin(0.12);
	
	// find optimal division point on each side vs. energy
	std::vector<TGraph*> gOptDiv;
	for(Side s = EAST; s <= WEST; ++s) {
		// separate by energy bins
		TH3F* hTpII = (TH3F*)(anodeTypeID[s][TYPE_II_EVENT]->h[GV_OPEN]);
		unsigned int nx = hTpII->GetNbinsX();
		gOptDiv.push_back(new TGraph(nx));
		gOptDiv[s]->SetMarkerColor(2+2*s);
		std::vector<TH2F*> v2 = sliceTH3(*hTpII,X_DIRECTION);
		std::vector<TH2F*> v3 = sliceTH3(*(TH3F*)(anodeTypeID[s][TYPE_III_EVENT]->h[GV_OPEN]),X_DIRECTION);
		for(unsigned int ne = 1; ne <= nx; ne++) {
			// energy window under consideration
			double e0 = hTpII->GetXaxis()->GetBinLowEdge(ne);
			double e1 = hTpII->GetXaxis()->GetBinUpEdge(ne);
			std::string eRange = itos(e0)+"-"+itos(e1)+" keV";
			
			// plot 2D distributions
			v2[ne]->SetMarkerColor(2);
			v2[ne]->SetTitle(("Type II/III Backscatters, "+eRange).c_str());
			v2[ne]->Draw();
			v3[ne]->SetMarkerColor(4);
			v3[ne]->Draw("Same");
			OM.printCanvas(sideSubst("PrimVSec_%c_",s)+itos(ne));
			
			if(false) {
				//////////////////////////
				// study separation efficiency by angle
				//////////////////////////
				int nAngles = 50;
				std::vector<TH1*> aProj2, aProj3;
				TGraph gAngleSep(nAngles);
				TGraph gDivPt(nAngles);
				double amin = 0;
				double omin = 1.0;
				double divmin = 0;
				for(int a = 0; a < nAngles; a++) {
					double th = a*M_PI/nAngles-M_PI/2;
					TH1* ap2 = projectTH2(*v2[ne],50,cos(th),sin(th));
					ap2->SetLineColor(2);
					ap2->SetTitle(("Optimized Type II/III separation, "+eRange).c_str());
					ap2->GetXaxis()->SetTitle("MWPC Energy [keV]");
					aProj2.push_back(ap2);
					TH1* ap3 = projectTH2(*v3[ne],50,cos(th),sin(th));
					ap3->SetLineColor(4);
					ap3->SetTitle(("Optimized Type II/III separation, "+eRange).c_str());
					ap3->GetXaxis()->SetTitle("MWPC Energy [keV]");
					aProj3.push_back(ap3);
					double o,xdiv;
					if(ap2->GetMaximumBin() <= ap3->GetMaximumBin())
						histoverlap(*ap2, *ap3, o, xdiv);
					else
						histoverlap(*ap3, *ap2, o, xdiv);
					double cTot = ap2->Integral()+ap3->Integral();
					gAngleSep.SetPoint(a,th,o/cTot);
					gDivPt.SetPoint(a,th,xdiv);
					//if(o/cTot < omin) {
					if(a==nAngles/2) {
						omin = o/cTot;
						amin = a;
						divmin = xdiv;
					}
				}
				gAngleSep.Draw("A*");
				gAngleSep.GetXaxis()->SetRangeUser(-M_PI/2.,M_PI/2);
				gAngleSep.GetXaxis()->SetTitle("Projection angle [radians]");
				gAngleSep.GetYaxis()->SetRangeUser(0.,0.5);
				gAngleSep.GetYaxis()->SetTitle("Misidentified event fraction");
				gAngleSep.GetYaxis()->SetTitleOffset(1.5);
				gAngleSep.SetTitle(("Type II/III Separation Cut Optimization, "+eRange).c_str());
				gAngleSep.Draw("A*");
				OM.printCanvas(sideSubst("OptCut_%c_",s)+itos(ne));
				
				gOptDiv[s]->SetPoint(ne-1, 0.5*(e0+e1), divmin);
				
				std::vector<TH1*> hToPlot;
				hToPlot.push_back(aProj2[amin]);
				hToPlot.push_back(aProj3[amin]);
				drawSimulHistos(hToPlot);
				drawVLine(divmin,OM.defaultCanvas);
				OM.printCanvas(sideSubst("OptSeparation_%c_",s)+itos(ne));
				for(int a = 0; a < nAngles; a++) {
					delete aProj2[a];
					delete aProj3[a];
				}
			} else {
				//////////////////////////
				// use only primary side MWPC energy
				//////////////////////////
				
				TH1D* ap2 = v2[ne]->ProjectionX();
				ap2->SetLineColor(2);
				ap2->SetTitle(("Optimized Type II/III separation, "+eRange).c_str());
				ap2->GetXaxis()->SetTitle("MWPC Energy [keV]");
				TH1D* ap3 = v3[ne]->ProjectionX();
				ap3->SetLineColor(4);
				ap3->SetTitle(("Optimized Type II/III separation, "+eRange).c_str());
				ap3->GetXaxis()->SetTitle("MWPC Energy [keV]");
				
				// optimize separation point
				TH1* hSep = histsep(*ap2,*ap3);
				hSep->Scale(1.0/(ap2->Integral()+ap3->Integral()));
				hSep->SetTitle(("Type II/III MWPC cut optimization, "+eRange).c_str());
				hSep->GetYaxis()->SetRangeUser(0.,1.);
				hSep->GetYaxis()->SetTitle("Misidentified Fraction");
				hSep->GetYaxis()->SetTitleOffset(1.2);
				// fit minimum region with cubic for refined minimum
				int bsep = hSep->GetMinimumBin();
				double c = hSep->GetBinCenter(bsep);
				TF1 minFit("minFit","pol3",c-1,c+1);
				hSep->Fit(&minFit,"WR");
				double p1 = 3*minFit.GetParameter(3);
				double p2 = 2*minFit.GetParameter(2);
				double p3 = 1*minFit.GetParameter(1);
				double del = sqrt(p2*p2-4*p1*p3);
				double x1 = (-p2+del)/(2*p1);
				double x2 = (-p2-del)/(2*p1);
				double xdiv = minFit.Eval(x1)<minFit.Eval(x2)?x1:x2;
				if(1 < xdiv && xdiv < 9)
					gOptDiv[s]->SetPoint(ne-1, 0.5*(e0+e1), xdiv);
				else
					gOptDiv[s]->SetPoint(ne-1,0,0);
				
				hSep->Draw();
				drawVLine(xdiv,OM.defaultCanvas);
				OM.printCanvas(sideSubst("SepErr_%c_",s)+itos(ne));
				
				std::vector<TH1*> hToPlot;
				hToPlot.push_back(ap2);
				hToPlot.push_back(ap3);
				drawSimulHistos(hToPlot);
				drawVLine(xdiv,OM.defaultCanvas);
				OM.printCanvas(sideSubst("OptSeparation_%c_",s)+itos(ne));
				
				delete ap2;
				delete ap3;
				delete hSep;
			}
			
		}
	}
	
	TF1 cutFit("cutFit","[0]+[1]*exp(-x/[2])",25,600);
	//cutFit.SetParameter(0,3.0);
	//cutFit.SetParameter(1,4.0);
	cutFit.SetParameter(0,5.0);
	cutFit.SetParLimits(0, 1.0, 8.0);
	cutFit.SetParameter(1,7.0);
	cutFit.SetParLimits(1,1.0,15.0);
	//cutFit.SetParameter(2,150);
	cutFit.SetParameter(2,300);
	cutFit.SetParLimits(2,75,1000);
	TGraph* gCombo = combine_graphs(gOptDiv);
	printf("\n\n///////////// Separation Cut Fit ////////////////\n\n");
	gCombo->Fit(&cutFit,"RBME");
	printf("%.2f + %.2f*exp(-Escint/%.1f)\n",cutFit.GetParameter(0),cutFit.GetParameter(1),cutFit.GetParameter(2));
	gCombo->Draw("A*");
	gCombo->GetXaxis()->SetRangeUser(0.,700.);
	gCombo->GetXaxis()->SetTitle("Scintillator Energy [keV]");
	gCombo->GetYaxis()->SetRangeUser(0.,10.0);
	gCombo->GetYaxis()->SetTitle("Optimum MWPC Cut [keV]");
	gCombo->SetTitle("Type II/III Separation Cut");
	gCombo->Draw("AP");
	gOptDiv[EAST]->Draw("*");
	gOptDiv[WEST]->Draw("*");
	OM.printCanvas("OptCut");
	
	
	// normalized MWPC Type II/III probability extraction
	TF1 efficfit("efficfit",&crudeEffic,-1.0,4.0,4);
	std::vector<TGraph*> efficParams[2];
	for(Side s = EAST; s <= WEST; ++s) {
		TH2F* hTpII = (TH2F*)(anodeNormCoords[s][TYPE_II_EVENT]->h[GV_OPEN]);
		unsigned int nx = hTpII->GetNbinsX();
		std::vector<TH1F*> v2 = sliceTH2(*hTpII,X_DIRECTION,true);
		std::vector<TH1F*> v3 = sliceTH2(*(TH2F*)(anodeNormCoords[s][TYPE_III_EVENT]->h[GV_OPEN]),X_DIRECTION,true);
		for(int n=0; n<efficfit.GetNpar(); n++)
			efficParams[s].push_back(new TGraph(nx));
		for(unsigned int ne = 1; ne <= nx; ne++) {
			// energy window under consideration
			double e0 = hTpII->GetXaxis()->GetBinLowEdge(ne);
			double e1 = hTpII->GetXaxis()->GetBinUpEdge(ne);
			if(0.5*(e0+e1)>700) continue;
			std::string eRange = itos(e0)+"-"+itos(e1)+" keV";
			printf("--- %s ---\n",eRange.c_str());
			
			// backscatter counts and mis-ID'd events
			Stringmap m;
			int b0 = v2[ne]->FindBin(0);
			int nb = v2[ne]->GetNbinsX();
			m.insert("side",sideWords(s));
			m.insert("e0",e0);
			m.insert("e1",e1);
			m.insert("II",v2[ne]->Integral(0,b0));
			m.insert("IIx",v2[ne]->Integral(b0,nb+1));
			m.insert("IIIx",v3[ne]->Integral(0,b0));
			m.insert("III",v3[ne]->Integral(b0,nb+1));
			OM.qOut.insert("23counts",m);
			
			// sum and plot
			v2[ne]->SetLineColor(2);
			v3[ne]->SetLineColor(4);
			TH1F* hBoth = (TH1F*)v2[ne]->Clone("sumIIandIII");
			hBoth->Add(v3[ne]);
			hBoth->SetLineColor(1);
			hBoth->SetTitle(("Normalized MWPC for "+eRange).c_str());
			hBoth->Draw();
			v2[ne]->Draw("Same");
			v3[ne]->Draw("Same");
			OM.printCanvas(sideSubst("NormMWPC_%c_",s)+itos(ne));
			
			// divide histograms
			TGraphAsymmErrors* gPIII = new TGraphAsymmErrors(hBoth->GetNbinsX());
			hBoth->GetSumw2()->Set(NULL);
			v3[ne]->GetSumw2()->Set(NULL);
			gPIII->BayesDivide(v3[ne],hBoth);
			
			
			efficfit.SetParameter(0,0.1);
			efficfit.SetParLimits(0,0.0,0.5);
			//efficfit.SetParameter(1,0.25);
			//efficfit.SetParLimits(1,0.01,10.0);
			efficfit.FixParameter(1,0.15);
			efficfit.SetParameter(2,0.9);
			efficfit.SetParLimits(2,0.5,1.0);
			//efficfit.SetParameter(3,0.2);
			//efficfit.SetParLimits(3,0.01,10.0);
			efficfit.FixParameter(3,0.20);
			gPIII->Fit(&efficfit,"RBME");
			
			gPIII->Draw("AP");
			gPIII->SetTitle(("Type III probability, "+eRange).c_str());
			gPIII->GetXaxis()->SetTitle("MWPC Normalized to II/III Cut");
			gPIII->GetYaxis()->SetTitle("Probability of being Type III");
			gPIII->GetYaxis()->SetRangeUser(0.,1.);
			gPIII->GetYaxis()->SetTitleOffset(1.2);
			gPIII->Draw("AP");
			OM.printCanvas(sideSubst("ProbIII_%c_",s)+itos(ne));
			
			for(int n=0; n<efficfit.GetNpar(); n++)
				efficParams[s][n]->SetPoint(ne, 0.5*(e0+e1), efficfit.GetParameter(n));
			
			delete gPIII;
			delete hBoth;
		}
	}
	
	for(int n=0; n<efficfit.GetNpar(); n++) {
		std::vector<TGraph*> sgraphs;
		for(Side s = EAST; s <= WEST; ++s) {
			efficParams[s][n]->SetMarkerColor(2+2*s);
			sgraphs.push_back(efficParams[s][n]);
		}
		TGraph* gCombo = combine_graphs(sgraphs);
		
		if(n==0) {
			TF1 fPol("fPol","pol2",25,625);
			gCombo->Fit(&fPol,"R");
			printf("Low Asympt: %.3f+%.3g*Escint+%.3g*Escint*Escint\n",fPol.GetParameter(0),fPol.GetParameter(1),fPol.GetParameter(2));
		}
		if(n==2) {
			TF1 expFit("expFit","[0]-[1]*exp(-x/[2])",25,625);
			expFit.SetParameter(0,0.9);
			expFit.SetParameter(1,0.2);
			expFit.SetParameter(2,100);
			gCombo->Fit(&expFit,"RBME");
			printf("High Asympt: %.3f-%.3f*exp(-Escint/%.1f)\n",expFit.GetParameter(0),expFit.GetParameter(1),expFit.GetParameter(2));
		}
		
		gCombo->Draw("AP");
		gCombo->SetTitle("Fit Parameter Value");
		gCombo->GetXaxis()->SetTitle("Scintillator Energy [keV]");
		gCombo->GetXaxis()->SetRangeUser(0.,700.);
		gCombo->GetYaxis()->SetTitle("Fit parameter value");
		gCombo->GetYaxis()->SetTitleOffset(1.2);
		gCombo->Draw("AP");
		for(Side s = EAST; s <= WEST; ++s)
			sgraphs[s]->Draw("*");
		OM.printCanvas("ProbFitParam_"+itos(n));
	}
}

