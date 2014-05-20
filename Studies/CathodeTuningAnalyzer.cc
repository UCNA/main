#include "CathodeTuningAnalyzer.hh"
#include "GraphUtils.hh"
#include "GraphicsUtils.hh"
#include <TObjArray.h>
#include <TProfile.h>

//---------------------------------------------------

Stringmap cathseg2sm(const CathodeSeg& c) {
	Stringmap m;
	m.insert("side",sideSubst("%c",c.s));
	m.insert("plane",c.d==X_DIRECTION?"x":"y");
	m.insert("i",c.i);
	m.insert("position",c.pos);
	m.insert("n_exp",c.n_exp);
	m.insert("n_obs",c.n_obs);
	m.insert("dndx_lo",c.dndx_lo);
	m.insert("dndx_hi",c.dndx_hi);
	return m;
}

//---------------------------------------------------

CathodeGainPlugin::CathodeGainPlugin(RunAccumulator* RA): AnalyzerPlugin(RA,"cathodegain") {
	for(Side s = EAST; s <= WEST; ++s) {
		for(AxisDirection d = X_DIRECTION; d <= Y_DIRECTION; ++d) {
			for(unsigned int c=0; c < kMaxCathodes; c++) {
				TH2F hCathNormTemplate((std::string("hCathNorm_")+(d==X_DIRECTION?"x_":"y_")+itos(c)).c_str(),
									   "Normalized Cathode",601,-6.0,6.0,100,-0.05,0.4);
				cathNorm[s][d][c] = registerFGBGPair(hCathNormTemplate,AFP_OTHER,s);
				cathNorm[s][d][c]->setAxisTitle(X_DIRECTION,"normalized position");
				cathNorm[s][d][c]->setAxisTitle(Y_DIRECTION,"cathode charge fraction");
			}
			TProfile pTemplate((std::string("pCathGainfact_")+(d==X_DIRECTION?"x_":"y_")).c_str(),"Cathode Gain Factors",kMaxCathodes,-0.5,kMaxCathodes-0.5);
			prevGain[s][d] = registerFGBGPair(pTemplate,AFP_OTHER,s);
			prevGain[s][d]->doSubtraction = false;
		}
	}
}

void CathodeGainPlugin::fillCoreHists(ProcessedDataScanner& PDS, double weight) {
	const Side s = PDS.fSide;
	if(!(PDS.fPID == PID_BETA && (s==EAST||s==WEST) && PDS.fType == TYPE_0_EVENT)) return;
	double x = PDS.wires[s][X_DIRECTION].center;
	double y = PDS.wires[s][Y_DIRECTION].center;
	double eta_cc = PDS.ActiveCal->ccloud_eta[s]->eval(s,0,x,y,true);
	if(PDS.mwpcEnergy[s]>15.) return;	// avoid clipped events
	unsigned int n;
	float c;
	for(AxisDirection d = X_DIRECTION; d <= Y_DIRECTION; ++d) {
		PDS.ActiveCal->toLocal(s,d,PDS.wires[s][d].center,n,c);
		smassert(n<kMaxCathodes);
		for(unsigned int i=0; i<PDS.ActiveCal->nWires(s,d); i++) {
			double gi_gcc = PDS.ActiveCal->getCathCCloudGain(s,d,i);
			double c_enorm = PDS.ActiveCal->cathseg_energy_norm[s][d][i];
			((TH2F*)cathNorm[s][d][i]->h[currentGV])->Fill(c+(int(n)-int(i)), PDS.cathodes[s][d][i] / (c_enorm*eta_cc*PDS.mwpcEnergy[s]), weight);
			((TProfile*)prevGain[s][d]->h[currentGV])->Fill(i,gi_gcc,weight);
		}
	}
}

void CathodeGainPlugin::calculateResults() {
	
	TF1 fLowCX("fLowCX","pol2",-0.6,-0.4);
	TF1 fHighCX("fHiCX","pol2",0.4,0.6);
	TF1 fCathCenter("fCathCenter","gaus",-0.3,0.3);
	fCathCenter.SetLineColor(1);
	fLowCX.SetLineColor(0);
	fHighCX.SetLineColor(0);
	
	for(Side s = EAST; s <= WEST; ++s) {
		for(AxisDirection d = X_DIRECTION; d <= Y_DIRECTION; ++d) {
			for(unsigned int c=0; c<kMaxCathodes; c++) {
				
				// cathode/anode vs normalized position TH2
				TH2* hCathNorm =  (TH2*)cathNorm[s][d][c]->h[GV_OPEN];
				
				// fit each position slice with a gaussian
				// (better than TProfile for rejecting far-from-center events)
				TObjArray sls;
				hCathNorm->FitSlicesY(0, 0, -1, 0, "QNR", &sls);
				sls.SetOwner(kFALSE);
				for(int i=0; i<sls.GetEntriesFast(); i++) {
					myA->addObject(sls[i]);
					slicefits[s][d][c].push_back((TH1D*)sls[i]);
				}
				
				// fit gaussian centers vs normalized position, at crossover points and center
				slicefits[s][d][c][1]->Fit(&fLowCX,"QR");
				slicefits[s][d][c][1]->Fit(&fCathCenter,"QR+");
				slicefits[s][d][c][1]->Fit(&fHighCX,"QR+");
				
				// write fits output
				std::string fitnm = "cathnorm_fit_"+std::string(d==X_DIRECTION?"x":"y")+"_"+itos(c);
				AnaNumber AN(fitnm+"_prevgain");
				AN.s = s;
				
				AN.value = ((TProfile*)prevGain[s][d]->h[GV_OPEN])->GetBinContent(c);
				myA->uploadAnaNumber(AN, GV_OPEN, AFP_OTHER);
				
				AN.name = fitnm + "_mid";
				for(int n=0; n<3; n++) {
					AN.n = n;
					AN.value = fCathCenter.GetParameter(n);
					AN.err = fCathCenter.GetParameter(n);
					myA->uploadAnaNumber(AN, GV_OPEN, AFP_OTHER);
				}
				
				AN.name = fitnm + "_lo";
				for(int n=0; n<3; n++) {
					AN.n = n;
					AN.value = fLowCX.GetParameter(n);
					AN.err = fLowCX.GetParError(n);
					myA->uploadAnaNumber(AN, GV_OPEN, AFP_OTHER);
				}
				
				AN.name = fitnm + "_hi";
				for(int n=0; n<3; n++) {
					AN.n = n;
					AN.value = fHighCX.GetParameter(n);
					AN.err = fHighCX.GetParError(n);
					myA->uploadAnaNumber(AN, GV_OPEN, AFP_OTHER);
				}
			}
		}
	}
}

void CathodeGainPlugin::makePlots() {
	myA->defaultCanvas->SetLeftMargin(0.14);
	myA->defaultCanvas->SetRightMargin(0.04);
	for(Side s = EAST; s <= WEST; ++s) {
		for(AxisDirection d = X_DIRECTION; d <= Y_DIRECTION; ++d) {
			for(unsigned int c=0; c<kMaxCathodes; c++) {
				std::string cname = sideSubst("%c",s)+(d==X_DIRECTION?"x":"y")+itos(c+1);
				TH2* hCathNorm =  (TH2*)cathNorm[s][d][c]->h[GV_OPEN];
				hCathNorm->SetMinimum(0);
				hCathNorm->GetYaxis()->SetTitleOffset(1.55);
				hCathNorm->SetTitle(("Cathode "+cname).c_str());
				hCathNorm->GetXaxis()->SetRangeUser(-3, 3);
				hCathNorm->Draw("Col");
				if(slicefits[s][d][c].size()>=2) {
					slicefits[s][d][c][1]->SetLineColor(0);
					slicefits[s][d][c][1]->Draw("Same");
				}
				printCanvas("Cathodes/Center_"+cname);
			}
		}
	}
}

//---------------------------------------------------

CathodeTweakPlugin::CathodeTweakPlugin(RunAccumulator* RA): AnalyzerPlugin(RA,"wirechamber") {
	TH2F hPositionsTemplate("hPositions","Event Positions",200,-60,60,200,-60,60);
	hPositionsTemplate.GetXaxis()->SetTitle("x position [mm]");
	hPositionsTemplate.GetYaxis()->SetTitle("y position [mm]");
	TH2F hPositionsRawTemplate("hPositionsRaw","Event Raw Positions",200,-60,60,200,-60,60);
	hPositionsRawTemplate.GetXaxis()->SetTitle("x position [mm]");
	hPositionsRawTemplate.GetYaxis()->SetTitle("y position [mm]");
	for(Side s = EAST; s <= WEST; ++s) {
		hitPos[s] = registerFGBGPair(hPositionsTemplate,AFP_OTHER,s);
		hitPosRaw[s] = registerFGBGPair(hPositionsRawTemplate,AFP_OTHER,s);
		for(AxisDirection d = X_DIRECTION; d <= Y_DIRECTION; ++d) {
			for(unsigned int c=0; c < kMaxCathodes; c++) {
				TH2F hCathHitposTemplate((std::string("hCathHitpos_")+(d==X_DIRECTION?"x_":"y_")+itos(c)).c_str(),
										 "Cathode Hit Positions",64,-0.5,0.5,1000,0,1000);
				cathHitpos[s][d][c] = registerFGBGPair(hCathHitposTemplate,AFP_OTHER,s);
				cathHitpos[s][d][c]->setAxisTitle(X_DIRECTION,"Normalized Raw Position");
				cathHitpos[s][d][c]->setAxisTitle(Y_DIRECTION,"Energy [keV]");
			}
		}
	}
}

void CathodeTweakPlugin::fillCoreHists(ProcessedDataScanner& PDS, double weight) {
	const Side s = PDS.fSide;
	if(!(PDS.fType == TYPE_0_EVENT && PDS.fPID == PID_BETA && (s==EAST||s==WEST))) return;
	for(AxisDirection d = X_DIRECTION; d <= Y_DIRECTION; ++d) {
		unsigned int n;
		float c;
		PDS.ActiveCal->toLocal(s,d,PDS.wires[s][d].rawCenter,n,c);
		smassert(n<kMaxCathodes);
		((TH2F*)cathHitpos[s][d][n]->h[currentGV])->Fill(c,PDS.scints[s].energy.x,weight);
	}
	((TH2F*)(hitPos[s]->h[currentGV]))->Fill(PDS.wires[s][X_DIRECTION].center,PDS.wires[s][Y_DIRECTION].center,weight);
	((TH2F*)(hitPosRaw[s]->h[currentGV]))->Fill(PDS.wires[s][X_DIRECTION].rawCenter,PDS.wires[s][Y_DIRECTION].rawCenter,weight);
}

void CathodeTweakPlugin::makePlots() {
	
	PMTCalibrator PCal(16000);
	
	for(Side s = EAST; s <= WEST; ++s) {
		
		myA->defaultCanvas->SetRightMargin(0.14);
		myA->defaultCanvas->SetLeftMargin(0.10);
		
		// normalize to rate per cm^2
		double ascale = 10.*10.;
		
		TH2F* hPos = (TH2F*)hitPos[s]->h[GV_OPEN]->Clone();
		hPos->Scale(ascale/(hPos->GetXaxis()->GetBinWidth(1)*hPos->GetYaxis()->GetBinWidth(1)*myA->getTotalTime(AFP_OTHER, GV_OPEN)[BOTH]));
		hPos->SetMaximum(0.2);
		hPos->Draw("COL Z");
		PCal.drawWires(s, X_DIRECTION, myA->defaultCanvas, 0, X_DIRECTION);
		PCal.drawWires(s, Y_DIRECTION, myA->defaultCanvas, 0, Y_DIRECTION);
		drawCircle(50, 1, 2);
		printCanvas(sideSubst("hPos_%c",s));
		
		
		TH2F* hRaw = (TH2F*)hitPosRaw[s]->h[GV_OPEN]->Clone();
		hRaw->Scale(ascale/(hRaw->GetXaxis()->GetBinWidth(1)*hRaw->GetYaxis()->GetBinWidth(1)*myA->getTotalTime(AFP_OTHER, GV_OPEN)[BOTH]));
		hRaw->SetMaximum(0.2);
		hRaw->Draw("COL Z");
		PCal.drawWires(s, X_DIRECTION, myA->defaultCanvas, 0, X_DIRECTION);
		PCal.drawWires(s, Y_DIRECTION, myA->defaultCanvas, 0, Y_DIRECTION);
		drawCircle(50, 1, 2);
		printCanvas(sideSubst("hPosRaw_%c",s));
		
		// profiles
		myA->defaultCanvas->SetRightMargin(0.04);
		myA->defaultCanvas->SetLeftMargin(0.12);
		for(AxisDirection d = X_DIRECTION; d <= Y_DIRECTION; ++d) {
			
			TH1* pPos = d==X_DIRECTION?hPos->ProjectionX():hPos->ProjectionY();
			pPos->GetYaxis()->SetTitle("event rate [Hz/mm]");
			pPos->GetYaxis()->SetTitleOffset(1.5);
			pPos->Scale(1./ascale);
			pPos->GetYaxis()->SetRangeUser(0,0.34);
			pPos->Draw();
			pPos->Draw("HIST SAME");
			PCal.drawWires(s, d, myA->defaultCanvas,1);
			printCanvas(sideSubst("hPos_%c",s)+(d==X_DIRECTION?"X":"Y"));
			
			TH1* pRaw = d==X_DIRECTION?hRaw->ProjectionX():hRaw->ProjectionY();
			pRaw->GetYaxis()->SetTitle("event rate [Hz/mm]");
			pRaw->GetYaxis()->SetTitleOffset(1.5);
			pRaw->Scale(1./ascale);
			pRaw->GetYaxis()->SetRangeUser(0,0.34);
			pRaw->Draw();
			pRaw->Draw("HIST SAME");
			PCal.drawWires(s, d, myA->defaultCanvas,1);
			printCanvas(sideSubst("hPosRaw_%c",s)+(d==X_DIRECTION?"X":"Y"));
			
			delete pPos;
			delete pRaw;
		}
		
		delete hPos;
		delete hRaw;
	}
}

//-----------------------------------------------------------

CathodeTuningAnalyzer::CathodeTuningAnalyzer(OutputManager* pnt, const std::string& nm, const std::string& inflName):
OctetAnalyzer(pnt,nm,inflName) {
	addPlugin(myCG = new CathodeGainPlugin(this));
	addPlugin(myCT = new CathodeTweakPlugin(this));
	//addPlugin(myGS = new GravitySpectrometerPlugin(this));
}

SimCathodeTuningAnalyzer::SimCathodeTuningAnalyzer(OutputManager* pnt, const std::string& nm, const std::string& inflName):
RunAccumulator(pnt,nm,inflName) {
	addPlugin(myCT = new CathodeTweakPlugin(this));
}

//-----------------------------------------------------------

void processCathTweak(CathodeTweakPlugin& CTDat, CathodeTweakPlugin& CTSim) {
	
	RunNum r1 = CTDat.myA->runCounts.counts.begin()->first;
	RunNum r2 = CTDat.myA->runCounts.counts.rbegin()->first;
	OutputManager OM("MWPCCal",getEnvSafe("UCNA_ANA_PLOTS")+"/WirechamberCal/"+itos(r1)+"-"+itos(r2)+"/");
	OM.qOut.transfer(CTDat.myA->qOut,"runcal");
	
	// construct fourier series TF1 fit function
	const unsigned int nterms = 3;
	std::string fser = "[0]";
	for(unsigned int n=1; n<nterms; n++)
		fser += " + ["+itos(2*n-1)+"]*sin(2*pi*"+itos(n)+"*x) + ["+itos(2*n)+"]*cos(2*pi*"+itos(n)+"*x)";
	printf("Fitter '%s'\n",fser.c_str());
	TF1 fFourier("fFourier",fser.c_str(),-0.5,0.5);
	
	for(Side s = EAST; s <= WEST; ++s) {
		for(AxisDirection d = X_DIRECTION; d <= Y_DIRECTION; ++d) {
			for(unsigned int c=0; c<kMaxCathodes; c++) {
				
				TH2* hCathHitposDat =  (TH2*)CTDat.cathHitpos[s][d][c]->h[GV_OPEN];
				TH2* hCathHitposSim =  (TH2*)CTSim.cathHitpos[s][d][c]->h[GV_OPEN];
				
				// cathode event fractions, simulation vs data
				CathodeSeg cs;
				cs.s = s;
				cs.d = d;
				cs.i = c;
				cs.n_obs = hCathHitposDat->Integral();
				cs.n_exp = hCathHitposSim->Integral();
				unsigned int bxm0 = hCathHitposDat->GetXaxis()->FindBin(-0.499);
				unsigned int bxm1 = hCathHitposDat->GetXaxis()->FindBin(-0.4);
				unsigned int bxp0 = hCathHitposDat->GetXaxis()->FindBin(0.4);
				unsigned int bxp1 = hCathHitposDat->GetXaxis()->FindBin(0.499);
				unsigned int ny = hCathHitposDat->GetNbinsY();
				cs.dndx_lo = hCathHitposDat->Integral(bxm0,bxm1,1,ny)/(hCathHitposDat->GetXaxis()->GetBinUpEdge(bxm1)-hCathHitposDat->GetXaxis()->GetBinLowEdge(bxm0));
				cs.dndx_hi = hCathHitposDat->Integral(bxp0,bxp1,1,ny)/(hCathHitposDat->GetXaxis()->GetBinUpEdge(bxp1)-hCathHitposDat->GetXaxis()->GetBinLowEdge(bxp0));
				OM.qOut.insert("cathseg",cathseg2sm(cs));
				
				if(cs.n_exp<10000) continue;
				
				// determine variable energy bin size re-binning based on available statistics
				std::vector<TH1F*> hEnDats = sliceTH2(*hCathHitposDat,Y_DIRECTION);
				std::vector<TH1F*> hEnSims = sliceTH2(*hCathHitposSim,Y_DIRECTION);
				std::vector<float> enCounts;
				for(unsigned int i=0; i<hEnDats.size(); i++)
					enCounts.push_back(hEnDats[i]->Integral());
				std::vector<unsigned int> part = equipartition(enCounts, 80);
				
				// fit each energy range with Fourier series
				std::vector<TH1*> hToPlot;
				for(unsigned int i = 0; i < part.size()-1; i++) {
					// accumulate over energy range
					TH1F* hEnDat = hEnDats[part[i]];
					TH1F* hEnSim = hEnSims[part[i]];
					float eavg = hEnDat->Integral()*hCathHitposDat->GetYaxis()->GetBinCenter(part[i]+1);
					for(unsigned int j=part[i]+1; j<part[i+1]; j++) {
						hEnDat->Add(hEnDats[j]);
						hEnSim->Add(hEnSims[j]);
						eavg += hEnDats[j]->Integral()*hCathHitposDat->GetYaxis()->GetBinCenter(j+1);
						delete(hEnDats[j]);
						delete(hEnSims[j]);
					}
					eavg /= hEnDat->Integral();
					//printf("Accumulating %u-%u (%.1fkeV) of %i\n",part[i],part[i+1],eavg,(int)hEnDats.size());
					//hEnDat->Rebin();
					//hEnSim->Rebin();
					hEnSim->Scale(hEnDat->Integral()/hEnSim->Integral());
					hEnDat->Divide(hEnSim);
					fixNaNs(hEnDat);
					hEnDat->SetLineColor(i+1);
					OM.addObject(hEnDat);
					OM.addObject(hEnSim);
					
					// fourier components fit
					if(c==1)
						fFourier.SetRange(-0.1,0.5);
					else if(c==kMaxCathodes-2)
						fFourier.SetRange(-0.5,0.1);
					else
						fFourier.SetRange(-0.5,0.5);
					fFourier.SetParameter(0,1.0);
					for(unsigned int n=1; n<=2*nterms; n++)
						fFourier.SetParameter(n,0.);
					fFourier.FixParameter(3,0.);
					fFourier.SetLineColor(i+1);
					hEnDat->Fit(&fFourier,"QR");
					hEnDat->SetMinimum(0);
					hEnDat->SetMaximum(1.5);
					hToPlot.push_back(hEnDat);
					
					// record results
					std::vector<double> terms;
					std::vector<double> dterms;
					double c0 = fFourier.GetParameter(0);
					for(unsigned int n=0; n<2*nterms-1; n++) {
						terms.push_back(fFourier.GetParameter(n)/c0);
						dterms.push_back(fFourier.GetParError(n)/c0);
					}
					Stringmap ff;
					ff.insert("side",sideSubst("%c",s));
					ff.insert("plane",d==X_DIRECTION?"x":"y");
					ff.insert("cathode",c);
					ff.insert("elo",hCathHitposDat->GetYaxis()->GetBinLowEdge(part[i]+1));
					ff.insert("ehi",hCathHitposDat->GetYaxis()->GetBinUpEdge(part[i+1]));
					ff.insert("eavg",eavg);
					ff.insert("terms",vtos(terms));
					ff.insert("dterms",vtos(dterms));
					//ff.display();
					OM.qOut.insert("hitdist",ff);
				}
				drawSimulHistos(hToPlot);
				OM.printCanvas(sideSubst("CathTweak/Positions_%c",s)+(d==X_DIRECTION?"x":"y")+itos(c));
			}
		}
	}
	
	OM.write();
	OM.setWriteRoot(true);
}

//void processWirechamberCal(RunNum r0, RunNum r1, unsigned int nrings) {
//	std::string basePath = getEnvSafe("UCNA_ANA_PLOTS")+"/PositionMaps/";
//	OutputManager OM("NameUnused",basePath);
//	std::string readname = itos(r0)+"-"+itos(r1)+"_"+itos(nrings);
//	CathodeTuningAnalyzer WCdat(&OM, "NameUnused", basePath+"/Xenon_"+readname+"/Xenon_"+readname);
//	CathodeTuningAnalyzer WCsim(&OM, "NameUnused", basePath+"/SimXe_"+readname+"/SimXe_"+readname);
//	processWirechamberCal(WCdat,WCsim);
//}
