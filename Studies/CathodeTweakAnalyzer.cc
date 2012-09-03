#include "CathodeTweakAnalyzer.hh"
#include "GraphUtils.hh"
#include "GraphicsUtils.hh"
#include <TObjArray.h>

CathodeTweakAnalyzer::CathodeTweakAnalyzer(RunAccumulator* RA): AnalyzerPlugin(RA,"wirechamber") {
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

void CathodeTweakAnalyzer::fillCoreHists(ProcessedDataScanner& PDS, double weight) {
	const Side s = PDS.fSide;
	if(!(PDS.fType == TYPE_0_EVENT && PDS.fPID == PID_BETA && (s==EAST||s==WEST))) return;
	for(AxisDirection d = X_DIRECTION; d <= Y_DIRECTION; ++d) {
		unsigned int n;
		float c;
		PDS.ActiveCal->toLocal(s,d,PDS.wires[s][d].rawCenter,n,c);
		assert(n<kMaxCathodes);
		((TH2F*)cathHitpos[s][d][n]->h[currentGV])->Fill(c,PDS.scints[s].energy.x,weight);
	}
	((TH2F*)(hitPos[s]->h[currentGV]))->Fill(PDS.wires[s][X_DIRECTION].center,PDS.wires[s][Y_DIRECTION].center,weight);
	((TH2F*)(hitPosRaw[s]->h[currentGV]))->Fill(PDS.wires[s][X_DIRECTION].rawCenter,PDS.wires[s][Y_DIRECTION].rawCenter,weight);
}

void CathodeTweakAnalyzer::makePlots() {
	for(Side s = EAST; s <= WEST; ++s) {
		TH2F* hPos = (TH2F*)hitPos[s]->h[GV_OPEN]->Clone();
		hPos->Scale(10*M_PI*52*52/(hPos->GetXaxis()->GetBinWidth(1)*hPos->GetYaxis()->GetBinWidth(1)*hPos->Integral()));
		hPos->SetMaximum(14);
		hPos->Draw("COL Z");
		printCanvas(sideSubst("hPos_%c",s));
		delete hPos;
		
		TH2F* hRaw = (TH2F*)hitPosRaw[s]->h[GV_OPEN]->Clone();
		hRaw->Scale(10*M_PI*52*52/(hRaw->GetXaxis()->GetBinWidth(1)*hRaw->GetYaxis()->GetBinWidth(1)*hRaw->Integral()));
		hRaw->SetMaximum(14);
		hRaw->Draw("COL Z");
		printCanvas(sideSubst("hPosRaw_%c",s));
		delete hRaw;
	}
}

//-----------------------------------------------------------

MWPCTuningAnalyzer::MWPCTuningAnalyzer(OutputManager* pnt, const std::string& nm, const std::string& inflName):
RunAccumulator(pnt,nm,inflName) {
	addPlugin(myCG = new CathodeGainAnalyzer(this));
	addPlugin(myCT = new CathodeTweakAnalyzer(this));
}

//-----------------------------------------------------------

void processWirechamberCal(MWPCTuningAnalyzer& WCdat, MWPCTuningAnalyzer& WCsim) {
	
	RunNum r1 = WCdat.runCounts.counts.begin()->first;
	RunNum r2 = WCdat.runCounts.counts.rbegin()->first;
	OutputManager OM("MWPCCal",getEnvSafe("UCNA_ANA_PLOTS")+"/WirechamberCal/"+itos(r1)+"-"+itos(r2)+"/");
	OM.qOut.transfer(WCdat.qOut,"runcal");
	TF1 fCathCenter("fCathCenter","gaus",-0.5,0.5);
	fCathCenter.SetLineColor(2);
	
	const unsigned int nterms = 3;
	std::string fser = "[0]";
	for(unsigned int n=1; n<nterms; n++)
		fser += " + ["+itos(2*n-1)+"]*sin(2*pi*"+itos(n)+"*x) + ["+itos(2*n)+"]*cos(2*pi*"+itos(n)+"*x)";
	printf("Fitter '%s'\n",fser.c_str());
	TF1 fFourier("fFourier",fser.c_str(),-0.5,0.5);
	
	for(Side s = EAST; s <= WEST; ++s) {
		for(AxisDirection d = X_DIRECTION; d <= Y_DIRECTION; ++d) {
			for(unsigned int c=0; c<kMaxCathodes; c++) {
				////////////////////
				// cathode normalization plots
				////////////////////
				
				TH2* hCathNorm =  (TH2*)WCdat.myCG->cathNorm[s][d][c]->h[GV_OPEN];
				TObjArray sls;
				hCathNorm->FitSlicesY(0, 0, -1, 0, "QNR", &sls);
				sls.SetOwner(kFALSE);
				std::vector<TH1D*> slicefits;
				for(int i=0; i<sls.GetEntriesFast(); i++) {
					OM.addObject(sls[i]);
					slicefits.push_back((TH1D*)sls[i]);
				}
				slicefits[1]->Fit(&fCathCenter,"QR");
				CathodeSeg cs;
				cs.s = s;
				cs.d = d;
				cs.i = c;
				cs.height = float_err(fCathCenter.GetParameter(0),fCathCenter.GetParError(0));
				cs.center = float_err(fCathCenter.GetParameter(1),fCathCenter.GetParError(1));
				cs.width = float_err(fCathCenter.GetParameter(2),fCathCenter.GetParError(2));
				int bm = slicefits[1]->GetMaximumBin();
				cs.max = float_err(slicefits[1]->GetBinContent(bm),slicefits[1]->GetBinError(bm));
				hCathNorm->Draw("Col");
				slicefits[1]->Draw("Same");
				OM.printCanvas(sideSubst("Cathodes/Center_%c",s)+(d==X_DIRECTION?"x":"y")+itos(c));
				
				
				////////////////////
				// cathode positioning distributions by energy
				////////////////////
				
				TH2* hCathHitposDat =  (TH2*)WCdat.myCT->cathHitpos[s][d][c]->h[GV_OPEN];
				TH2* hCathHitposSim =  (TH2*)WCsim.myCT->cathHitpos[s][d][c]->h[GV_OPEN];
				
				// cathode event fractions
				double nDat = hCathHitposDat->Integral();
				double nSim = hCathHitposSim->Integral();
				cs.fill_frac = nDat/nSim;
				OM.qOut.insert("cathseg",cathseg2sm(cs));
				
				// determine energy re-binning
				std::vector<TH1F*> hEnDats = sliceTH2(*hCathHitposDat,Y_DIRECTION);
				std::vector<TH1F*> hEnSims = sliceTH2(*hCathHitposSim,Y_DIRECTION);
				std::vector<float> enCounts;
				for(unsigned int i=0; i<hEnDats.size(); i++)
					enCounts.push_back(hEnDats[i]->Integral());
				std::vector<unsigned int> part = equipartition(enCounts, 80);
				
				// fit each energy range
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
				OM.printCanvas(sideSubst("Cathodes/Positions_%c",s)+(d==X_DIRECTION?"x":"y")+itos(c));
			}
		}
	}
	
	OM.write();
	OM.setWriteRoot(true);
}

void processWirechamberCal(RunNum r0, RunNum r1, unsigned int nrings) {
	std::string basePath = getEnvSafe("UCNA_ANA_PLOTS")+"/PositionMaps/";
	OutputManager OM("NameUnused",basePath);
	std::string readname = itos(r0)+"-"+itos(r1)+"_"+itos(nrings);
	MWPCTuningAnalyzer WCdat(&OM, "NameUnused", basePath+"/Xenon_"+readname+"/Xenon_"+readname);
	MWPCTuningAnalyzer WCsim(&OM, "NameUnused", basePath+"/SimXe_"+readname+"/SimXe_"+readname);
	processWirechamberCal(WCdat,WCsim);
}
