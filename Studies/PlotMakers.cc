#include "PlotMakers.hh"
#include <TH2F.h>
#include <TProfile.h>
#include <TGraphErrors.h>
#include <TStyle.h>
#include <map>
#include <cassert>
#include <vector>
#include "GraphicsUtils.hh"
#include "QFile.hh"
#include "PathUtils.hh"
#include "ReSource.hh"
#include "Types.hh"
#include "DataSource.hh"
#include "KurieFitter.hh"
#include "G4toPMT.hh"
#include "BetaSpectrum.hh"

void processCalibrationsFile(const std::string& fInName, const std::string& fOutName) {
	
	QFile qIn(fInName,true);
	QFile qOut(fOutName,false);
	std::map<RunNum,PMTCalibrator*> PCals;
	
	std::vector<Stringmap> srs = qIn.retrieve("SourceLine");
	for(std::vector<Stringmap>::iterator sit = srs.begin(); sit != srs.end(); sit++) {
		
		RunNum rn = RunNum(sit->getDefault("runNum",0));
		std::string sname = sit->getDefault("side","N");
		Side s = NONE;
		if(sname == "E")
			s = EAST;
		else if(sname == "W")
			s = WEST;
		unsigned int t = (unsigned int)(sit->getDefault("tube",1000));
		float x = sit->getDefault("x",0);
		float y = sit->getDefault("y",0);
		float adc = sit->getDefault("adc",0);
		float dadc = sit->getDefault("dadc",0);
		float evis = sit->getDefault("evis",0);
		float evisExp = sit->getDefault("evisExp",0);
		
		std::map<RunNum,PMTCalibrator*>::iterator it = PCals.find(rn);
		if(it == PCals.end()) {
			PCals.insert(std::make_pair(rn,new PMTCalibrator(rn, CalDBSQL::getCDB())));
			it = PCals.find(rn);
			it->second->printSummary();
		}
		
		if(t==nBetaTubes) {
			sit->insert("eta",it->second->combEta(s,x,y));
		}
		if(t<nBetaTubes) {
			sit->insert("eta",it->second->eta(s,t,x,y));
			sit->insert("lvis",it->second->linearityCorrector(s, t, adc, 0));
			sit->insert("dlinearity",it->second->dLinearity(s, t, adc, 0)*dadc);
		}
		
		qOut.insert("SourceLine",*sit);
	}
	
	qOut.commit();
}


void plotGMScorrections(const std::vector<RunNum>& runs, const std::string& foutPath) {
	QFile fout(foutPath+"/GMS_Plot.txt",false);
	for(unsigned int i=0; i<runs.size(); i++) {
		PMTCalibrator LC(runs[i],CalDBSQL::getCDB());
		Stringmap m = LC.calSummary();
		printf("----- %i -----\n",runs[i]);
		m.display();
		fout.insert("gmsinfo",m);
	}
	fout.commit();
}

void etaPlot(OutputManager& OM, PositioningCorrector* P, bool normalize, float axisRange) {
	
	for(Side s = EAST; s<=WEST; ++s) {
		for(int t=0; t<=nBetaTubes; t++) {
			if(!P->eval(s, t, 0, 0))
				continue;
			TH2F interpogrid("Interpogrid",(sideSubst("%c ",s)+itos(t+1)+" Interpolated").c_str(),200,-60,60,200,-60,60);
			interpogrid.SetAxisRange(0,axisRange,"Z");
			float x,y;
			float r0 = P->getSectors(s,t).r;
			float rscale = 1.2;
			for(int nx=1; nx<=interpogrid.GetNbinsX(); nx++) {
				for(int ny=1; ny<=interpogrid.GetNbinsY(); ny++) {
					x = interpogrid.GetXaxis()->GetBinCenter(nx);
					y = interpogrid.GetYaxis()->GetBinCenter(ny);
					if(x*x+y*y <= r0*r0*rscale)
						interpogrid.SetBinContent(nx,ny,P->eval(s, t, x, y, normalize));
					else
						interpogrid.SetBinContent(nx,ny,0.0);
				}
			}
			interpogrid.Draw("COL Z");
			drawSectors(P->getSectors(s,t),6);
			OM.printCanvas(sideSubst("Posmap_%c_",s)+itos(t));
		}
	}
}

void etaGradPlot(OutputManager& OM, PositioningCorrector* P) {
	
	for(Side s = EAST; s<=WEST; ++s) {
		for(int t=0; t<nBetaTubes; t++) {
			TH2F interpogrid("Interpogrid",(sideSubst("%c ",s)+itos(t+1)+" Interpolated").c_str(),200,-60,60,200,-60,60);
			interpogrid.SetAxisRange(0,10,"Z");
			float x,y;
			float r0 = P->getSectors(s,t).r;
			float gradsum = 0;
			float nn = 0;
			for(int nx=1; nx<interpogrid.GetNbinsX(); nx++) {
				for(int ny=1; ny<interpogrid.GetNbinsY(); ny++) {
					x = interpogrid.GetXaxis()->GetBinCenter(nx);
					y = interpogrid.GetYaxis()->GetBinCenter(ny);
					if(x*x+y*y <= r0*r0*1.001) {
						float e0 = P->eval(s, t, x, y, true);
						float e1 = P->eval(s, t, x+0.01, y, true);
						float e2 = P->eval(s, t, x, y+0.01, true);
						float gd = sqrt( (e1-e0)*(e1-e0)+(e2-e0)*(e2-e0) )/0.01/e0;
						gradsum += gd*gd;
						nn += 1.0;
						interpogrid.SetBinContent(nx,ny,100.0*gd);
					} else {
						interpogrid.SetBinContent(nx,ny,0.0);
					}
				}
			}
			gradsum = sqrt(gradsum/nn);
			printf("**** %c%i RMS gradient = %g 1/mm ****\n",sideNames(s),t,gradsum);
			interpogrid.Draw("COL Z");
			drawSectors(P->getSectors(s,t),6);
			OM.printCanvas(sideSubst("Posmap_Grad_%c",s)+itos(t));
		}
	}
}


void dumpPosmap(std::string basepath, unsigned int pnum) {
	
	PositioningCorrector* PCor;
	PMTCalibrator* PCal = NULL;
	
	if(pnum < 5000) {
		PCor = CalDBSQL::getCDB()->getPositioningCorrectorByID(pnum);
	} else {
		PCal = new PMTCalibrator(pnum,CalDBSQL::getCDB());
		PCor = CalDBSQL::getCDB()->getPositioningCorrector(pnum);
	}
	
	if(!PCor) {
		printf("***ERROR*** Position map %i not found.\n",pnum);
		return;
	}
	
	makePath(basepath);
	QFile fout(basepath+"/Posmap_"+itos(pnum)+".txt",false);
	
	for(Side s = EAST; s<=WEST; ++s) {
		for(int t=0; t<=nBetaTubes; t++) {
			
			if(t==nBetaTubes && !PCal)
				continue;
			
			SectorCutter& psects = PCor->getSectors(s,0);
			
			for(unsigned int n=0; n<psects.nSectors(); n++) {
				
				float x,y;
				psects.sectorCenter(n,x,y);
				
				Stringmap m;
				m.insert("side",ctos(sideNames(s)));
				m.insert("tube",itos(t));
				m.insert("sector",n);
				m.insert("x",x);
				m.insert("y",y);
				if(t<nBetaTubes)
					m.insert("light",PCor->eval(s,t,x,y,true)); 
				else
					m.insert("light",PCal->nPE(s, nBetaTubes, 1000, x, y, 0));
				m.insert("energy",1.0);
				
				fout.insert("PosmapPoint",m);
			}
		}
		
	}
	fout.commit();
}

void npePlot(OutputManager& OM, PMTCalibrator* PCal, float e0, float s0, bool dumbsum) {
	
	for(Side s = EAST; s<=WEST; ++s) {
		
		TH2F interpogrid[nBetaTubes+1];
		for(unsigned int t=0; t<=nBetaTubes; t++) {
			interpogrid[t] = TH2F((sideSubst("%c",s)+itos(t+1)+"_nPE").c_str(),"nPE",200,-60,60,200,-60,60);
			interpogrid[t].SetAxisRange(0,300,"Z");
			if(t==nBetaTubes)
				interpogrid[t].SetAxisRange(0,600,"Z");
		}
		
		float x,y,npe;
		float r0 = PCal->P->getSectors(s,0).r;
		float rscale = 1.2;
		float npesum = 0;
		float gradsum = 0;
		float nn = 0;
		
		for(int nx=1; nx<=interpogrid[0].GetNbinsX(); nx++) {
			for(int ny=1; ny<=interpogrid[0].GetNbinsY(); ny++) {
				x = interpogrid[0].GetXaxis()->GetBinCenter(nx);
				y = interpogrid[0].GetYaxis()->GetBinCenter(ny);
				for(unsigned int t=0; t<=nBetaTubes; t++) {
					if(x*x+y*y <= r0*r0*rscale) {
						npe = PCal->nPE(s,t,e0*s0,x,y,0)/s0;
						if(t==nBetaTubes && dumbsum)
							npe = PCal->pmtSumPE(s,e0*s0,x,y,0)/s0;
						if(t==nBetaTubes && x*x+y*y<50*50) {
							float npe1 = PCal->nPE(s,t,e0*s0,x+0.01,y,0)/s0;
							float npe2 = PCal->nPE(s,t,e0*s0,x,y+.01,0)/s0;
							if(dumbsum) {
								npe1 = PCal->pmtSumPE(s,e0*s0,x+0.01,y,0)/s0;
								npe2 = PCal->pmtSumPE(s,e0*s0,x,y+.01,0)/s0;
							}
							float gd = sqrt( (1/npe1-1/npe)*(1/npe1-1/npe)+(1/npe2-1/npe)*(1/npe2-1/npe) )/0.01*npe;
							gradsum += gd*gd;
							npesum += npe;
							nn += 1.0;
						}
					} else {
						npe = 0;
					}
					interpogrid[t].SetBinContent(nx,ny,npe);
				}
			}
		}
		
		gradsum = sqrt(gradsum/nn/2);
		printf("**** %c RMS gradient sqrt < |Grad K|^2 / K^2 / 2 > = %g 1/mm ****\n",sideNames(s),gradsum);
		printf("**** %c average nPE = %g ****\n",sideNames(s),npesum/nn);
		
		
		for(unsigned int t=0; t<=nBetaTubes; t++) {
			OM.defaultCanvas->SetRightMargin(0.125);
			OM.defaultCanvas->SetBottomMargin(0.125);
			interpogrid[t].Draw("COL Z");
			drawSectors(PCal->P->getSectors(s,0),6);
			if(dumbsum && t==nBetaTubes)
				OM.printCanvas(sideSubst("nPE_%c_dumbsum",s));
			else
				OM.printCanvas(sideSubst("nPE_%c",s)+itos(t));			
		}
	}
}

// turn two profiles against the same axis into an x-y plot
TGraphErrors* correlateProfiles(TProfile* x, TProfile* y) {
	unsigned int nbins = x->GetNbinsX();
	assert(nbins == (unsigned int)y->GetNbinsX());
	TGraphErrors* g = new TGraphErrors(nbins-2);
	unsigned int n = 0;
	for(unsigned int i=1; i<nbins-1; i++) {
		if(x->GetBinContent(i) > 25) {
			g->SetPoint(n,x->GetBinContent(i),y->GetBinContent(i));
			g->SetPointError(n,x->GetBinError(i),y->GetBinError(i));
			n++;
		} else {
			g->RemovePoint(n);
		}
		
	}
	return g;
}

void SimSpectrumInfo(Sim2PMT& S, OutputManager& OM) {
	double emax = 1000;
	double nsegs = 100;
	double nbins = nsegs*5;
	
	// set up histograms
	std::vector<TH1*> hEnergy[2][TYPE_III_EVENT][2];	// energy [side][type][vis/true] for each input bin
	TProfile* inputEnergy[2][TYPE_III_EVENT][2];		// input bins [side][type][vis/true]
	TProfile* inputTotal = (TProfile*)OM.addObject(new TProfile("inputTotal","Input Visible Energy Total",nsegs,0,emax));
	for(Side s = EAST; s <= WEST; ++s) {
		for(unsigned int t = TYPE_0_EVENT; t <= TYPE_II_EVENT; t++) {
			inputEnergy[s][t][false] = (TProfile*)OM.addObject(new TProfile((std::string("inputEvis_%c",s)+itos(t)).c_str(),"Input Visible Energy",nsegs,0,emax));
			inputEnergy[s][t][true] = (TProfile*)OM.addObject(new TProfile((std::string("inputEtrue_%c",s)+itos(t)).c_str(),"Input True Energy",nsegs,0,emax));
			for(unsigned int i=0; i<nsegs; i++) {
				hEnergy[s][t][false].push_back(OM.registeredTH1F(sideSubst("hEvis_%c",s)+itos(t)+"_"+itos(i),"Visible Energy",nbins,0,emax));
				hEnergy[s][t][true].push_back(OM.registeredTH1F(sideSubst("hEtrue_%c",s)+itos(t)+"_"+itos(i),"True Energy",nbins,0,emax));
				hEnergy[s][t][false].back()->GetXaxis()->SetTitle("Energy [keV]");
				hEnergy[s][t][true].back()->GetXaxis()->SetTitle("Energy [keV]");
				hEnergy[s][t][true].back()->GetYaxis()->SetTicks("");
				hEnergy[s][t][true].back()->GetYaxis()->LabelsOption("");
				hEnergy[s][t][true].back()->GetYaxis()->SetTitle("Counts");
				hEnergy[s][t][false].back()->GetYaxis()->SetTicks("");
				hEnergy[s][t][false].back()->GetYaxis()->LabelsOption("");
				hEnergy[s][t][false].back()->GetYaxis()->SetTitle("Counts");
			}
		}
	}
	
	// scan data
	S.startScan();
	while(S.nextPoint()) {
		inputTotal->Fill(S.ePrim,S.ePrim);
		if(S.fPID != PID_BETA || S.fType > TYPE_II_EVENT || S.radius(S.fSide)>50.)
			continue;
		inputEnergy[S.fSide][S.fType][true]->Fill(S.ePrim,S.ePrim);
		int trueSeg = inputEnergy[S.fSide][S.fType][true]->FindBin(S.ePrim)-1;
		inputEnergy[S.fSide][S.fType][false]->Fill(S.getEnergy(),S.getEnergy());
		int visSeg = inputEnergy[S.fSide][S.fType][false]->FindBin(S.getEnergy())-1;
		if(0 <= trueSeg && trueSeg < nsegs)
			hEnergy[S.fSide][S.fType][false][trueSeg]->Fill(S.getEnergy());
		if(0 <= visSeg && visSeg < nsegs)
			hEnergy[S.fSide][S.fType][true][visSeg]->Fill(S.ePrim);
	}
	
	// collect histogram data
	for(unsigned int i=0; i<nsegs; i++) {
		for(Side s = EAST; s <= WEST; ++s) {
			for(unsigned int t = TYPE_0_EVENT; t <= TYPE_II_EVENT; t++) {
				Stringmap m;
				
				m.insert("segment",i);
				m.insert("side",sideSubst("%s",s));
				m.insert("type",t);
				
				m.insert("etrue_input",inputEnergy[s][t][true]->GetBinContent(i+1));
				m.insert("etrue_input_w",inputEnergy[s][t][true]->GetBinError(i+1));
				m.insert("etrue_input_n",inputEnergy[s][t][true]->GetBinEntries(i+1));
				m.insert("etrue_avg",hEnergy[s][t][true][i]->GetMean());
				m.insert("etrue_rms",hEnergy[s][t][true][i]->GetRMS());
				m.insert("etrue_counts",hEnergy[s][t][true][i]->Integral());
				
				m.insert("evis_input",inputEnergy[s][t][false]->GetBinContent(i+1));
				m.insert("evis_input_w",inputEnergy[s][t][false]->GetBinError(i+1));
				m.insert("evis_avg",hEnergy[s][t][false][i]->GetMean());
				m.insert("evis_rms",hEnergy[s][t][false][i]->GetRMS());
				m.insert("evis_counts",hEnergy[s][t][false][i]->Integral());
				
				m.insert("total_counts",inputTotal->GetBinEntries(i+1));
				
				OM.qOut.insert("spectrumInfo",m);
			}
		}
	}
	
	// plots
	OM.defaultCanvas->SetLogy(false);
	OM.defaultCanvas->cd();
	gStyle->SetOptStat("");
	for(Side s = EAST; s <= WEST; ++s) {
		for(unsigned int t = TYPE_0_EVENT; t <= TYPE_II_EVENT; t++) {
			for(unsigned int n=0; n<=1; n++) {
				for(unsigned int i=0; i<nsegs; i++) {
					hEnergy[s][t][n][i]->Scale(1.0/hEnergy[s][t][n][i]->GetBinWidth(1));
					hEnergy[s][t][n][i]->SetLineColor(2+(i%6));
				}
				double hmax = drawSimulHistos(hEnergy[s][t][n]);
				TH1D* hCounts = inputEnergy[s][t][1-n]->ProjectionX("_px","B");
				hCounts->Scale(hmax/hCounts->GetMaximum());
				hCounts->SetLineWidth(3);
				hCounts->Draw("SAME");
				//for(unsigned int i=0; i<nsegs; i++)
				//	drawVLine(hEnergy[s][t][n][i]->GetMean(),OM.defaultCanvas,2+(i%6));
				OM.printCanvas(sideSubst("Energies_%c_Type_",s)+itos(t)+"_"+itos(n));
			}
		}
	}
}


void ProcessLineSims(const std::string& baseDir, unsigned int nLines, float eMin, float eMax) {
	/*
	 OutputManager OM("EQuench2ETrue","../PostPlots/EQuench2ETrue/");
	 
	 // make listing of files to consider
	 std::vector<std::string> infnames;
	 std::vector<double> energies;
	 char tmp[1024];
	 for(unsigned int i=0; i<nLines; i++) {
	 energies.push_back(exp(log(eMin)+i*log(eMax/eMin)/(nLines-1.)));
	 sprintf(tmp,"%s_%.1fkeV/analyzed*.root",baseDir.c_str(),energies.back());
	 infnames.push_back(tmp);
	 }
	 
	 // set up PMT generators to use
	 PMTCalibrator PCal(15916,CalDBSQL::getCDB());
	 PMTGenerator* PGens[] = {new PMTGenerator(EAST), new PMTGenerator(WEST)};
	 for(Side s = EAST; s <= WEST; ++s)
	 PGens[s]->setCalibrator(&PCal);
	 
	 unsigned int i=0;
	 for(std::vector<std::string>::iterator it = infnames.begin(); it != infnames.end(); it++) {
	 
	 // load re-simulated data with PMT effects
	 G4toPMT g2p;
	 g2p.setGenerators(PGens[EAST],PGens[WEST]);
	 g2p.addFile(*it);
	 float e0 = energies[i];
	 float sigmaEst = PCal.energyResolution(EAST,nBetaTubes,e0,0,0,0);
	 AsymHists AH(&OM,std::string("SimLine_")+dtos(energies[i]),"",100,e0+5.*sigmaEst);
	 assert(false); // TODO not using eTrue
	 AH.loadSimData(g2p,0,AFP_OTHER);
	 AH.calcAsym(true);
	 AH.makePlots();
	 
	 for(Side s = EAST; s <= WEST; ++s) {
	 for(unsigned int tp = TYPE_0_EVENT; tp <= TYPE_III_EVENT; tp++) {
	 for(unsigned int afp = AFP_OFF; afp <= AFP_ON; afp++) {
	 
	 TH1F* hLine = AH.hTypes[afp][1][s][tp];
	 
	 Stringmap m;
	 m.insert("etrue",energies[i]);
	 m.insert("side",ctos(sideNames(s)));
	 m.insert("afp",itos(afp));
	 m.insert("type",itos(tp));
	 m.insert("nCounts",hLine->Integral());
	 m.insert("evis_mean",hLine->GetMean());
	 m.insert("evis_rms",hLine->GetRMS());
	 
	 // fit peak
	 float peak_est = hLine->GetBinCenter(hLine->GetMaximumBin());
	 sigmaEst = PCal.energyResolution(s,nBetaTubes,peak_est,0,0,0);
	 TF1 gausFit("gasufit","gaus",peak_est-1.5*sigmaEst,peak_est+1.5*sigmaEst);
	 if(!hLine->Fit(&gausFit,"QR")) {
	 m.insert("fitHeight",gausFit.GetParameter(0));
	 m.insert("evis_mp",gausFit.GetParameter(1));
	 m.insert("evis_fitwidth",gausFit.GetParameter(2));
	 m.insert("d_evis_mp",gausFit.GetParError(1));
	 }
	 
	 printf("\n------- %c %i %i %.2f -------\n",sideNames(s),tp,afp,energies[i]);
	 m.display();
	 OM.qOut.insert("linesim",m);
	 }
	 }
	 }
	 i++;
	 }
	 
	 OM.write();
	 OM.setWriteRoot(true);
	 */
}




void XeEndpointStudy() {
	// load data
	OutputManager OM("XeEndpoint","../PostPlots/XeEndpoint");
	std::vector<RunNum> xeRunNums;
	for(RunNum rn = 16070; rn <= 16077; rn++)
		xeRunNums.push_back(rn);
	PostAnalyzer XeRuns;
	XeRuns.addRuns(xeRunNums);
	
	// set up histograms
	TH1F* hXe[2][nBetaTubes+1];
	for(Side s = EAST; s <= WEST; ++s) {
		for(unsigned int t=0; t<=nBetaTubes; t++) {
			hXe[s][t] = OM.registeredTH1F(sideSubst("XeSpectrum_%c",s)+itos(t),"Xe Spectrum",150,0,1500);
			hXe[s][t]->Sumw2();
			hXe[s][t]->SetLineColor(nBetaTubes+1-t);
		}
	}
	
	// scan data and fill histograms
	XeRuns.startScan();
	while (XeRuns.nextPoint()) {
		for(Side s = EAST; s <= WEST; ++s) {
			if(!(XeRuns.fSide==s && XeRuns.fPID==PID_BETA && XeRuns.fType==TYPE_0_EVENT) || XeRuns.radius(s) > 30.0)
				continue;
			hXe[s][nBetaTubes]->Fill(XeRuns.scints[s].energy.x);
			for(unsigned int t=0; t<nBetaTubes; t++)
				hXe[s][t]->Fill(XeRuns.scints[s].tuben[t].x);
		}
	}
	
	// fits
	for(Side s = EAST; s <= WEST; ++s) {
		for(unsigned int t=0; t<=nBetaTubes; t++) {
			TGraphErrors* tgData = NULL;
			float_err ep = kurieIterator(hXe[s][t],850.0,&tgData,915.,350,700);
			if(tgData) {
				tgData->SetMarkerColor(2);
				tgData->Draw("AP");
				OM.printCanvas(sideSubst("Kurie_Plot_%c",s)+itos(t));
				delete(tgData);
			}
			Stringmap M;
			M.insert("side",ctos(sideNames(s)));
			M.insert("tube",t);
			M.insert("ep",ep.x);
			M.insert("dep",ep.err);
			OM.qOut.insert("endpointFit",M);
		}
	}
	
	// draw histograms
	for(Side s = EAST; s <= WEST; ++s) {
		hXe[s][nBetaTubes]->Draw();
		for(unsigned int t=0; t<nBetaTubes; t++) {
			hXe[s][t]->Draw("Same");
		}
		OM.printCanvas(sideSubst("Spectra_%c",s));
	}
	
	OM.write();
	OM.setWriteRoot(true);
}

void makeCorrectionsFile(const std::string& fout) {
	QFile Q;
	double Z = 1.;
	for(double e = 0.5; e < 800; e+=1.) {
		Stringmap m;
		double W = e/m_e+1.;
		m.insert("energy",e);
		m.insert("W",W);
		m.insert("F0m1",WilkinsonF0(Z,W)-1.0);
		m.insert("L0m1",WilkinsonL0(Z,W)-1.0);
		m.insert("RVm1",WilkinsonRV(W)-1.0);
		m.insert("RAm1",WilkinsonRA(W)-1.0);
		m.insert("BiRWM",Bilenkii_1958_11(W));
		m.insert("VCm1",WilkinsonVC(Z,W)-1.0);
		m.insert("ACm1",WilkinsonAC(Z,W)-1.0);
		m.insert("Qm1",WilkinsonQ(Z,W)-1.0);
		m.insert("g",Wilkinson_g(W));
		m.insert("hmg",shann_h_minus_g(W));
		m.insert("RWM",WilkinsonACorrection(W));
		m.insert("S0",plainPhaseSpace(W));
		m.insert("S",correctedBetaSpectrum(e));
		m.insert("dSm1",spectrumCorrectionFactor(e)-1.0);
		m.insert("A0",plainAsymmetry(e,0.5));
		m.insert("A",correctedAsymmetry(e,0.5));
		m.insert("dAm1",asymmetryCorrectionFactor(e)-1);
		Q.insert("spectrumPoint",m);
	}
	Q.commit(fout);
}
