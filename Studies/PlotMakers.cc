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
#include "KurieFitter.hh"
#include "G4toPMT.hh"
#include "BetaSpectrum.hh"
#include <TColor.h>

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
				hEnergy[s][t][true].back()->GetYaxis()->SetTitle("Counts");
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
				OM.printCanvas(sideSubst("Energies_%c_Type_",s)+itos(t)+"_"+itos(n));
			}
		}
	}
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


//-------------------------------------------------------------//

PosPlotter::PosPlotter(OutputManager* O): rscale(1.0), nbin(100), OM(O) {
	gStyle->SetOptStat("");
	OM->defaultCanvas->SetRightMargin(0.125);
	OM->defaultCanvas->SetBottomMargin(0.125);
}

void PosPlotter::startScan(TH2* h) {
	hCurrent = h;
	nx=ny=0;
	x=y=1.e6;
}

bool PosPlotter::nextPoint() {
	do {
		nx = (nx%hCurrent->GetNbinsX())+1;
		ny += (nx==1);
		x = hCurrent->GetXaxis()->GetBinCenter(nx);
		y = hCurrent->GetYaxis()->GetBinCenter(ny);
	} while(x*x+y*y > r0*r0*rscale && ny <= hCurrent->GetNbinsY());
	return ny <= hCurrent->GetNbinsY();
}

TH2F* PosPlotter::makeHisto(const std::string& nm, const std::string& title) {
	TH2F* interpogrid = OM->registeredTH2F(nm,title,nbin,-60,60,nbin,-60,60);
	interpogrid->GetXaxis()->SetTitle("x Position [mm]");
	interpogrid->GetYaxis()->SetTitle("y Position [mm]");
	return interpogrid;
}

void PosPlotter::npePlot(PMTCalibrator* PCal) {
	
	SectorCutter& Sects = PCal->P->getSectors(EAST,0);
	r0 = Sects.r;
	
	float e0 = 1000;
	float s0 = 0.5;
	
	for(Side s = EAST; s<=WEST; ++s) {
		for(unsigned int t=0; t<=nBetaTubes; t++) {
			
			TH2F* interpogrid = makeHisto(sideSubst("%c",s)+itos(t+1)+"_nPE",
										  (t==nBetaTubes?sideSubst("%s",s):sideSubst("%c",s)+itos(t+1))+" PE per MeV");
			interpogrid->SetAxisRange(0,250,"Z");
			if(t==nBetaTubes)
				interpogrid->SetAxisRange(0,500,"Z");
			
			startScan(interpogrid);
			while(nextPoint())
				interpogrid->SetBinContent(nx,ny,PCal->nPE(s,t,e0*s0,x,y,0)/s0);
			
			interpogrid->Draw("COL Z");
			//drawSectors(Sects,6);
			OM->printCanvas(sideSubst("nPE_%c",s)+itos(t),".eps");			
		}
	}
}

void PosPlotter::npeGradPlot(PMTCalibrator* PCal) {
	
	SectorCutter& Sects = PCal->P->getSectors(EAST,0);
	r0 = Sects.r;
	
	for(Side s = EAST; s<=WEST; ++s) {
		for(int t=0; t<=nBetaTubes; t++) {
			TH2F* interpogrid = makeHisto(sideSubst("%c",s)+itos(t+1)+"_Grad",
										  (t==nBetaTubes?sideSubst("%s",s):sideSubst("%c",s)+itos(t+1))+" Light Transport Gradient");
			interpogrid->SetAxisRange(0,10,"Z");
			
			float gradsum = 0;
			float nn = 0;
			startScan(interpogrid);
			while(nextPoint()) {
				float e0 = PCal->nPE(s,t,500,x,y,0);
				float e1 = PCal->nPE(s,t,500,x+0.01,y,0);
				float e2 = PCal->nPE(s,t,500,x,y+0.01,0);
				float gd = sqrt( (e1-e0)*(e1-e0)+(e2-e0)*(e2-e0) )/0.01/e0;
				gradsum += gd*gd;
				nn += 1.0;
				interpogrid->SetBinContent(nx,ny,100.0*gd);
			}
			gradsum = sqrt(gradsum/nn);
			printf("**** %c%i RMS gradient = %g 1/mm ****\n",sideNames(s),t,gradsum);
			interpogrid->Draw("COL Z");
			drawSectors(Sects,6);
			OM->printCanvas(sideSubst("Posmap_Grad_%c",s)+itos(t));
		}
	}
}

void PosPlotter::etaPlot(PositioningCorrector* P, double axisRange) {
	
	SectorCutter& Sects = P->getSectors(EAST,0);
	r0 = Sects.r;
	
	for(Side s = EAST; s<=WEST; ++s) {
		for(unsigned int t=0; t<nBetaTubes; t++) {
			
			if(!P->eval(s, t, 0, 0))
				continue;
			
			TH2F* interpogrid = makeHisto(sideSubst("%c ",s)+itos(t+1), sideSubst("%s PMT ",s)+itos(t+1)+" Light Transport");
			interpogrid->SetAxisRange(0,axisRange,"Z");
			
			startScan(interpogrid);
			while(nextPoint())
				interpogrid->SetBinContent(nx,ny,P->eval(s, t, x, y, true));
			
			interpogrid->Draw("COL Z");
			drawSectors(Sects,6);
			OM->printCanvas(sideSubst("Posmap_%c_",s)+itos(t));
		}
	}
}
