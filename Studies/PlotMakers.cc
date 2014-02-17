#include "PlotMakers.hh"
#include <TH2F.h>
#include <TProfile.h>
#include <TGraphErrors.h>
#include <TStyle.h>
#include <TRandom.h>
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
#include "LinHistCombo.hh"
#include "PostOfficialAnalyzer.hh"
#include "BetaDecayAnalyzer.hh"
#include "AnodeGainMapPlugin.hh"
#include "AsymmetryCorrections.hh"
#include <TColor.h>


TRandom3 plRndSrc;

void dumpCalInfo(const std::vector<RunNum>& runs, QFile& Qout) {
	for(unsigned int i=0; i<runs.size(); i++) {
		PMTCalibrator LC(runs[i]);
		Stringmap m = LC.calSummary();
		printf("----- %i -----\n",runs[i]);
		m.display();
		Qout.insert("runcal",m);
	}
	Qout.commit();
}

void dumpPosmap(QFile& qOut, PositioningCorrector& PCor, PMTCalibrator* PCal) {

	for(Side s = EAST; s<=WEST; ++s) {
		
		unsigned int nMaps = PCor.getNMaps(s);
		for(unsigned int t=0; t<nMaps; t++) {
			
			const SectorCutter& psects = PCor.getSectors(s,t);
			
			Stringmap sm = SCtoSM(psects);
			sm.insert("side",sideSubst("%c",s));
			sm.insert("tube",itos(t));
			qOut.insert("SectorCutter",sm);
			
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
					m.insert("light",PCor.eval(s,t,x,y,true));
				if(PCal && nMaps>=nBetaTubes)
					m.insert("nPE",PCal->nPE(s, t, 1000, x, y, 0));
				m.insert("energy",1.0);
				
				qOut.insert("PosmapPoint",m);
			}
		}
	}
}


void dumpPosmap(std::string basepath, unsigned int pnum) {
	
	PositioningCorrector* PCor;
	PMTCalibrator* PCal = NULL;
	
	if(pnum < 5000) {
		PCor = CalDBSQL::getCDB()->getPositioningCorrectorByID(pnum);
	} else {
		PCal = new PMTCalibrator(pnum);
		PCor = CalDBSQL::getCDB()->getPositioningCorrector(pnum);
	}
	
	if(!PCor) {
		printf("***ERROR*** Position map %i not found.\n",pnum);
		return;
	}
	
	makePath(basepath);
	QFile fout(basepath+"/Posmap_"+itos(pnum)+".txt",false);
	dumpPosmap(fout, *PCor, PCal);
	fout.commit();
}


void SimSpectrumInfo(Sim2PMT& S, OutputManager& OM) {
	double emax = 1000;
	double nsegs = 100;
	double nbins = nsegs*5;
	
	// set up histograms
	std::vector<TH1*> hEnergy[2][TYPE_III_EVENT+1][2];	// energy [side][type][vis/true] for each input bin
	TProfile* inputEnergy[2][TYPE_III_EVENT+1][2];		// input bins [side][type][vis/true]
	TProfile* inputTotal = (TProfile*)OM.addObject(new TProfile("inputTotal","Input Visible Energy Total",nsegs,0,emax));
	for(Side s = EAST; s <= WEST; ++s) {
		for(unsigned int t = TYPE_0_EVENT; t <= TYPE_III_EVENT; t++) {
			inputEnergy[s][t][false] = (TProfile*)OM.addObject(new TProfile((sideSubst("inputEvis_%c",s)+itos(t)).c_str(),"Input Visible Energy",nsegs,0,emax));
			inputEnergy[s][t][true] = (TProfile*)OM.addObject(new TProfile((sideSubst("inputEtrue_%c",s)+itos(t)).c_str(),"Input True Energy",nsegs,0,emax));
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
		if(S.fPID != PID_BETA || S.fType > TYPE_III_EVENT || S.radius(S.fSide)>50.)
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
			for(unsigned int t = TYPE_0_EVENT; t <= TYPE_III_EVENT; t++) {
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
		for(unsigned int t = TYPE_0_EVENT; t <= TYPE_III_EVENT; t++) {
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

void makeCorrectionsFile(int A, int Z, double Endpt, double M2_F, double M2_GT) {
	
	std::string fout = getEnvSafe("UCNA_ANA_PLOTS")+"/SpectrumCorrection/SpectrumCorrection_";
	fout += itos(A)+"_"+itos(Z)+"_"+dtos(Endpt)+".txt";
	
	BetaSpectrumGenerator BSG(A,Z,Endpt);
	BSG.M2_F = M2_F;
	BSG.M2_GT = M2_GT;
	if(A==1 && Z==1) { BSG.M2_F = 1; BSG.M2_GT = 3; } // neutron case
		
	QFile Q;
	Stringmap sm;
	sm.insert("A",A);
	sm.insert("Z",Z);
	sm.insert("endpt",Endpt);
	sm.insert("W0",BSG.W0);
	sm.insert("R",BSG.R);
	sm.insert("M0",BSG.M0);
	sm.insert("M2_F",BSG.M2_F);
	sm.insert("M2_GT",BSG.M2_GT);
	Q.insert("decayInfo",sm);
	
	for(double e = 0.5; e < Endpt; e+=1.) {
		Stringmap m;
		double W = e/m_e+1.;
		m.insert("energy",e);
		m.insert("W",W);
		m.insert("beta",beta(e));
		m.insert("F0m1",WilkinsonF0(Z,W,BSG.R)-1.0);
		m.insert("L0m1",WilkinsonL0(Z,W,BSG.R)-1.0);
		m.insert("RVm1",WilkinsonRV(W,BSG.W0,BSG.M0)-1.0);
		m.insert("RAm1",WilkinsonRA(W,BSG.W0,BSG.M0)-1.0);
		m.insert("Rm1",CombinedR(W,M2_F,M2_GT,BSG.W0,BSG.M0)-1.0);
		m.insert("BiRWM",Bilenkii59_RWM(W));
		m.insert("VCm1",WilkinsonVC(Z,W,BSG.W0,BSG.R)-1.0);
		m.insert("ACm1",WilkinsonAC(Z,W,BSG.W0,BSG.R)-1.0);
		m.insert("Cm1",CombinedC(Z,W,M2_F,M2_GT,BSG.W0,BSG.R)-1.0);
		m.insert("Qm1",WilkinsonQ(Z,W,BSG.W0,BSG.M0)-1.0);
		m.insert("g",Wilkinson_g_a2pi(W,BSG.W0,BSG.M0));
		m.insert("gS",Sirlin_g_a2pi(e,Endpt));
		m.insert("hmg",shann_h_minus_g_a2pi(W,BSG.W0));
		m.insert("h",shann_h_a2pi(e,Endpt));
		m.insert("RWM",WilkinsonACorrection(W));
		m.insert("S0",plainPhaseSpace(W,BSG.W0));
		m.insert("S",BSG.decayProb(e));
		m.insert("dSm1",BSG.spectrumCorrectionFactor(W)-1.0);
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
	interpogrid->GetYaxis()->SetTitleOffset(1.2);
	return interpogrid;
}

void PosPlotter::npePlot(PMTCalibrator* PCal) {
	
	const SectorCutter& Sects = PCal->P->getSectors(EAST,0);
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
			OM->printCanvas(sideSubst("nPE_%c",s)+itos(t));
		}
	}
}

void PosPlotter::npeGradPlot(PMTCalibrator* PCal) {
	
	const SectorCutter& Sects = PCal->P->getSectors(EAST,0);
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
			printf("**** %c%i RMS gradient/light = %g 1/mm ****\n",sideNames(s),t,gradsum);
			interpogrid->Draw("COL Z");
			drawSectors(Sects,6);
			OM->printCanvas(sideSubst("Posmap_Grad_%c",s)+itos(t));
		}
	}
}

void PosPlotter::diffPlot(const PositioningCorrector& P1, const PositioningCorrector& P2, double zRange) {
	const SectorCutter& Sects = P1.getSectors(EAST,0);
	r0 = Sects.r;
	
	for(Side s = EAST; s<=WEST; ++s) {
		assert(P1.getNMaps(s)==P2.getNMaps(s));
		for(unsigned int t=0; t<P1.getNMaps(s); t++) {
			std::string hTitle = sideSubst("Position Maps %c",s)+itos(t+1)+" Difference";
			TH2F* interpogrid = makeHisto(sideSubst("%c ",s)+itos(t+1), hTitle);
			interpogrid->SetAxisRange(-zRange,zRange,"Z");
			startScan(interpogrid);
			unsigned int npts = 0;
			double sx = 0;
			double sxx = 0;
			while(nextPoint()) {
				double z1 = P1.eval(s, t, x, y, true);
				double dx = 100.0*(z1-P2.eval(s, t, x, y, true))/z1;
				sx += dx;
				sxx += dx*dx;
				npts++;
				dx = dx>zRange ? zRange:(dx<-zRange ? -zRange:dx);
				interpogrid->SetBinContent(nx,ny,dx);
				
			}
			sx /= npts;
			sxx /= npts;
			printf("%c%i RMS %% difference: %.3f\n",sideNames(s),t,sqrt(sxx-sx*sx));
			interpogrid->Draw("COL Z");
			drawSectors(Sects,6);
			OM->printCanvas(sideSubst("PosmapDiff_%c",s)+itos(t));
		}
	}
}

void PosPlotter::npeDiffPlot(const PMTCalibrator& P1, const PMTCalibrator& P2) {
	const SectorCutter& Sects = P1.P->getSectors(EAST,0);
	r0 = Sects.r;
	for(Side s = EAST; s<=WEST; ++s) {
		TH2F* interpogrid = makeHisto(sideSubst("npe_diff_%c",s),"nPE difference");
		interpogrid->SetAxisRange(-4,4,"Z");
		unsigned int npts = 0;
		double sx = 0;
		double sxx = 0;
		startScan(interpogrid);
		while(nextPoint()) {
			float n1 = P1.nPE(s,nBetaTubes,1000,x,y);
			float n2 = P2.nPE(s,nBetaTubes,1000,x,y);
			float pdif = 100.*(n1-n2)/n1;
			sx += pdif;
			sxx += pdif*pdif;
			npts++;
			interpogrid->SetBinContent(nx,ny,n1-n2);
		}
		sx /= npts;
		sxx /= npts;
		printf("%c RMS nPE %% difference: %.3f\n",sideNames(s),sqrt(sxx-sx*sx));
		interpogrid->Draw("COL Z");
		drawSectors(Sects,6);
		OM->printCanvas(sideSubst("PosmapEtaDiff_%c",s));
	}
}

void PosPlotter::etaPlot(PositioningCorrector* P, double z0, double z1) {
	
	const SectorCutter& Sects = P->getSectors(EAST,0);
	r0 = Sects.r;
	
	for(Side s = EAST; s<=WEST; ++s) {
		for(unsigned int t=0; t<P->getNMaps(s); t++) {
			std::string hTitle = sideSubst("%s PMT ",s)+itos(t+1)+" Light Transport";
			if(P->getNMaps(s) == 1) {
				P->setNormAvg();
				hTitle = sideSubst("%s MWPC Charge Scaling", s);
			}
			TH2F* interpogrid = makeHisto(sideSubst("%c ",s)+itos(t+1), hTitle);
			interpogrid->SetAxisRange(z0,z1,"Z");
			
			startScan(interpogrid);
			while(nextPoint())
				interpogrid->SetBinContent(nx,ny,P->eval(s, t, x, y, true));
			
			interpogrid->Draw("COL Z");
			drawSectors(Sects,6);
			OM->printCanvas(sideSubst("Posmap_%c_",s)+itos(t));
		}
	}
}

//-------------------------------------------------------------//

void showSimSpectrum(const std::string& nm, OutputManager& OM, NucDecayLibrary& NDL, PMTCalibrator& PCal) {
	double emax = 1600;
	int nbins = 1600;
	NucDecaySystem& NDS = NDL.getGenerator(nm);
	NDS.display(true);
	PMTGenerator PGen;
	PGen.setCalibrator(&PCal);
	PGen.setSide(EAST);
	
	TH1F* hSpec = OM.registeredTH1F("hSpec","",nbins,0,emax);
	TH1F* hBlurd = OM.registeredTH1F("hBlurd","",nbins/5,0,emax);

	std::vector<NucDecayEvent> v;
	unsigned int npts = 0;
	unsigned int nToGen = 5e5;
	while(npts<nToGen) {
		if(!(npts%(nToGen/20))) { printf("*"); fflush(stdout); npts++; }
		v.clear();
		NDS.genDecayChain(v);
		double eQ = 0;
		for(std::vector<NucDecayEvent>::iterator it = v.begin(); it < v.end(); it++) {
			if(it->d == D_ELECTRON) {
				hSpec->Fill(it->E);
				if(plRndSrc.Uniform()<0.5)
					eQ += it->E;
			}
		}
		if(eQ) {
			ScintEvent sevt = PGen.generate(eQ);
			if(PGen.triggered()) {
				hBlurd->Fill(sevt.energy.x);
				npts++;
			}
		}
	}
	printf(" Done.\n");
	OM.defaultCanvas->SetLogy(true);
	hSpec->Draw();
	OM.printCanvas(nm+"_GenSpectrum");
	OM.defaultCanvas->SetLogy(false);
	hBlurd->Draw();
	OM.printCanvas(nm+"_DetResponse");
}

void compareXenonSpectra() {
	
	std::string isotlist[] = {
		"Xe125_1-2+",	"Xe127_1-2+",	"Xe129_11-2-",
		"Xe131_11-2-",	"Xe133_3-2+",	"Xe133_11-2-",
		"Xe135_3-2+",	"Xe135_11-2-",	"Xe137_7-2-"};
	std::vector<std::string> isots(isotlist,isotlist+9);
	
	OutputManager OM("XeIsots",getEnvSafe("UCNA_ANA_PLOTS")+"/test/XeIsots");
	NucDecayLibrary NDL(getEnvSafe("UCNA_AUX")+"/NuclearDecays",1e-6);
	NDL.BEL.display();
	
	PMTCalibrator PCal(16000);
	gStyle->SetOptStat("");
	
	for(unsigned int n=0; n<isots.size(); n++) {
		printf("\n\n---------------------- %s ---------------------------\n",isots[n].c_str());
		showSimSpectrum(isots[n],OM,NDL,PCal);
	}
	
	OM.setWriteRoot(true);
	OM.write();
}

void decomposeXenon(RunNum rn, bool includeFast) {
	
	gStyle->SetOptStat("");
	
	OutputManager OM("XeDecomp_"+itos(rn),getEnvSafe("UCNA_ANA_PLOTS")+"/test/XeDecomp");
	NucDecayLibrary NDL(getEnvSafe("UCNA_AUX")+"/NuclearDecays",1e-6);
	
	std::vector<std::string> isots;
	isots.push_back("Xe125_1-2+");
	//isots.push_back("Xe127_1-2+");	// unlikely 36d HL, rare parent
	isots.push_back("Xe129_11-2-");
	isots.push_back("Xe131_11-2-");
	isots.push_back("Xe133_3-2+");
	isots.push_back("Xe133_11-2-");
	isots.push_back("Xe135_3-2+");
	if(includeFast) {
		isots.push_back("Xe135_11-2-");
		isots.push_back("Xe137_7-2-");
	}
	
	double emax = 1200;
	int nbins = 240;
	double fidrad = 50.;
	
	// fill true energy spectrum
	TH1F* hSpec = OM.registeredTH1F("hSpec","Data Spectrum",nbins,0,emax);
	hSpec->Sumw2();
	hSpec->SetLineColor(2);
	PostOfficialAnalyzer POA(true);
	POA.addRun(rn);
	POA.startScan();
	while(POA.nextPoint()) {
		Side s = POA.fSide;
		if(POA.fPID == PID_BETA && POA.fType <= TYPE_III_EVENT && (s==EAST || s==WEST) && POA.radius(s)<fidrad)
			hSpec->Fill(POA.getErecon());
	}
	
	SectorCutter SC(5,fidrad);
	
	// simulate each isotope
	std::vector<TH1F*> hIsot;
	LinHistCombo LHC;
	for(unsigned int i=0; i<isots.size(); i++) {
		std::string g4dat = getEnvSafe("UCNA_SIM_PLOTS") + "/20120917_";
		G4SegmentMultiplier g2p(SC);
		g2p.addFile(g4dat + isots[i] + "/analyzed_*.root");
		assert(g2p.getnFiles());
		g2p.setCalibrator(*POA.ActiveCal);
		TH1F* hSim = OM.registeredTH1F(isots[i]+"_Sim","",nbins,0,emax);
		hSim->GetXaxis()->SetTitle("Energy [keV]");
		g2p.startScan(hSpec->GetEntries());
		double simFactor = (isots[i]=="Xe135_3-2+")?0.75:0.25;
		while(hSim->GetEntries() < simFactor*hSpec->GetEntries()) {
			g2p.nextPoint();
			Side s = g2p.fSide;
			if(g2p.fPID == PID_BETA && g2p.fType <= TYPE_III_EVENT && (s==EAST || s==WEST) && g2p.radius(s)<50.)
				hSim->Fill(g2p.getErecon());
		}
		hIsot.push_back((TH1F*)hSim->Clone());
		LHC.addTerm(hSim);
	}
	
	
	// fit composition
	double w = hSpec->GetBinWidth(1);
	hSpec->Scale(w/hSpec->Integral());
	LHC.getFitter()->SetLineColor(4);
	LHC.getFitter()->SetLineWidth(1);
	LHC.getFitter()->SetNpx(nbins);
	LHC.forceNonNegative();
	LHC.Fit(hSpec,75,emax);
	TH1F* hSim = OM.registeredTH1F("hSim","Simulated Spectrum",nbins,0,emax);
	Stringmap m;
	m.insert("run",rn);
	m.insert("counts",hSpec->GetEntries());
	for(unsigned int i=0; i<isots.size(); i++) {
		hIsot[i]->Scale(LHC.coeffs[i]);
		hIsot[i]->SetLineStyle(3);
		hSim->Add(hIsot[i]);
		m.insert("name_"+itos(i),isots[i]);
		m.insert("prop_"+itos(i),hIsot[i]->Integral()/w);
		m.insert("err_"+itos(i),LHC.dcoeffs[i]/LHC.coeffs[i]);
	}
	m.display();
	OM.qOut.insert("xedecomp",m);
	
	// draw
	hSpec->SetMinimum(0);
	hSpec->SetTitle("Xenon spectrum decomposition");
	hSpec->GetXaxis()->SetTitle("Energy [keV]");
	hSpec->Draw();
	for(unsigned int i=0; i<isots.size(); i++)
		hIsot[i]->Draw("Same");
	OM.printCanvas("XeDecomp_"+itos(rn));
	
	OM.write();
	//OM.setWriteRoot(true);
}

void lowStatsTest() {
	int muMax=100;
	unsigned int nBins = 10000;
	TF1 cFit("cFit","pol0",0,nBins);
	
	for(double mu = 0; mu <= muMax; mu++) {
		TH1F hCounts("hCounts","hCounts",nBins,-0.5,nBins+0.5);
		hCounts.Sumw2();
		for(unsigned int i=0; i<nBins; i++) {
			int n = plRndSrc.Poisson(mu);
			//hCounts.SetBinContent(i+1,n+1.-(n>0?1./n:0));
			//hCounts.SetBinError(i+1,sqrt(n+1));
			hCounts.SetBinContent(i+1,n);
			hCounts.SetBinError(i+1,sqrt(n));
		}
		hCounts.Fit(&cFit,"Q");
		double hint = hCounts.Integral()/nBins;
		printf("mu = %g:\tobs = %.1f\t%.1f\t%.1f\t%.1f\n",mu,cFit.GetParameter(0),cFit.GetParameter(0)-mu,hint,hint-mu);
	}
}

//-------------------------------------------------------------//

void separate23(std::string simName) {
	OutputManager OM("23Separation",getEnvSafe("UCNA_ANA_PLOTS")+"/test/23Separation_"+simName+"/");
	RunAccumulator RA(&OM, "RunAccumulator", getEnvSafe("UCNA_ANA_PLOTS")+"/OctetAsym_Offic_"+simName+"/OctetAsym_Offic_"+simName);
	WirechamberSimBackscattersPlugin* WS = new WirechamberSimBackscattersPlugin(&RA);
	RA.addPlugin(WS);
	WS->make23SepInfo(OM);
	OM.write();
}

//-------------------------------------------------------------//

void refitXeAnode(std::string datname) {
	OutputManager OM("NameUnused",getEnvSafe("UCNA_ANA_PLOTS")+"/test/");
	RunAccumulator RA(&OM, "RunAccumulator", datname);
	AnodeGainMapPlugin* AP = new AnodeGainMapPlugin(&RA,0);
	RA.addPlugin(AP);
	AP->genAnodePosmap();
}

//-------------------------------------------------------------//

void calcAnalysisChoices(OutputManager& OM, const std::string& inflname) {
	for(AnalysisChoice ac = ANCHOICE_A; ac <= ANCHOICE_K; ++ac) {
		OctetAnalyzer OA(&OM, "Anchoice_"+ctos(choiceLetter(ac)), inflname);
		AsymmetryPlugin* AA = new AsymmetryPlugin(&OA);
		AA->anChoice = ac;
		OA.addPlugin(AA);
		OA.calculateResults();
		OA.makePlots();
		OA.write();
	}
}

//-------------------------------------------------------------//

void paperDataPlot() {
	
	OutputManager OM("PlotData",getEnvSafe("UCNA_ANA_PLOTS")+"/Paper/");
	
	std::vector<std::string> rPaths;
	rPaths.push_back(getEnvSafe("UCNA_ANA_PLOTS")+"/OctetAsym_Offic/Range_0-16");
	rPaths.push_back(getEnvSafe("UCNA_ANA_PLOTS")+"/OctetAsym_Offic/Range_17-1000");
	
	std::vector<double> rWeights;
	rWeights.push_back(0.59);
	rWeights.push_back(0.41);
	
	std::vector<AsymmetryPlugin*> rAsyms;
	
	double tWeight = 0;
	for(unsigned int n=0; n<rPaths.size(); n++) {
		OctetAnalyzer* OAdat = new OctetAnalyzer(&OM, "DataCorrector", rPaths[n]+"/OctetAsym_Offic");
		AsymmetryPlugin* AAdat = new AsymmetryPlugin(OAdat);
		OAdat->addPlugin(AAdat);
		AAdat->anChoice = ANCHOICE_C;
		doFullCorrections(*AAdat,OM);
		
		// weighted asymmetry sum
		AAdat->hAsym->Scale(rWeights[n]);
		for(GVState gv = GV_CLOSED; gv <= GV_OPEN; ++gv)
			AAdat->hSuperSum[GV_OPEN]->SetBit(TH1::kIsAverage);
		AAdat->hCxn->SetBit(TH1::kIsAverage);
		if(n>0) {
			AAdat->hAsym->Add(rAsyms.back()->hAsym);
			AAdat->hCxn->Add(AAdat->hCxn,rAsyms.back()->hCxn,rWeights[n],tWeight);
			for(GVState gv = GV_CLOSED; gv <= GV_OPEN; ++gv)
				AAdat->hSuperSum[gv]->Add(AAdat->hSuperSum[gv],rAsyms.back()->hSuperSum[gv],rWeights[n],tWeight);
		}
		tWeight += rWeights[n];
		rAsyms.push_back(AAdat);
	}
	
	AsymmetryPlugin* AA = rAsyms.back();
	
	// simulation data for comparison
	OctetAnalyzer* OAsim = new OctetAnalyzer(&OM, "DataCorrector", getEnvSafe("UCNA_ANA_PLOTS")+"/OctetAsym_Offic_Sim0823_4x/OctetAsym_Offic_Sim0823_4x");
	AsymmetryPlugin* AAsim = new AsymmetryPlugin(OAsim);
	OAsim->addPlugin(AAsim);
	AAsim->anChoice = ANCHOICE_C;
	doFullCorrections(*AAsim,OM);
	AAsim->hSuperSum[GV_OPEN]->Scale(AA->hSuperSum[GV_OPEN]->Integral()/AAsim->hSuperSum[GV_OPEN]->Integral());
	
	// save data points to file
	for(int b = 1; b <= AA->hAsym->GetNbinsX(); b++) {
		Stringmap m;
		double e0 = AA->hAsym->GetBinCenter(b);
		if(e0>800) continue;
		m.insert("KE",e0);
		m.insert("A0",AA->hAsym->GetBinContent(b));
		m.insert("dA0",AA->hAsym->GetBinError(b));
		m.insert("corr",AA->hCxn->GetBinContent(b));
		m.insert("dcorr",AA->hCxn->GetBinError(b));
		m.insert("ssfg",AA->hSuperSum[GV_OPEN]->GetBinContent(b));
		m.insert("ssMC",AAsim->hSuperSum[GV_OPEN]->GetBinContent(b));
		m.insert("ssbg",AA->hSuperSum[GV_CLOSED]->GetBinContent(b));
		OM.qOut.insert("asymplot",m);
	}
	
	OM.write();
}


