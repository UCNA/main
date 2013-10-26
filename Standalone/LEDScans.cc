#include "LEDScans.hh"
#include "LED2PMT.hh"
#include "SMExcept.hh"
#include "GraphicsUtils.hh"
#include <TProfile.h>
#include <climits>
#include <cmath>
#include <TLatex.h>
#include "RollingWindow.hh"
#include "PMTGenerator.hh"
#include "GraphUtils.hh"


std::string LEDScanScanner::locateRun(RunNum r) {
	printf("Looking for run %i...\n",r);
	std::string fname = getEnvSafe("UCNAOUTPUTDIR")+"/hists/spec_"+itos(r)+".root";
	if(fileExists(fname))
		return fname;
	printf("*** Replay of run %i not found! ***\n",r);
	return "";
}

void LEDScanScanner::setReadpoints() {	
	// clock
	Tch->SetBranchAddress("Time",&runClock);
	Tch->SetBranchAddress("Sis00",&SIS00);
	for(Side s=EAST; s<=WEST; ++s) {
		// beta scintillators
		Tch->SetBranchAddress(sideSubst("Scint%c",s).c_str(),&scints[s]);
		// MWPC totals
		Tch->SetBranchAddress(sideSubst("Anode%c",s).c_str(),&mwpcs[s].anode);
	}
}


//----------------------------------------------------------------

void AveragerUnit::clear() {
	pts.clear();
	RW.clear();
}

void SweepAverager::start() {
	//printf("Starting SweepAverager for %i units...\n",(int)AUs.size());
	for(std::vector<AveragerUnit*>::iterator it = AUs.begin(); it != AUs.end(); it++)
		(*it)->clear();
	ncached = 0;
}

void SweepAverager::addPt() {
	for(std::vector<AveragerUnit*>::iterator it = AUs.begin(); it != AUs.end(); it++) {
		(*it)->RW.addCount(0,(*it)->x);
		(*it)->pts.push_back((*it)->x);
	}
	ncached++;
	if(ncached > npts) {
		doWithAverage();
		for(std::vector<AveragerUnit*>::iterator it = AUs.begin(); it != AUs.end(); it++)
			(*it)->pts.pop_front();
		ncached--;
	}
}

void SweepAverager::stop() {
	//printf("Completed sweep; processing cached points...\n");
	while(ncached--) {
		for(std::vector<AveragerUnit*>::iterator it = AUs.begin(); it != AUs.end(); it++) {
			(*it)->RW.popExcess();
			(*it)->pts.pop_front();
		}
		doWithAverage();
	}
}

//----------------------------------------------------------------

LEDAnalyzer::LEDAnalyzer(std::string nm, std::string bp): OutputManager(nm,bp), SweepAverager(50), pass(0) {
	
	float binw = 1.0;
	emax = 1500-0.5*binw;
	wmax = 3*5*sqrt(emax);
	float emin = -10+0.5*binw; // minimum energy for plot ranges
	unsigned int nbinsx = (emax-emin)/binw; // number of energy bins for total energy
	unsigned int nbinsy = 200; // number of energy bins for energy deviation
	
	// averagers
	Eavg = addUnit();
	clockAvg = addUnit();
	for(Side s = EAST; s <= WEST; ++s)
		for(unsigned int t=0; t<nBetaTubes; t++)
			Tavg[s][t] = addUnit();

	// set up profiles
	pE0 = new TProfile("EnergyCorrector","Energy averaging error",nbinsx,emin,emax);
	addObject(pE0);
	pE0->GetXaxis()->SetTitle("Energy equivalent [keV]");
	pE0->GetYaxis()->SetTitle("deviation from average [keV]");
	hPed = registeredTH1F("hPed", "Pedestal Events", 250, -5, 20);
	hAvgEnergy = registeredTH1F("hAvgEnergy", "Corrected Smoothed Energy", nbinsx*5, emin, emax);
	
	double rebin = 10;
	hE8 = registeredTH2F("hE8","corrected energy spread",nbinsx/rebin,emin,emax,nbinsy,-wmax/3,wmax/3);
	pE8 = new TProfile("pE8","corrected energy spread",nbinsx/rebin,emin,emax);
	
	for(unsigned int t1=0; t1<nBetaTubes; t1++) {
		for(unsigned int t2=0; t2<nBetaTubes; t2++) {
			corrsProfEW[t1][t2] = new TProfile(("CorrEW_"+itos(t1)+"_"+itos(t2)).c_str(),"PMT Correlation",nbinsx/rebin,emin,emax);
			corrsProfEW[t1][t2]->GetXaxis()->SetTitle("Event Energy [keV]");
			corrsProfEW[t1][t2]->GetYaxis()->SetTitle("PMT deviation from mean [keV]");
			corrsProfEW[t1][t2]->GetYaxis()->SetTitleOffset(1.25);
		}
	}
	
	for(Side s = EAST; s <= WEST; ++s) {
		for(unsigned int t1=0; t1<nBetaTubes; t1++) {
			pEi[s][t1] = new TProfile((sideSubst("EnergyCorrector_%c",s)+itos(t1)).c_str(),"EnergyCorrector",nbinsx/rebin,emin,emax);
			addObject(pEi[s][t1]);
			for(unsigned int t2=0; t2<=t1; t2++) {
				corrs[s][t1][t2] = registeredTH2F(sideSubst("corr_%c_",s)+itos(t1)+"_"+itos(t2),"PMT Correlation",nbinsx/rebin,emin,emax,nbinsy,-wmax,wmax);
				corrsProf[s][t1][t2] = corrs[s][t1][t2]->ProfileX();
				corrsProf[s][t1][t2]->GetXaxis()->SetTitle("Event Energy [keV]");
				corrsProf[s][t1][t2]->GetYaxis()->SetTitle("PMT deviation from mean [keV]");
				corrsProf[s][t1][t2]->GetYaxis()->SetTitleOffset(1.25);
			}
			hLightBal[s][t1] = registeredTH1F(sideSubst("lightBal_%c_",s)+itos(t1),"LED light coupling",200,0,2);
		}
		for(unsigned int t=0; t<=nBetaTubes; t++) {
			hLowE[s][t][false] = registeredTH1F(sideSubst("lowE_%c_",s)+itos(t),"Low-energy events",120,-20,100);
			hLowE[s][t][false]->GetXaxis()->SetTitle("Event Energy [keV]");
			hLowE[s][t][true] = registeredTH1F(sideSubst("lowE_trig_%c_",s)+itos(t),"Triggered Low-energy events",120,-20,100);
			hLowE[s][t][true]->GetXaxis()->SetTitle("Event Energy [keV]");
		}
		h2foldProb[s][false] = registeredTH1F(sideSubst("p2fold_%c_",s),"Trigger Probability",50,0,1);
		h2foldProb[s][false]->GetXaxis()->SetTitle("expected trigger probability");
		h2foldProb[s][true] = registeredTH1F(sideSubst("p2fold_trig_%c_",s),"Trigger Probability",50,0,1);
		MLPdat[s]=NULL;
	}
}

std::vector<unsigned int> LEDAnalyzer::findJumps(ProcessedDataScanner& LSS, float emin, float emax, Side s) {
	assert(s<=WEST);
	printf("Locating LED scan transition points...\n");
	std::vector<unsigned int> jumps;
	bool isHi = false;
	LSS.startScan();
	while(LSS.nextPoint()) {
		if(isHi && LSS.scints[s].energy.x && LSS.scints[s].energy.x <= emin) {
			//printf("+Jump @%i (%.1f)",LSS.getCurrentEvent(),LSS.scints[s].energy.x);
			jumps.push_back(LSS.getCurrentEvent());
			isHi = false;
		}
		isHi |= LSS.scints[s].energy.x >= emax;
	}
	printf("Located %i LED scanner jumps in %.1f minutes (~%.1fs intervals)\n",
		   (int)jumps.size(),LSS.totalTime[BOTH]/60.,LSS.totalTime[BOTH]/jumps.size());
	return jumps;
}

void LEDAnalyzer::doWithAverage() {
	
	double E0 = Eavg->getAvgExcl();
	if(E0>emax) return;
	double E0dev = Eavg->getP()-E0;
	if(fabs(E0dev)>wmax/3.) return;
	
	if(pass==0) {
		// averaged fine correction profiles
		pE0->Fill(E0,E0dev);
		hPed->Fill(E0);
		for(Side s = EAST; s <= WEST; ++s) {
			for(unsigned int t=0; t<nBetaTubes; t++) {
				double tdev = Tavg[s][t]->getP()-Tavg[s][t]->getAvgExcl();
				if(fabs(tdev)>wmax) continue;
				pEi[s][t]->Fill(E0,tdev);
			}
		}
	} else if(pass==1) {
		// deviations, individual and correlated, from corrected energy
		double Efx = E0 + gE0->Eval(E0) - EZero;
		double E8dev = Eavg->getP()-(Efx + EZero);
		if(fabs(E8dev)>wmax/3.) return;
		hAvgEnergy->Fill(Efx);
		hE8->Fill(Efx,E8dev);
		pE8->Fill(Efx,E8dev);
		
		// PMT deviations from mean
		double tdev[2][nBetaTubes];
		for(Side s = EAST; s <= WEST; ++s)
			for(unsigned int t=0; t<nBetaTubes; t++)
				tdev[s][t] = Tavg[s][t]->getP() - (Tavg[s][t]->getAvgExcl()+gEi[s][t]->Eval(E0));
				
		// PMT width and same-side correlations
		for(Side s = EAST; s <= WEST; ++s) {
			for(unsigned int t1=0; t1<nBetaTubes; t1++) {
				if(fabs(tdev[s][t1])>wmax) continue;
				corrs[s][t1][t1]->Fill(Efx,tdev[s][t1]);
				corrsProf[s][t1][t1]->Fill(Efx,tdev[s][t1]);
				for(unsigned int t2=0; t2<t1; t2++) {
					if(fabs(tdev[s][t2])>wmax) continue;
					corrs[s][t1][t2]->Fill(Efx,tdev[s][t1]+tdev[s][t2]);
					corrsProf[s][t1][t2]->Fill(Efx,tdev[s][t1]+tdev[s][t2]);
				}
			}
		}
		
		// East-West PMT correlations
		for(unsigned int t1=0; t1<nBetaTubes; t1++)
			for(unsigned int t2=0; t2<nBetaTubes; t2++)
				corrsProfEW[t1][t2]->Fill(Efx,tdev[EAST][t1]+tdev[WEST][t2]);
				
	} else if(pass==2) {
		// example time-domain plots!
		gRawVT->SetPoint(nVT, clockAvg->getP(),Eavg->getP());
		gAvgVT->SetPoint(nVT,clockAvg->getP(),E0);
		gCorVT->SetPoint(nVT, clockAvg->getP(),E0 + gE0->Eval(E0));
		nVT++;
	}
	
}

void LEDAnalyzer::TimePlot(ProcessedDataScanner& LSS, unsigned int jn) {
	
	unsigned int pad = 5; // number of events to skip at edges of jumps
	unsigned int estart = jumps[jn]+pad;
	unsigned int eend = jumps[jn+1]-pad;
	unsigned int npts = eend-estart+1;

	gRawVT = new TGraph(npts);
	gAvgVT = new TGraph(npts);
	gCorVT = new TGraph(npts);
	nVT = 0;
	pass = 2;
	ScanSeg(LSS, estart, eend);
	while(gRawVT->GetN()>(int)nVT) {
		gRawVT->RemovePoint(nVT);
		gAvgVT->RemovePoint(nVT);
		gCorVT->RemovePoint(nVT);
	}
	
	gRawVT->Draw("AP");
	gRawVT->SetTitle("Smoothed average");
	gRawVT->GetXaxis()->SetTitle("Time [s]");
	gRawVT->GetYaxis()->SetTitle("Energy [keV]");
	gRawVT->GetYaxis()->SetTitleOffset(1.4);
	gRawVT->Draw("AP");
	gAvgVT->SetLineWidth(3);
	gAvgVT->SetLineStyle(2);
	gAvgVT->Draw("L");
	printCanvas("Segments/Segment_"+itos(jn)+"a");
	
	gRawVT->SetTitle("Corrected estimated LED output");
	gRawVT->Draw("AP");
	gCorVT->SetLineWidth(2);
	gCorVT->Draw("L");
	printCanvas("Segments/Segment_"+itos(jn)+"b");
}

void LEDAnalyzer::ScanSeg(ProcessedDataScanner& LSS, unsigned int eStart, unsigned int eEnd) {
	start();
	LSS.gotoEvent(eStart);
	do {
		// skip clipped events!
		bool isClipped = false;
		for(Side s=EAST; s<=WEST; ++s)
			for(unsigned int t=0; t<nBetaTubes; t++)
				isClipped |= LSS.scints[s].adc[t] > LSS.ActiveCal->getClipThreshold(s,t)-300;
		if(isClipped) continue;

		Eavg->x = LSS.getEnergy()*0.5;
		clockAvg->x = LSS.runClock[EAST];
		for(Side s = EAST; s <= WEST; ++s)
			for(unsigned int t=0; t<nBetaTubes; t++)
				Tavg[s][t]->x = LSS.scints[s].tuben[t].x;
		addPt();
	} while (LSS.nextPoint() && LSS.getCurrentEvent() <= eEnd);
	stop();
}

void LEDAnalyzer::dataPass(ProcessedDataScanner& LSS) {
	unsigned int pad = 5; // number of events to skip at edges of jumps
	for(unsigned int j=0; j<jumps.size()-1; j++) {
		if(jumps[j+1] < jumps[j]+1000) continue; // skip shortened segments
		ScanSeg(LSS,jumps[j]+pad,jumps[j+1]-pad);
	}
}

void LEDAnalyzer::ScanData(ProcessedDataScanner& LSS) {
	jumps = findJumps(LSS);
	if(jumps.size()<2) return;
	printf("Scanning data over %i intervals...\n",int(jumps.size())-1);
	
	// first pass to determine energy corrections
	pass = 0;
	dataPass(LSS);
	gE0 = TProf2TGraph(*pE0,10);
	for(Side s = EAST; s <= WEST; ++s)
		for(unsigned int t=0; t<nBetaTubes; t++)
			gEi[s][t] = TProf2TGraph(*pEi[s][t]);
	
	defaultCanvas->cd();
	defaultCanvas->SetLeftMargin(0.12);
	defaultCanvas->SetRightMargin(0.03);
	
	// find "true 0"
	EZero = hPed->GetBinCenter(hPed->GetMaximumBin());
	TF1 g1("g1","gaus");
	hPed->Fit(&g1,"Q","",EZero-5./sqrt(npts),EZero+5./sqrt(npts));
	EZero = g1.GetParameter(1);
	hPed->Draw();
	printCanvas("AvgPedestal");
	
	// make some time plots to check
	TimePlot(LSS,1);
	TimePlot(LSS,10);
	TimePlot(LSS,100);
	
	// process data with corrected average energy
	pass = 1;
	dataPass(LSS);
	
	hAvgEnergy->Scale(1.0/hAvgEnergy->GetBinWidth(1)/1000.0);
	hAvgEnergy->GetYaxis()->SetTitle("kCounts/keV");
	hAvgEnergy->GetXaxis()->SetTitle("LED output [keV equivalent]");
	hAvgEnergy->SetMaximum(50);
	hAvgEnergy->Draw();
	printCanvas("AvgEnergy");
	
	hE8->Draw();
	pE8->Draw("Same");
	printCanvas("E8Spread");
	
	printf(" Done.\n\n");
}


void LEDAnalyzer::CalcTrigEffic() {
	defaultCanvas->cd();
	for(Side s = EAST; s <= WEST; ++s) {
		for(unsigned int t = 0; t <= nBetaTubes; t++) {
			hLowE[s][t][false]->Draw();
			hLowE[s][t][true]->Draw("Same");
			printCanvas(sideSubst("TrigEffic/LowE_%c",s)+itos(t));
			TGraphAsymmErrors gEffic(hLowE[s][t][false]->GetNbinsX());
			gEffic.BayesDivide(hLowE[s][t][true],hLowE[s][t][false],"w");
			gEffic.SetMinimum(-0.10);
			gEffic.SetMaximum(1.10);
			gEffic.Draw("AP");
			gEffic.SetTitle("Trigger Efficiency");
			gEffic.GetXaxis()->SetTitle("Energy [keV]");
			gEffic.GetYaxis()->SetTitle("Efficiency");
			gEffic.Draw("AP");
			printCanvas(sideSubst("TrigEffic/Trigger_Efficiency_%c",s)+itos(t));
		}
		defaultCanvas->SetLogy(true);
		h2foldProb[s][true]->SetLineStyle(2);
		h2foldProb[s][false]->Draw();
		h2foldProb[s][true]->Draw("Same");
		printCanvas(sideSubst("TrigEffic/P2fold_%c",s));
		defaultCanvas->SetLogy(false);
		TGraphAsymmErrors gEffic(h2foldProb[s][false]->GetNbinsX());
		gEffic.BayesDivide(h2foldProb[s][true],h2foldProb[s][false],"w");
		gEffic.SetMinimum(0);
		gEffic.SetMaximum(1);
		gEffic.Draw("AP");
		gEffic.SetTitle("Trigger Efficiency");
		gEffic.GetXaxis()->SetTitle("expected 2fold efficiency");
		gEffic.GetXaxis()->SetRangeUser(0,1);
		gEffic.GetYaxis()->SetTitle("observed 2fold efficiency");
		gEffic.GetYaxis()->SetTitleOffset(1.25);
		gEffic.Draw("AP");
		printCanvas(sideSubst("TrigEffic/P2fold_Effic_%c",s));
	}
}

void LEDAnalyzer::CalcLightBal() {
	defaultCanvas->cd();
	for(Side s = EAST; s <= WEST; ++s) {
		Stringmap m;
		m.insert("side",sideSubst("%c",s));
		for(unsigned int t = 0; t < nBetaTubes; t++) {
			hLightBal[s][t]->SetLineStyle(t+1);
			float c = hLightBal[s][t]->GetBinCenter(hLightBal[s][t]->GetMaximumBin());
			TF1 g1("g1","gaus");
			g1.SetLineColor(2);
			hLightBal[s][t]->Fit(&g1,"Q","",c-0.15,c+0.15);
			m.insert(itos(t),g1.GetParameter(1));
		}
		for(unsigned int t = 0; t < nBetaTubes; t++)
			hLightBal[s][t]->Draw(t?"Same":"");
		printCanvas(sideSubst("LightBal_%c",s));
		qOut.insert("lightBal",m);
	}
}

TGraphErrors* TProfWidthSquared(const TProfile& P, unsigned int minpts = 10) {
	TGraphErrors* g = new TGraphErrors(P.GetNbinsX()-2);
	unsigned int ig = 0;
	for(int i=0; i<P.GetNbinsX()-2; i++) {
		if(P.GetBinEntries(i+1)<minpts) continue;
		double w = P.GetBinError(i+1);
		double N = P.GetBinEntries(i+1);
		g->SetPoint(ig,P.GetBinCenter(i+1),w*w*N);
		g->SetPointError(ig,0.0,2*w*w*sqrt(N));
		ig++;
	}
	while(g->GetN()>(int)ig) g->RemovePoint(ig);
	return g;
}

double CalcCorrelation(TGraph* wi, TGraph* wj, TGraph* wij, double x) {
	double si = wi->Eval(x);
	double sj = wj->Eval(x);
	return (wij->Eval(x)-si-sj)/(2*sqrt(si*sj));
}

TGraphErrors* CorrelationGraph(TGraphErrors* wi, TGraphErrors* wj, TGraphErrors* wij) {
	TGraphErrors* g = new TGraphErrors(wij->GetN());
	for(int i=0; i<wij->GetN(); i++) {
		double x,si,sj,sij;
		wi->GetPoint(i,x,si);
		wj->GetPoint(i,x,sj);
		wij->GetPoint(i,x,sij);
		g->SetPoint(i,x,(sij-si-sj)*0.5); // /(2*sqrt(si*sj)));
		
		double ei = wi->GetErrorY(i);
		double ej = wj->GetErrorY(i);
		//double eij = wij->GetErrorY(i);
		g->SetPointError(i,0,sqrt(ei*ei+ej*ej)*0.5);
	}
	return g;
}


void LEDAnalyzer::CalcCorrelations() {

	defaultCanvas->SetLogz();
	
	pE0->Draw();
	printCanvas("EnergyCorrection");
	
	defaultCanvas->SetLeftMargin(0.15);
	
	TGraphErrors* gE8W = TProfWidthSquared(*pE8);
	
	TF1 lineFit("lineFit","pol1",50,emax);
	gE8W->Fit(&lineFit,"R");
	
	gE8W->Draw("APZ");
	gE8W->SetMarkerStyle(7);
	gE8W->SetTitle("Spread in E_{8} from LED scan");
	gE8W->GetXaxis()->SetTitle("Energy [keV]");
	gE8W->GetYaxis()->SetTitle("Width^{2} [keV^{2}]");
	gE8W->GetYaxis()->SetTitleOffset(1.75);
	gE8W->Draw("APZ");
	//lineFit.SetLineStyle(2);
	//lineFit.SetRange(0,emax);
	//lineFit.Draw("Same");
	printCanvas("E8width");
	
	double w20 = gE8W->Eval(0);
	double w2l = lineFit.Eval(0);
	printf("Ped width kev^2: %.1f; Extrapolated 0 width kev^2: %.1f\n",w20,w2l);
	
	// PMT width (and PMT sum width) graphs
	TGraphErrors* gCP[2][nBetaTubes][nBetaTubes];
	TGraphErrors* gCPEW[nBetaTubes][nBetaTubes];
	for(Side s = EAST; s <= WEST; ++s) {
		for(unsigned int t=0; t<nBetaTubes; t++) {
			gCP[s][t][t] = TProfWidthSquared(*corrsProf[s][t][t]);
			lineFit.SetLineStyle(1);
			lineFit.SetRange(50,emax);
			gCP[s][t][t]->Fit(&lineFit,"QR");
			printf("...%i... Ped width kev^2: %.1f; Extrapolated 0 width kev^2: %.1f\n",t,gCP[s][t][t]->Eval(0),lineFit.Eval(0));
			gCP[s][t][t]->Draw("APZ");
			gCP[s][t][t]->SetMarkerStyle(7);
			gCP[s][t][t]->SetTitle((sideSubst("PMT %c",s)+itos(t+1)).c_str());
			gCP[s][t][t]->GetXaxis()->SetTitle("Energy [keV]");
			gCP[s][t][t]->GetYaxis()->SetTitle("Width^{2} [keV^{2}]");
			gCP[s][t][t]->GetYaxis()->SetTitleOffset(1.8);
			gCP[s][t][t]->Draw("APZ");
			printCanvas(sideSubst("Width_%c",s)+itos(t+1));
		}
	}
	for(unsigned int t1=0; t1<nBetaTubes; t1++) {
		for(unsigned int t2=0; t2<nBetaTubes; t2++) {
			gCPEW[t1][t2] = TProfWidthSquared(*corrsProfEW[t1][t2]);
		}
	}
	
	// pedestal correlations
	printf("\n\n************ Pedestal Correlations ****************\n");
	for(Side s = EAST; s <= WEST; ++s) {
		for(unsigned int t1=0; t1<nBetaTubes; t1++) {
			for(unsigned int t2=0; t2<t1; t2++) {
				gCP[s][t1][t2] = TProfWidthSquared(*corrsProf[s][t1][t2]);
				double s2i = gCP[s][t1][t1]->Eval(0);
				double s2j = gCP[s][t2][t2]->Eval(0);
				double s2ij = gCP[s][t1][t2]->Eval(0);
				double c = (s2ij-s2i-s2j)/(2*sqrt(s2i*s2j));
				printf("%c %i,%i %.2f,%.2f,%.2f c = %.2f\n",sideNames(s),t1,t2,s2i,s2j,s2ij,c);
			}
		}
	}
	printf("************ Cross-Side Correlations ****************\n");
	for(unsigned int t1=0; t1<nBetaTubes; t1++) {
		for(unsigned int t2=0; t2<nBetaTubes; t2++) {
			double s2i = gCP[EAST][t1][t1]->Eval(0);
			double s2j = gCP[WEST][t2][t2]->Eval(0);
			double s2ij = gCPEW[t1][t2]->Eval(0);
			double c = (s2ij-s2i-s2j)/(2*sqrt(s2i*s2j));
			printf("E%i:W%i %.2f,%.2f,%.2f c = %.2f\n",t1,t2,s2i,s2j,s2ij,c);
		}
	}
	
	// full correlations, and average
	TGraphErrors* gAvg = NULL;
	for(Side s = EAST; s <= WEST; ++s) {
		for(unsigned int t1=0; t1<nBetaTubes; t1++) {
			for(unsigned int t2=0; t2<t1; t2++) {
				TGraphErrors* g = CorrelationGraph(gCP[s][t1][t1],gCP[s][t2][t2],gCP[s][t1][t2]);
				g->Draw("AP");
				g->SetTitle((sideSubst("PMT %c",s)+itos(t1)+":"+itos(t2)).c_str());
				g->SetMarkerStyle(7);
				g->GetXaxis()->SetTitle("Energy [keV]");
				//g->GetYaxis()->SetTitle("Correlation");
				g->GetYaxis()->SetTitle("Correlated Width^{2} [keV^{2}]");
				g->GetYaxis()->SetTitleOffset(1.4);
				g->GetYaxis()->SetRangeUser(0,800);
				g->Draw("AP");
				if(!gAvg) gAvg = g;
				else accumPoints(*gAvg,*g);
				printCanvas(sideSubst("Correlations/%c_",s)+itos(t1)+itos(t2));
			}
		}
	}
	for(unsigned int t1=0; t1<nBetaTubes; t1++) {
		for(unsigned int t2=0; t2<nBetaTubes; t2++) {
			TGraphErrors* g = CorrelationGraph(gCP[EAST][t1][t1],gCP[WEST][t2][t2],gCPEW[t1][t2]);
			g->Draw("AP");
			g->SetTitle(("PMT E"+itos(t1)+":W"+itos(t2)).c_str());
			g->SetMarkerStyle(7);
			g->GetXaxis()->SetTitle("Energy [keV]");
			//g->GetYaxis()->SetTitle("Correlation");
			g->GetYaxis()->SetTitle("Correlated Width^{2} [keV^{2}]");
			g->GetYaxis()->SetTitleOffset(1.5);
			g->GetYaxis()->SetRangeUser(0,400);
			g->Draw("AP");
			if(!gAvg) gAvg = g;
			else accumPoints(*gAvg,*g);
			printCanvas("Correlations/Cross_E"+itos(t1)+"_W"+itos(t2));
			accumPoints(*gAvg,*g);
		}
	}
	gAvg->Draw("AP");
	gAvg->SetTitle("average correlated width");
	gAvg->SetMarkerStyle(7);
	gAvg->GetXaxis()->SetTitle("Energy [keV]");
	gAvg->GetYaxis()->SetTitle("Correlated Width^{2} [keV^{2}]");
	gAvg->GetYaxis()->SetTitleOffset(1.5);
	gAvg->GetYaxis()->SetRangeUser(0,300);
	gAvg->Draw("APZ");
	printCanvas("Correlations/LED_Width_Avg");
	
	// fix LED width...
	scale(*gAvg,-1.0);
	accumPoints(*gE8W,*gAvg,false,true);
	lineFit.SetLineStyle(1);
	lineFit.SetLineWidth(2);
	gE8W->Fit(&lineFit,"R");
	gE8W->Draw("APZ");
	gE8W->SetMarkerStyle(7);
	gE8W->SetTitle("Spread in E_{8} from LED scan, minus LED width");
	gE8W->GetXaxis()->SetTitle("Energy [keV]");
	gE8W->GetYaxis()->SetTitle("Width^{2} [keV^{2}]");
	gE8W->GetYaxis()->SetTitleOffset(1.75);
	gE8W->Draw("APZ");
	//lineFit.SetLineStyle(2);
	//lineFit.SetRange(0,emax);
	//lineFit.Draw("Same");
	printCanvas("E8width_MinusLED");

	
	return;
	/*
	for(Side s = EAST; s <= WEST; ++s) {
		float s0[nBetaTubes];	//< zero-signal noise
		float k[nBetaTubes];	//< noise proportionality
				
				TF1 corrFit("corrFit","([0]*[1]*[2]+[5]*sqrt([3]*[4])*(x-[6])*(x-[6]>0))/sqrt(([1]*[1]+[3]*(x-[6])*(x-[6]>0))*([2]*[2]+[4]*(x-[6])*(x-[6]>0)))",-10.,emax);
				
				corrFit.SetParameter(0,0.5);
				corrFit.SetParameter(5,0.02);
				
				corrFit.FixParameter(1,s0[t1]);
				corrFit.FixParameter(2,s0[t2]);
				corrFit.FixParameter(6,0);
				
				corrFit.FixParameter(3,k[t1]);
				corrFit.FixParameter(4,k[t2]);
				corrFit.SetLineColor(2);
				grat.Fit(&corrFit,"BR");
				grat.SetMinimum(-0.2);
				grat.SetMaximum(0.5);
				grat.SetTitle((sideSubst("PMT correlation %c",s)+itos(t1+1)+":"+itos(t2+1)).c_str());
				grat.Draw("ALP");
				grat.GetXaxis()->SetTitle("Event Energy [keV]");
				grat.GetYaxis()->SetTitle("PMT Correlation");
				grat.GetYaxis()->SetTitleOffset(1.25);
				grat.Draw("ALP");
				printCanvas(sideSubst("CorRat_%c_",s)+itos(t1)+"+"+itos(t2));
				{
					Stringmap m;
					m.insert("t1",t1);
					m.insert("t2",t2);
					m.insert("side",sideSubst("%c",s));
					m.insert("c0",corrFit.GetParameter(0));
					m.insert("dc0",corrFit.GetParError(0));
					m.insert("cc",corrFit.GetParameter(5));
					m.insert("dcc",corrFit.GetParError(5));
					m.insert("E0",EZero);
					m.display();
					qOut.insert("PMTCorr",m);
				}
			}
		}
	}
	*/
}



		/*
		LSS.gotoEvent(jumps[j]+pad+avgW-1);
		for(unsigned int i=avgW; i<=jumps[j+1]-pad-(jumps[j]+pad)-avgW; i++) {
			LSS.nextPoint();
			for(Side s = EAST; s <= WEST; ++s) {
			
				// skip clipped events
				bool isClipped = false;
				for(unsigned int t=0; t<nBetaTubes; t++)
					isClipped |= LSS.scints[s].adc[t] > LSS.ActiveCal->getClipThreshold(s,t)-300;
				if(isClipped) continue;
				evtE[s] = LSS.scints[s].energy.x;
				for(unsigned int t=0; t<=nBetaTubes; t++) {
					float tubeE = (t==nBetaTubes ? evtE[s] : LSS.scints[s].tuben[t].x);
					hLowE[s][t][false]->Fill(tubeE);
					if(LSS.Sis00_2fold(s)) hLowE[s][t][true]->Fill(tubeE);
					if(t<nBetaTubes && evtE[s]>200)
						hLightBal[s][t]->Fill(tubeE/evtE[s]);
				}
				// trigger probs
				for(unsigned int t=0; t<4; t++) TP[s].tubeProbs[t] = LSS.probTrig(s,t);
				float p2fold = TP[s].calcProb();
				h2foldProb[s][false]->Fill(p2fold);
				is2fold[s] = LSS.Sis00_2fold(s);
				if(is2fold[s]) h2foldProb[s][true]->Fill(p2fold);
				if(MLPdat[s] && p2fold>0.001 && p2fold < 0.999) {
					TriggerProbMLP::condition(TP[s].tubeProbs);
					is2fold[s] -= 0.5;
					ewt = evtE[s]<15?0.1:1;
					MLPdat[s]->Fill();
				}
				
				// calculate deviation from average
				float Eavg = ledAvg[s][nBetaTubes][i];
				for(unsigned int t1=0; t1<nBetaTubes; t1++) {
					float t1dev = LSS.scints[s].tuben[t1].x-ledAvg[s][t1][i];
					if(fabs(t1dev)>wmax) continue;
					corrs[s][t1][t1]->Fill(Eavg,t1dev);
					corrsProf[s][t1][t1]->Fill(Eavg,t1dev);
					for(unsigned int t2=0; t2<t1; t2++) {
						float t2dev = LSS.scints[s].tuben[t2].x-ledAvg[s][t2][i];
						if(fabs(t2dev)>wmax) continue;
						corrs[s][t1][t2]->Fill(Eavg,t1dev+t2dev);
						corrsProf[s][t1][t2]->Fill(Eavg,t1dev+t2dev);
					}
				}
			}		
		}
		*/

//-------------------------------------------------------

#include <TMultiLayerPerceptron.h>

void LEDAnalyzer::buildPerceptronTree() {
	for(Side s = EAST; s <= WEST; ++s) {
			MLPdat[s] = (TTree*)addObject(new TTree(sideSubst("MLPdat%c",s).c_str(),"Trigger efficiency learning set"));
			//MLPdat[s]->Branch("prob",&TP[s].tubeProbs,"p1/D:p2/D:p3/D:p4/D:p21/D:p31/D:p32/D:p41/D:p42/D:p43/D");
			MLPdat[s]->Branch("prob",&TP[s].tubeProbs,"p1/D:p2/D:p3/D:p4/D");
			MLPdat[s]->Branch("trig",&is2fold[s],"t/F");
			MLPdat[s]->Branch("w",&ewt,"w/F");
	}
}

void makeMLPfit(OutputManager& OM) {
	TFile f((OM.basePath+"/PMTCorr.root").c_str(),"UPDATE");
	for(Side s = EAST; s <= WEST; ++s) {
		TTree* MLPdat = (TTree*)f.Get((sideSubst("MLPdat%c",s)).c_str());
		//TMultiLayerPerceptron TMLP("p1,p2,p3,p4,p21,p31,p32,p41,p42,p43:12:4:12:t", MLPdat);
		TMultiLayerPerceptron TMLP("p1,p2,p3,p4:4:12:t", MLPdat);
		TMLP.SetEventWeight("w");
		TMLP.Train(55,"text");
		TMLP.Write(sideSubst("TrigMLP_%c",s).c_str());
	}
}

//-------------------------------------------------------

#include "G4toPMT.hh"

void spectrumGenTest() {
	G4toPMT G2P;
	G2P.addFile("/home/mmendenhall/geant4/output/LivPhys_495_neutronBetaUnpol_geomC/analyzed_1.root");
	RunNum rn = 14077;
	PMTCalibrator PCal(rn);
	G2P.setCalibrator(PCal);
	
	unsigned int nToSim = 1000000;
	G2P.resetSimCounters();
	G2P.startScan(true);
	G2P.setAFP(AFP_OFF);
	
	OutputManager OM("ThisNameIsNotUsedAnywhere",getEnvSafe("UCNA_ANA_PLOTS")+"/test/");
	std::vector<TH1*> hSpec;
	for(Side s=EAST; s<=WEST; ++s) {
		hSpec.push_back(OM.registeredTH1F(sideSubst("spec%c",s),"Energy Spectrum",100,0,1000));
		hSpec[s]->SetLineColor(2+2*s);
	}
	
	while(true) {
		G2P.nextPoint();
		for(Side s = EAST; s <= WEST; ++s) {
			EventType tp = G2P.fType;
			if(tp>=TYPE_IV_EVENT || !G2P.passesPositionCut(s) || G2P.fSide != s)
				continue;
			hSpec[G2P.fSide]->Fill(G2P.getErecon(),G2P.physicsWeight);
		}
		if(G2P.getCurrentEvent()>=nToSim)
			break;
	}
	
	drawSimulHistos(hSpec);
	OM.printCanvas("BetaSpectrum");
}

//-------------------------------------------------------


int main(int argc, char *argv[]) {

	if(true) {
		{
			LEDAnalyzer LA("PMTCorr",getEnvSafe("UCNA_ANA_PLOTS")+"/PMTCorrDatNew");
			LEDScanScanner LSS;
			std::vector<RunNum> bruns = CalDBSQL::getCDB()->findRuns("run_type = 'Asymmetry' AND gate_valve = 'Open'",16139,16216); // selectRuns(16193, 16216, "beta");
			//std::vector<RunNum> bruns = selectRuns(16097, 16216, "beta");
			LSS.addRuns(bruns);
			//LSS.addRun(16194);
			//LA.buildPerceptronTree();
			LA.ScanData(LSS);
			LA.CalcCorrelations();
			//LA.CalcTrigEffic();
			//LA.CalcLightBal();
			LA.write();
			//LA.setWriteRoot(true);
		}
		//OutputManager OM("PMTCorr",getEnvSafe("UCNA_ANA_PLOTS")+"/PMTCorrDatNew");
		//makeMLPfit(OM);
	}
	if(false) {
		LED2PMT L2P;

		// load MultiLayerPerceptron trigger prob fit
		TFile f((getEnvSafe("UCNA_ANA_PLOTS")+"/PMTCorrDatNew/PMTCorr.root").c_str(),"READ");
		TMultiLayerPerceptron* TMLP[2];
		TriggerProbMLP* TProb[2];
		for(Side s = EAST; s <= WEST; ++s) {
			TMLP[s] = (TMultiLayerPerceptron*)f.Get(sideSubst("TrigMLP_%c",s).c_str());
			assert(TMLP[s]);
			TMLP[s]->Print();
			TProb[s] = new TriggerProbMLP(TMLP[s]);
			L2P.PGen[s].setTriggerProb(TProb[s]);
		}

		L2P.nToSim = 60000;
		L2P.PGen[EAST].crosstalk = L2P.PGen[WEST].crosstalk  = 0;
		L2P.PGen[EAST].pedcorr = L2P.PGen[WEST].pedcorr  = 0.2;
		L2P.PGen[EAST].setLightbal(EAST, 0.928815, 0.859369, 1.20333, 1.03771);
		L2P.PGen[WEST].setLightbal(WEST, 0.935987, 1.04881, 1.07959, 0.941108);
		PMTCalibrator PCal(16194);
		L2P.setCalibrator(PCal);
		LEDAnalyzer LA("PMTCorrSim",getEnvSafe("UCNA_ANA_PLOTS")+"/PMTCorrSim_C0.2");
		for(unsigned int i=0; i<9; i++)
			LA.ScanData(L2P);
		LA.CalcCorrelations();
		LA.CalcTrigEffic();
		LA.CalcLightBal();
		LA.write();
	}

	return 0;
}

