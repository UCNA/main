#include "LEDScans.hh"
#include "SMExcept.hh"
#include "GraphicsUtils.hh"
#include <TProfile.h>
#include <climits>
#include <cmath>
#include <TLatex.h>
#include "RollingWindow.hh"

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
	
	for(Side s=EAST; s<=WEST; ++s) {
		// beta scintillators
		Tch->SetBranchAddress(sideSubst("Scint%c",s).c_str(),&scints[s]);
		
		// MWPC totals
		Tch->SetBranchAddress(sideSubst("Anode%c",s).c_str(),&anode[s]);
		Tch->SetBranchAddress(sideSubst("CathMax%c",s).c_str(),&cathMax[s]);
	}
}

std::vector<unsigned int> LEDScanScanner::findJumps(float emin, float emax, Side s) {
	assert(s<=WEST);
	printf("Locating LED scan transition points...\n");
	std::vector<unsigned int> jumps;
	bool isHi = false;
	startScan();
	while(nextPoint()) {
		if(isHi && scints[s].tuben[0].x && scints[s].tuben[0].x <= emin) {
			jumps.push_back(currentEvent);
			isHi = false;
		}
		isHi |= scints[s].tuben[0].x >= emax;
	}
	printf("Located %i LED scanner jumps in %.1f minutes (~%.1fs intervals)\n",
		   (int)jumps.size(),totalTime[BOTH]/60.,totalTime[BOTH]/jumps.size());
	return jumps;
}

void LEDScanScanner::analyzeSegment(unsigned int estart, unsigned int eend, unsigned int w) {
	assert(estart <= eend && eend < nEvents);
		
	// clear previous data; set up rolling windows
	std::vector<RollingWindow> RWs[2];
	for(Side s = EAST; s <= WEST; ++s) {
		for(unsigned int t=0; t<=nBetaTubes; t++) {
			ledAvg[s][t].clear();
			ledRMS[s][t].clear();
			RWs[s].push_back(RollingWindow(2*w+1));
		}
	}
	
	// pre-fill rolling window leading edge
	gotoEvent(estart);
	do {
		for(Side s = EAST; s <= WEST; ++s)
			for(unsigned int t=0; t<=nBetaTubes; t++)
				RWs[s][t].addCount(0,t==nBetaTubes?scints[s].energy.x:scints[s].tuben[t].x);
	} while (nextPoint() && currentEvent <= eend && currentEvent <= estart+w);
	
	// store rolling window average over event range
	for(unsigned int i=estart; i<=eend; i++) {
		for(Side s = EAST; s <= WEST; ++s) {
			for(unsigned int t=0; t<=nBetaTubes; t++) {
				float thisE = (t==nBetaTubes)?scints[s].energy.x:scints[s].tuben[t].x;
				if(currentEvent && currentEvent <= eend) { 
					RWs[s][t].addCount(0,thisE);
				} else if(eend-estart > w) {
					RWs[s][t].popExcess();
				}
				ledAvg[s][t].push_back(RWs[s][t].getAvg());
				ledRMS[s][t].push_back(RWs[s][t].getRMS());
			}
		}
		if(currentEvent && currentEvent <= eend) nextPoint();
	}
}


//----------------------------------------------------------------

void PMT_LED_Correlations(OutputManager& OM, LEDScanScanner& LSS) {
	
	unsigned int nbinsx = 500;
	unsigned int nbinsy = 200;
	float emax = 2000;
	float wmax = 500;
	float emin = -50;
	
	// set up profiles
	TH2F* corrs[2][nBetaTubes][nBetaTubes];
	TH1F* truezero[2][nBetaTubes+1];
	TProfile* corrsProf[2][nBetaTubes][nBetaTubes];
	for(Side s = EAST; s <= WEST; ++s) {
		for(unsigned int t1=0; t1<nBetaTubes; t1++) {			
			for(unsigned int t2=0; t2<=t1; t2++) {
				corrs[s][t1][t2] = OM.registeredTH2F(sideSubst("corr_%c_",s)+itos(t1)+"_"+itos(t2),"PMT Correlation",nbinsx,emin,emax,nbinsy,-wmax,wmax);
				corrsProf[s][t1][t2] = corrs[s][t1][t2]->ProfileX();
				corrsProf[s][t1][t2]->SetErrorOption("s");
			}
		}
		for(unsigned int t=0; t<=nBetaTubes; t++)
			truezero[s][t] = OM.registeredTH1F(sideSubst("zero_%c_",s)+itos(t),"Low-energy events",240,-20,100);
	}
	
	// collect data
	std::vector<unsigned int> jumps = LSS.findJumps();
	if(jumps.size()<2) return;
	printf("Scanning data over %i intervals...\n",int(jumps.size())-1);
	unsigned int pad = 5;
	unsigned int w = 50;
	for(unsigned int j=0; j<jumps.size()-1; j++) {
		if(jumps[j+1] < jumps[j]+2000) continue;
		LSS.analyzeSegment(jumps[j]+pad,jumps[j+1]-pad,w);
		LSS.gotoEvent(jumps[j]+pad+w-1);
		for(unsigned int i=w; i<=jumps[j+1]-jumps[j]-2*pad-w; i++) {
			LSS.nextPoint();
			for(Side s = EAST; s <= WEST; ++s) {
				// skip clipped events
				bool isClipped = false;
				for(unsigned int t=0; t<nBetaTubes; t++)
					isClipped |= LSS.scints[s].adc[t] > LSS.ActiveCal->getClipThreshold(s,t)-300;
				if(isClipped) continue;
				for(unsigned int t=0; t<=nBetaTubes; t++)
					truezero[s][t]->Fill(t==nBetaTubes?LSS.scints[s].energy.x:LSS.scints[s].tuben[t].x);
				// calculate deviation from average
				float Eavg = LSS.ledAvg[s][nBetaTubes][i];
				for(unsigned int t1=0; t1<nBetaTubes; t1++) {
					float t1dev = LSS.scints[s].tuben[t1].x-LSS.ledAvg[s][t1][i];
					if(fabs(t1dev)>1000) continue;
					corrs[s][t1][t1]->Fill(Eavg,t1dev);
					corrsProf[s][t1][t1]->Fill(Eavg,t1dev);
					for(unsigned int t2=0; t2<t1; t2++) {
						float t2dev = LSS.scints[s].tuben[t2].x-LSS.ledAvg[s][t2][i];
						if(fabs(t2dev)>1000) continue;
						corrs[s][t1][t2]->Fill(Eavg,t1dev+t2dev);
						corrsProf[s][t1][t2]->Fill(Eavg,t1dev+t2dev);
					}
				}
			}		
		}
	}
	printf(" Done.\n\n");
	
		
	// analyze and make plots
	TLatex lx;
	lx.SetTextColor(2);
	lx.SetTextSize(0.04);
	char tmp[1024];
	OM.defaultCanvas->cd();
	OM.defaultCanvas->SetLogz();
	for(Side s = EAST; s <= WEST; ++s) {
		float s0[nBetaTubes]; //< zero-signal noise
		float k[nBetaTubes]; //< noise proportionality
		float ezero = truezero[s][nBetaTubes]->GetBinCenter(truezero[s][nBetaTubes]->GetMaximumBin());
		truezero[s][nBetaTubes]->Draw();
		drawVLine(ezero, OM.defaultCanvas, 2);
		OM.printCanvas(sideSubst("LowEnergy_%c",s));
		for(unsigned int t1=0; t1<nBetaTubes; t1++) {
			// individual PMT spreads
			corrs[s][t1][t1]->Draw("Col");
			corrsProf[s][t1][t1]->SetLineColor(6);
			//corrsProf[s][t1][t1]->Draw("Same");
			OM.printCanvas(sideSubst("PMT_%c",s)+itos(t1));
			s0[t1] = corrsProf[s][t1][t1]->GetBinError(corrsProf[s][t1][t1]->FindBin(0.0));
			
			// individual PMT width graphs
			TGraph gtb(nbinsx+1);
			TF1 lineFit("lineFit","pol1",10,500);
			lineFit.SetLineColor(2);
			gtb.SetPoint(0,0,0);
			for(unsigned int i=1; i<=nbinsx; i++)
				gtb.SetPoint(i,corrsProf[s][t1][t1]->GetBinCenter(i),pow(corrsProf[s][t1][t1]->GetBinError(i),2.0));
			gtb.Fit(&lineFit,"R");
			gtb.SetTitle("PMT Width^2 vs Energy");
			gtb.Draw("ALP");			
			OM.printCanvas(sideSubst("Width_%c",s)+itos(t1));
			k[t1] = lineFit.GetParameter(1);
			
			// PMT correlation widths
			for(unsigned int t2=0; t2<t1; t2++) {
				TGraph gunc(nbinsx);
				gunc.SetLineColor(4);
				TGraph gcor(nbinsx);
				gcor.SetLineColor(2);
				TGraph grat(nbinsx);
				for(unsigned int i=0; i<nbinsx; i++) {
					float e0 = corrsProf[s][t1][t1]->GetBinCenter(i+1);
					float xc1 = corrsProf[s][t1][t1]->GetBinError(i+1);
					float xc2 = corrsProf[s][t2][t2]->GetBinError(i+1);
					float xunc = pow(xc1,2.0)+pow(xc2,2.0);
					float xcor = pow(corrsProf[s][t1][t2]->GetBinError(i+1),2.0);
					gunc.SetPoint(i,e0,xunc);
					gcor.SetPoint(i,e0,xcor);
					if(xcor>4 && xunc>4)
						grat.SetPoint(i,e0,(xcor-xunc)/(2*xc1*xc2));
					else
						grat.SetPoint(i,0,0);
				}
				gunc.Draw("ALP");
				gcor.Draw("LP");
				OM.printCanvas(sideSubst("Correlation_%c_",s)+itos(t1)+"+"+itos(t2));
				
				TF1 corrFit("corrFit","([0]*[1]*[2]+[5]*sqrt([3]*[4])*(x-[6]))/sqrt(([1]*[1]+[3]*(x-[6]))*([2]*[2]+[4]*(x-[6])))",ezero,emax);
				
				corrFit.SetParameter(0,0.5);
				corrFit.SetParameter(5,0.02);
				
				corrFit.FixParameter(1,s0[t1]);
				corrFit.FixParameter(2,s0[t2]);
				corrFit.FixParameter(6,ezero);
				
				corrFit.FixParameter(3,k[t1]);
				corrFit.FixParameter(4,k[t2]);
				corrFit.SetLineColor(2);
				grat.Fit(&corrFit,"BR");
				grat.SetMinimum(0);
				grat.SetMaximum(1.0);
				grat.Draw("ALP");
				sprintf(tmp,"c1=%.3f(%.3f), c2=%.3f(%.3f)",corrFit.GetParameter(0),corrFit.GetParError(0),corrFit.GetParameter(5),corrFit.GetParError(5));
				lx.DrawLatex(100,0.5,tmp);
				OM.printCanvas(sideSubst("CorRat_%c_",s)+itos(t1)+"+"+itos(t2));
			}
		}
	}	
	
}

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
			hSpec[G2P.fSide]->Fill(G2P.getEtrue(),G2P.physicsWeight);
		}
		if(G2P.nSimmed>=nToSim)
			break;
	}
	
	drawSimulHistos(hSpec);
	OM.printCanvas("BetaSpectrum");
}

