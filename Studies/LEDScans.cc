#include "LEDScans.hh"
#include "UCNAException.hh"
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
		   (int)jumps.size(),totalTime.t[BOTH]/60.,totalTime.t[BOTH]/jumps.size());
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
				if(currentEvent && currentEvent <= eend) { 
					RWs[s][t].addCount(0,t==nBetaTubes?scints[s].energy.x:scints[s].tuben[t].x);
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
	
	unsigned int nbinsx = 200;
	float emax = 2000;
	unsigned int nbinsy = 200;
	float wmax = 1000;
	
	float emin = -emax/(2*nbinsx);
	emax += -emin;
	nbinsx++;
	
	// set up profiles
	TH2F* corrs[2][nBetaTubes][nBetaTubes];
	TProfile* corrsProf[2][nBetaTubes][nBetaTubes];
	for(Side s = EAST; s <= WEST; ++s) {
		for(unsigned int t1=0; t1<nBetaTubes; t1++) {
			corrs[s][t1][t1] = OM.registeredTH2F(sideSubst("PMT_%c_",s)+itos(t1),"PMT Width",nbinsx,emin,emax,nbinsy,-wmax,wmax);
			for(unsigned int t2=0; t2<t1; t2++)
				corrs[s][t1][t2] = OM.registeredTH2F(sideSubst("corr_%c_",s)+itos(t1)+"_"+itos(t2),"PMT Correlation",nbinsx,emin,emax,nbinsy,-wmax,wmax);
		}
	}
	
	// collect data
	std::vector<unsigned int> jumps = LSS.findJumps();
	if(jumps.size()<2) return;
	printf("Scanning data over %i intervals...\n",int(jumps.size())-1);
	unsigned int pad = 5;
	unsigned int w = 25;
	for(unsigned int j=0; j<jumps.size()-1; j++) {
		if(jumps[j+1] < jumps[j]+20) continue;
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
				// calculate deviation from average
				float Eavg = LSS.ledAvg[s][nBetaTubes][i];
				for(unsigned int t1=0; t1<nBetaTubes; t1++) {
					float t1e = LSS.scints[s].tuben[t1].x;
					float t1avg = LSS.ledAvg[s][t1][i];
					corrs[s][t1][t1]->Fill(Eavg,t1e-t1avg);
					for(unsigned int t2=0; t2<t1; t2++)
						corrs[s][t1][t2]->Fill(Eavg,t1e-t1avg+LSS.scints[s].tuben[t2].x-LSS.ledAvg[s][t2][i]);
				}
			}		
		}
	}
	printf(" Done.\n\n");
	
		
	// analyze and make plots
	OM.defaultCanvas->cd();
	OM.defaultCanvas->SetLogz();
	for(Side s = EAST; s <= WEST; ++s) {
		for(unsigned int t1=0; t1<nBetaTubes; t1++) {
			// individual PMT spreads
			corrs[s][t1][t1]->Draw("Col");
			corrsProf[s][t1][t1] = corrs[s][t1][t1]->ProfileX();
			corrsProf[s][t1][t1]->SetErrorOption("s");
			corrsProf[s][t1][t1]->SetLineColor(6);
			corrsProf[s][t1][t1]->Draw("Same");
			OM.printCanvas(sideSubst("PMT_%c",s)+itos(t1));
			
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
			TLatex lx;
			lx.SetTextColor(2);
			lx.SetTextSize(0.04);
			char tmp[1024];
			sprintf(tmp,"E=0 Width: %.1fkeV",sqrt(lineFit.GetParameter(0)));
			lx.DrawLatex(100,lineFit.Eval(1500),tmp);			
			OM.printCanvas(sideSubst("Width_%c",s)+itos(t1));
			
			// PMT correlation widths
			for(unsigned int t2=0; t2<t1; t2++) {
				TGraph gunc(nbinsx);
				gunc.SetLineColor(4);
				TGraph gcor(nbinsx);
				gcor.SetLineColor(2);
				TGraph grat(nbinsx);
				for(unsigned int i=0; i<nbinsx; i++) {
					corrsProf[s][t1][t2] = corrs[s][t1][t2]->ProfileX();
					corrsProf[s][t1][t2]->SetErrorOption("s");
					float e0 = corrsProf[s][t1][t1]->GetBinCenter(i+1);
					float xc1 = corrsProf[s][t1][t1]->GetBinError(i+1);
					float xc2 = corrsProf[s][t2][t2]->GetBinError(i+1);
					float xunc = pow(xc1,2.0)+pow(xc2,2.0);
					float xcor = pow(corrsProf[s][t1][t2]->GetBinError(i+1),2.0);
					gunc.SetPoint(i,e0,xunc);
					gcor.SetPoint(i,e0,xcor);
					grat.SetPoint(i,e0,(xcor-xunc)/(2*xc1*xc2));
				}
				gunc.Draw("ALP");
				gcor.Draw("LP");
				OM.printCanvas(sideSubst("Correlation_%c_",s)+itos(t1)+"+"+itos(t2));
				grat.Draw("ALP");
				OM.printCanvas(sideSubst("CorRat_%c_",s)+itos(t1)+"+"+itos(t2));
			}
		}
	}	
	
}
