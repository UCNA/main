#include "PeakTracker.hh"
#include "GraphUtils.hh"
#include "BetaScint.hh"
#include "RunManager.hh"
#include "CoTracker.hh"
#include "AnalysisDB.hh"
#include "GraphicsUtils.hh"
#include "PathUtils.hh"
#include <TProfile.h>
#include <TH1D.h>
#include "QFile.hh"
#include "TSpectrumUtils.hh"
#include "PostOfficialAnalyzer.hh"
#include "EnergyCalibrator.hh"

std::vector<TGraphErrors*> mergeloader(const std::vector<RunNum>& druns, std::string sname, std::string infotype, bool domerged) {
	int tstart = AnalysisDB::ADB().startTime(druns.front());
	std::vector<TGraphErrors*> gs;
	std::vector<int> runstarts;
	unsigned int nruns = druns.size();
	for(unsigned int n=0; n<nruns; n++) {
		TGraphErrors* g = CalDBSQL::getCDB()->getRunMonitor(druns[n],sname,infotype);
		if(g) {
			gs.push_back(g);
			runstarts.push_back(AnalysisDB::ADB().startTime(druns[n],tstart));
		} else {
			printf("mergeloader: Missing data for %s in run %i!\n",sname.c_str(),druns[n]);
		}
	}
	if(domerged)
		gs.push_back(merge_plots(gs,runstarts));
	return gs;
}

void trackPeaks(RunNum startRun, RunNum endRun, bool makePlots, float tsplit) {
	
	// select runs with known time stamps
	char tmp[1024];
	sprintf(tmp,"run_number >= %i AND run_number <= %i AND start_time IS NOT null AND end_time IS NOT null ORDER BY run_number ASC",startRun,endRun);
	std::vector<RunNum> druns = AnalysisDB::ADB().findRuns(tmp);
	unsigned int nruns = druns.size();
	printf("%i runs found.\n",nruns);
	if(!nruns)
		return;
	
	// set plot paths, output canvas
	char plotpath[1024];
	sprintf(plotpath,"../PostPlots/PeakShifts/R_%i_%i/",druns.front(),druns.back());
	makePath(plotpath);
	TCanvas defaultCanvas;
	defaultCanvas.SetFillColor(0);
	defaultCanvas.SetCanvasSize(300,300);
	
	// collect and process data
	for(Side s = EAST; s <= WEST; s = nextSide(s) ) {
		
		printf("Collecting data...\n");
		
		// sensor names to track
		std::vector<std::string> snames;
		if(s == EAST) {
			snames.push_back("ADCRefELED");
			snames.push_back("ADCE3Beta");
			snames.push_back("ADCE4Beta");
			snames.push_back("ADCE1Beta");
			snames.push_back("ADCE2Beta");
		} else {
			snames.push_back("ADCRefWLED");
			snames.push_back("ADCW1Beta");
			snames.push_back("ADCW2Beta");
			snames.push_back("ADCW3Beta");
			snames.push_back("ADCW4Beta");
		}
		
		// graphs from each run
		std::vector<TGraphErrors*> ledGraphs[5];
		std::vector<TGraphErrors*> pedGraphs[5];
		std::vector<TGraphErrors*> ledCombo;
		std::vector<TGraphErrors*> pedCombo;
		
		for(unsigned int t=0; t<5; t++) {
			ledGraphs[t] = mergeloader(druns, snames[t], "GMS_peak", makePlots);
			if(makePlots) {
				pedGraphs[t] = mergeloader(druns, snames[t], "pedestal", makePlots);
				ledCombo.push_back(ledGraphs[t].back());
				pedCombo.push_back(pedGraphs[t].back());
			}
		}
		
		printf("\tDone.\n");
		
		if(makePlots) {
			printf("Making plots...\n");
			sprintf(tmp,"%s/LED_History_%c.pdf",plotpath,sideNames(s));
			drawTogether(ledCombo, 0, 4000, &defaultCanvas, tmp, "LED Peak History");
			sprintf(tmp,"%s/Ped_History_%c.pdf",plotpath,sideNames(s));
			drawTogether(pedCombo, 0, 4000, &defaultCanvas, tmp, "Pedestal History");
		} else {
			std::vector<RunNum> refBreaks;
			std::vector<int> runstarts;
			for(unsigned int n=0; n<nruns; n++)
				runstarts.push_back(AnalysisDB::ADB().startTime(druns[n]));
			std::vector<int> runends;
			for(unsigned int n=0; n<nruns; n++)
				runends.push_back(AnalysisDB::ADB().endTime(druns[n]));
			printf("\nScanning %i runs...\n",nruns);
			double x0,y0,x1,y1;
			for(unsigned int n=0; n<nruns-1; n++) {
				int nup = 0;
				int ndown = 0;
				bool refjump = false;
				for(unsigned int t=0; t<5; t++) {
					ledGraphs[t][n]->GetPoint(ledGraphs[t][n]->GetN()-1,x0,y0);
					ledGraphs[t][n+1]->GetPoint(0,x1,y1);
					if( (x1-x0)/x0 > 0.02 )
						nup++;
					if( (x1-x0)/x0 < -0.02 )
						ndown++;
					if(t==0 && fabs((x1-x0)/x0) > 0.03)
						refjump = true;
				}				 
				if( (refjump && nup<5 && ndown<5) || runstarts[n+1]-runends[n] > tsplit) {
					if(runstarts[n+1]-runends[n] > tsplit)
						printf("Split at %i due to time gap %i\n",n,runstarts[n+1]-runends[n]);
					if(refjump && nup<5 && ndown<5)
						printf("Split at %i due to reference jump",n);
					refBreaks.push_back(n);
				}
			}
			printf("Generating breaks file...\n");
			refBreaks.push_back(nruns-1);
			sprintf(tmp,"../SummaryData/RefBreaks_%c.txt",sideNames(s));
			QFile fbreak(tmp,false);
			RunNum rn0 = druns[0];
			for(unsigned int i=0; i<refBreaks.size(); i++) {
				Stringmap m;
				m.insert("run_start",itos(rn0));
				m.insert("run_end",itos(druns[refBreaks[i]]));
				fbreak.insert("ContigRange",m);
				if(refBreaks[i]+1<druns.size())
					rn0 = druns[refBreaks[i]+1];
			}
			
			fbreak.commit();
		}
	}
}


void trackGMS(RunNum startRun, RunNum endRun) {
	
	assert(false); //TODO make this work
	
	// select runs with known time stamps (e.g. GMS correctible)
	char tmp[1024];
	sprintf(tmp,"run_number >= %i AND run_number <= %i AND start_time IS NOT null AND end_time IS NOT null ORDER BY run_number ASC",startRun,endRun);
	std::vector<RunNum> druns = AnalysisDB::ADB().findRuns(tmp);
	unsigned int nruns = druns.size();
	printf("%i runs found.\n",nruns);
	if(!nruns)
		return;
	
	// set plot paths, output canvas
	char plotpath[1024];
	sprintf(plotpath,"../PostPlots/GMS/GMS_%i_%i/",druns.front(),druns.back());
	makePath(plotpath);
	TCanvas defaultCanvas;
	defaultCanvas.SetFillColor(0);
	defaultCanvas.SetCanvasSize(300,300);
	
	std::vector<float> time;
	std::vector<float> gmsCor[2][nBetaTubes];
	std::vector<float> ledBright[2];
		
	float r0_start = AnalysisDB::ADB().startTime(druns[0]);
	
	// collect data
	for(unsigned int r=0; r<nruns; r++) {
		PMTCalibrator PCal(druns[r],CalDBSQL::getCDB());
		PCal.printSummary();
		float rstart = AnalysisDB::ADB().startTime(druns[r]);
		float rend = AnalysisDB::ADB().endTime(druns[r]);
		unsigned int ndivs = (unsigned int)((rend-rstart)/600.0);
		if(ndivs<1)
			ndivs = 1;
		printf("Run %i in %i segments\n",druns[r],ndivs);
		for(unsigned int n=0; n<=ndivs; n++) {
			float t0 = (rend-rstart)*float(n)/float(ndivs);
			time.push_back((rstart-r0_start+t0)/3600);
			// TODO
			/*
			for(Side s = EAST; s <= WEST; s = nextSide(s) ) {
				for(unsigned int t=0; t<nBetaTubes; t++)
					gmsCor[s][t].push_back(PCal.gmsFactor(s, t, t0)/PCal.getGMS0(s,t));
				ledBright[s].push_back(PCal.ledBrightness(s, t0));
			}
			 */
		}
	}
	
	// make plots
	defaultCanvas.SetLogy(true);
	for(Side s = EAST; s <= WEST; s = nextSide(s) ) {
		
		for(unsigned int t=0; t<nBetaTubes; t++) {
			TGraph* gGMS = new TGraph(time.size());
			gGMS->SetMarkerColor(2+t);
			gGMS->SetLineColor(2+t);
			gGMS->SetMinimum(0.5);
			gGMS->SetMaximum(1.5);
			sprintf(tmp,"GMS Corrections History, Runs %i-%i",druns.front(),druns.back());
			gGMS->SetTitle(tmp);
			for(unsigned int i=0; i<time.size(); i++)
				gGMS->SetPoint(i,time[i],gmsCor[s][t][i]);
			if(t==0) {
				gGMS->Draw("APL");
				gGMS->GetXaxis()->SetTitle("Time [Hours]");
				gGMS->GetYaxis()->SetTitle("GMS Correction Factor");
				gGMS->Draw("APL");
			} else
				gGMS->Draw("PL");
		}
		
		sprintf(tmp,"%s/GMS_Tube_History_%c.pdf",plotpath,sideNames(s));
		defaultCanvas.Print(tmp);
		
		TGraph* gLED = new TGraph(time.size());
		gLED->SetMarkerColor(1);
		gLED->SetMinimum(0.05);
		gLED->SetMaximum(1.0);
		sprintf(tmp,"LED Brightness History, Runs %i-%i",druns.front(),druns.back());
		gLED->SetTitle(tmp);
		for(unsigned int i=0; i<time.size(); i++)
			gLED->SetPoint(i,time[i],ledBright[s][i]);
		gLED->Draw("APL");
		gLED->GetXaxis()->SetTitle("Time [Hours]");
		gLED->GetYaxis()->SetTitle("LED brightness (relative to Co60 source)");
		gLED->Draw("APL");
		
		sprintf(tmp,"%s/GMS_LED_History_%c.pdf",plotpath,sideNames(s));
		defaultCanvas.Print(tmp);
	}
}



// data structure for a TGraphErrors point
struct tgePoint {
	tgePoint(float_err a, float_err b): x(a), y(b) {}
	float_err x;
	float_err y;
	void setGraph(TGraphErrors& tg, unsigned int n) const {
		tg.SetPoint(n,x.x,y.x);
		tg.SetPointError(n,x.err,y.err);
	}
};

void trackCo60(RunNum startRun0, RunNum endRun0, Side s) {
	
	if(s==BOTH) {
		trackCo60(startRun0, endRun0, EAST);
		trackCo60(startRun0, endRun0, WEST);
		return;
	} else if (s==NONE) {
		return;
	}
	
	
	TCanvas* defaultCanvas = new TCanvas();
	defaultCanvas->SetFillColor(0);
	defaultCanvas->SetCanvasSize(300,300);
	
	char tmp[1024];
	
	// input ranges file
	sprintf(tmp,"../SummaryData/RefBreaks_%c.txt",sideNames(s));
	QFile fbreak(tmp,true);
	std::vector<Stringmap> ranges = fbreak.retrieve("ContigRange");
	
	// process each range
	for(unsigned int nrange = 0; nrange < ranges.size(); nrange++) {
		
		// select runs with known time stamps
		RunNum startRun = (RunNum)ranges[nrange].getDefault("run_start",0.0);
		RunNum endRun = (RunNum)ranges[nrange].getDefault("run_end",0.0);
		if(endRun < startRun0 || startRun > endRun0)
			continue;
		sprintf(tmp,"run_number >= %i AND run_number <= %i AND start_time IS NOT null AND end_time IS NOT null ORDER BY run_number ASC",startRun,endRun);
		std::vector<RunNum> druns = AnalysisDB::ADB().findRuns(tmp);
		unsigned int nruns = druns.size();
		printf("%i...%i: %i runs.\n",startRun,endRun,nruns);
		char plotpath[1024];
		sprintf(plotpath,"../PostPlots/PeakShifts/Co60_%c_%i_%i/",sideNames(s),startRun,endRun);
		makePath(plotpath);
		
		float segmentSize = 40000;
		if(s==WEST) {
			segmentSize = 100000;
			if( startRun == 9409 || startRun == 9490 || startRun == 9589 || startRun == 9717 )
				segmentSize = 400000;
			if( startRun == 9766 )
				segmentSize = 1000000;
		}
		
		
		// load run data from each run
		std::vector<float> adc;
		std::vector<float> time;
		for(unsigned int i=0; i<nruns; ++i) {
			unsigned int tstart = AnalysisDB::ADB().startTime(druns[i]);
			RunManager TM(druns[i],RUNMODE_FULLID);
			CoTracker CT(&TM,s);
			for(unsigned int e=0; e<CT.nEvents; e++) {
				if(CT.gmsCo(e) && !CT.isCrud(e)) {
					adc.push_back(CT.f_refadc[e]);
					time.push_back(CT.f_rclock[e]+tstart);
				}
			}
		}
		
		
		// combined data fit for peak estimates
		float plotmin = 1200;
		float plotmax = 3500;
		float pkMin,pkMax,sigma;
		if(s==EAST) {
			pkMin = 1200;
			pkMax = 2400;
			sigma = 25;
			
		} else {
			pkMin = 1900;
			pkMax = 3200;
			sigma = 30;
			if(startRun == 9589)
				sigma = 20;
		}
		
		if(startRun >= 12229) {
			plotmin = 500;
			plotmax = 2000;
			pkMin = 800;
			pkMax = 1400;
		}
		
		TH1F hRefCoCombined("gmsRefCo","GMS Co60 Reference",200,plotmin,plotmax);
		hRefCoCombined.Sumw2();
		hRefCoCombined.SetLineColor(4);
		for(unsigned int e = 0; e< adc.size(); e++)
			hRefCoCombined.Fill(adc[e]);
		std::vector<SpectrumPeak> expectedPeaks;
		expectedPeaks.push_back(SpectrumPeak(REF_CO60_2));
		expectedPeaks.push_back(SpectrumPeak(REF_CO60_1));
		defaultCanvas->cd();
		defaultCanvas->SetLogy(true);
		sprintf(tmp,"%s/Co60_%05i-%05i_%c.pdf",plotpath,startRun,endRun,sideNames(s));
		std::vector<SpectrumPeak> foundPeaks = fancyMultiFit(&hRefCoCombined, sigma, expectedPeaks, true, tmp, 1.5, pkMin, pkMax);
		defaultCanvas->SetLogy(false);
		if(foundPeaks.size()!=2) {
			printf("**** Warning: %i %c %i/2 Co60 peaks found! Aborting\n",startRun,sideNames(s),(int)foundPeaks.size());
			continue;
		}
		
		// divide and conquer!
		unsigned int nsegs = (unsigned int)(adc.size()/segmentSize);
		if(!nsegs) nsegs = 1;
		printf("Dividing %i events into %i segments.\n",(int)adc.size(),nsegs);
		sprintf(tmp,"%s/PeakData_%c.txt",plotpath,sideNames(s));
		FILE* peaksout = fopen(tmp,"w");
		sprintf(tmp,"%s/_ErrorNotes.txt",plotpath);
		FILE* errout = fopen(tmp,"w");
		TProfile segTimes("segTimes","Segment Timing",nsegs,0,nsegs);
		std::vector<tgePoint> centers[2];
		std::vector<tgePoint> widths[2];
		
		for(unsigned int n=0; n<nsegs; n++) {
			
			
			// fill spectrum histogram
			sprintf(tmp,"gmsRefCo_Seg_%i",n);
			TH1F hRefCo(tmp,"GMS Co60 Reference",200,plotmin,plotmax);
			hRefCo.Sumw2();
			hRefCo.SetLineColor(4);
			for(unsigned int e = (unsigned int)(adc.size()*(float(n)/float(nsegs))); e<(unsigned int)(adc.size()*(float(n+1)/float(nsegs))); e++) {
				hRefCo.Fill(adc[e]);
				segTimes.Fill(float(n)+0.5,time[e]);
			}
			
			// bacground subtract
			TH1F* bgspec = (TH1F*)TSpectrum().Background(&hRefCo);
			hRefCo.Add(bgspec,-1.0);
			delete(bgspec);
			
			// fit for peaks
			defaultCanvas->cd();
			MultiGaus mg = multiPeakFitter(&hRefCo, foundPeaks);
			hRefCo.Draw();
			sprintf(tmp,"%s/Co60_%05i-%05i_%i.pdf",plotpath,startRun,endRun,n);
			defaultCanvas->Print(tmp);
			
			// process results
			for(unsigned int p=0; p<2; p++) {
				float_err t(segTimes.GetBinContent(n+1),segTimes.GetBinError(n+1)*sqrt(segTimes.GetBinEntries(n+1)));
				centers[p].push_back(tgePoint(t,mg.getPar(3*p+1)));
				widths[p].push_back(tgePoint(t,mg.getPar(3*p+2)));
				sprintf(tmp,"%i\t%i\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\n",p,(int)t.x,t.err,
						mg.getPar(3*p+1).x, mg.getPar(3*p+1).err,
						fabs(mg.getPar(3*p+2).x), mg.getPar(3*p+2).err);
				printf("%s",tmp);
				if( mg.getPar(3*p+1).x > plotmax || mg.getPar(3*p+1).x < plotmin ||  mg.getPar(3*p+1).err > 15 || 
				   fabs(mg.getPar(3*p+2).x) > 150 || fabs(mg.getPar(3*p+2).x) < 25 ||  mg.getPar(3*p+1).err > 15 )
					fprintf(errout,"Badpoint\t%i\t%i\t%i\n",n,p,(int)t.x);
				fprintf(peaksout,"%s",tmp);
			}
		}
		
		// make graphs for interpolation; draw and save interpolated points
		TGraphErrors gCenters[2] = { TGraphErrors(centers[0].size()), TGraphErrors(centers[1].size()) };
		TGraphErrors gWidths[2] = { TGraphErrors(widths[0].size()), TGraphErrors(widths[1].size()) };
		for(unsigned int p=0; p<2; p++) {
			for(unsigned int i=0; i<centers[p].size(); i++) {
				centers[p][i].setGraph(gCenters[p],i);
				widths[p][i].setGraph(gWidths[p],i);
			}
		}
		
		gCenters[0].SetMinimum(plotmin);
		gCenters[0].SetMaximum(plotmax);
		gCenters[0].Draw("AP");
		gCenters[1].Draw("P");
		sprintf(tmp,"%s/_Co60.pdf",plotpath);
		defaultCanvas->Print(tmp);
		
		for(unsigned int i=0; i<nruns; i++) {
			for(unsigned int p=0; p<2; p++) {
				unsigned int tstart = AnalysisDB::ADB().startTime(druns[i]);
				unsigned int tend = AnalysisDB::ADB().endTime(druns[i]);
				fprintf(peaksout,"%i\t%i\t0.0\t%.1f\t0.0\t%.1f\t0.0\n",p,tstart,gCenters[p].Eval(tstart),gWidths[p].Eval(tstart));
				fprintf(peaksout,"%i\t%i\t0.0\t%.1f\t0.0\t%.1f\t0.0\n",p,tend,gCenters[p].Eval(tend),gWidths[p].Eval(tend));
			}
		}
		
		fclose(peaksout);
		fclose(errout);
	}	
}


void trackPeds(RunNum startRun, RunNum endRun) {
	
	// select runs with known time stamps
	char tmp[1024];
	sprintf(tmp,"run_number >= %i AND run_number <= %i AND start_time IS NOT null AND end_time IS NOT null ORDER BY run_number ASC",startRun,endRun);
	std::vector<RunNum> druns = AnalysisDB::ADB().findRuns(tmp);
	unsigned int nruns = druns.size();
	printf("%i runs found.\n",nruns);
	if(!nruns)
		return;
	
	// set plot paths, output canvas
	char plotpath[1024];
	sprintf(plotpath,"../PostPlots/PeakShifts/R_%i_%i/",druns.front(),druns.back());
	makePath(plotpath);
	TCanvas defaultCanvas;
	defaultCanvas.SetFillColor(0);
	defaultCanvas.SetCanvasSize(300,300);
	
	// collect and process data
	for(Side s = EAST; s <= WEST; s = nextSide(s) ) {
		
		std::vector<std::string> wxnames;
		std::vector<std::string> wynames;
		std::vector<TGraphErrors*> wxCombo;
		std::vector<TGraphErrors*> wyCombo;
		for(unsigned int i=1; i<=16; i++) {
			wxnames.push_back(std::string("MWPC")+sideNames(s)+'x'+itos(i));
			wynames.push_back(std::string("MWPC")+sideNames(s)+'y'+itos(i));
		}
		for(unsigned int i=0; i<16; i++) {
			wxCombo.push_back(mergeloader(druns, wxnames[i], "pedestal",true).back());
			wyCombo.push_back(mergeloader(druns, wynames[i], "pedestal",true).back());
		}
		defaultCanvas.SetLogy(true);
		sprintf(tmp,"%s/WX_Ped_History_%c.pdf",plotpath,sideNames(s));
		drawTogether(wxCombo, 200, 2000, &defaultCanvas, tmp, "X Wires Pedestal History");
		sprintf(tmp,"%s/WY_Ped_History_%c.pdf",plotpath,sideNames(s));
		drawTogether(wyCombo, 200, 2000, &defaultCanvas, tmp, "Y Wires Pedestal History");
		defaultCanvas.SetLogy(false);
	}
}
