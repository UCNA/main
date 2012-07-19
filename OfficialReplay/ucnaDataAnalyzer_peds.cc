#include "ucnaDataAnalyzer11b.hh"
#include "GraphicsUtils.hh"
#include <TProfile.h>
#include <TGraph.h>
#include <utility>

void ucnaDataAnalyzer11b::pedestalPrePass() {
	
	// get total run time; force pre-pass if missing
	wallTime = PCal.CDB->totalTime(rn);
	needsPeds |= !wallTime;
	if(needsPeds)
		printf("Total run time not found in DB; force pre-pass\n");
	else {
		// check if pedestal pre-pass needed
		printf("Checking for pedestals data...\n");
		for(Side s = EAST; s <= WEST; ++s) {
			for(unsigned int t=0; t<nBetaTubes; t++)
				needsPeds |= !PCal.checkPedestals(PCal.sensorNames[s][t]);
			for(unsigned int p = X_DIRECTION; p <= Y_DIRECTION; p++)
				for(unsigned int c=0; c<cathNames[s][p].size(); c++)
					needsPeds |= !PCal.checkPedestals(cathNames[s][p][c]);
			needsPeds |= !PCal.checkPedestals(sideSubst("MWPC%cAnode",s));
		}
	}
	if(!needsPeds) return;
	
	printf("Pre-pass for pedestals and run time...\n");
	
	// collect pedestal data points
	std::vector<float> pmtPeds[2][nBetaTubes];
	std::vector<float> pmtTimes[2];
	std::vector<float> anodePeds[2];
	std::vector<float> cathPeds[2][2][kMaxCathodes];
	std::vector<float> mwpcTimes[2];
	startScan();
	while (nextPoint()) {
		convertReadin();
		calibrateTimes();
		for(Side s = EAST; s <= WEST; ++s) {
			if(isUCNMon() || iSis00==(s==EAST?2:1)) {
				pmtTimes[s].push_back(fTimeScaler[BOTH]);
				for(unsigned int t=0; t<nBetaTubes; t++)
					pmtPeds[s][t].push_back(sevt[s].adc[t]);
			}
			if(isLED() || (isPulserTrigger() && !nFiring(s)) || isUCNMon()) {
				mwpcTimes[s].push_back(fTimeScaler[BOTH]);
				for(AxisDirection p = X_DIRECTION; p <= Y_DIRECTION; ++p)
					for(unsigned int c=0; c<kMaxCathodes; c++)
						cathPeds[s][p][c].push_back(fMWPC_caths[s][p][c]);
				anodePeds[s].push_back(fMWPC_anode[s].val);
			}
		}
	}
	
	// fit pedestals, save results
	for(Side s = EAST; s <= WEST; ++s) {
		for(unsigned int t=0; t<nBetaTubes; t++)
			monitorPedestal(pmtPeds[s][t],pmtTimes[s],PCal.sensorNames[s][t],50);
		for(AxisDirection p = X_DIRECTION; p <= Y_DIRECTION; ++p)
			for(unsigned int c=0; c<cathNames[s][p].size(); c++)
				monitorPedestal(cathPeds[s][p][c],mwpcTimes[s],cathNames[s][p][c],150);
		monitorPedestal(anodePeds[s],mwpcTimes[s],sideSubst("MWPC%cAnode",s),100);
	}
	
	// re-set for next scan
	wallTime = totalTime[BOTH];
}

void ucnaDataAnalyzer11b::monitorPedestal(const std::vector<float>& vdata, const std::vector<float>& vtime,
										  const std::string& mon_name, double graphWidth, float tmin, unsigned int cmin, bool isPed) {
	
	printf("Monitoring data '%s'\n",mon_name.c_str());
	// collect data
	unsigned int npts = vdata.size();
	assert(vtime.size()==npts); assert(npts);
	float t0 = vtime[0];
	float t1 = vtime.back();
	printf("\tfound %i points over %.2f minutes.\n",npts,(t1-t0)/60.0);
	// determine data division
	unsigned int ndivs = (unsigned int)((t1-t0)/tmin);
	if(npts/cmin < ndivs)
		ndivs = npts/cmin;
	if(!ndivs)
		ndivs = 1;
	printf("\tdividing into %i intervals.\n",ndivs);
	
	// estimate data range using TProfile
	TProfile* p = new TProfile("pmon","PeakMonitor",ndivs,0,npts+1);
	p->SetErrorOption("s");
	TProfile* pTime = new TProfile("pmonT","PeakMonitor_Time",ndivs,0,npts+1);
	for(unsigned int i=0; i<npts; i++) {
		pTime->Fill(i,vtime[i]);
		p->Fill(i,vdata[i]);
	}
	printf("\tEstimated peak in each interval.\n");
	
	// extract pedestal in each interval
	defaultCanvas->cd();
	TGraph* tg = new TGraph(ndivs>1?ndivs:2);
	TGraph* tgw = new TGraph(ndivs>1?ndivs:2);
	std::vector<double> times;
	std::vector<double> dtimes;
	std::vector<double> centers;
	std::vector<double> dcenters;
	std::vector<double> sigmas;
	std::vector<double> dwidths;
	std::vector<TH1*> hToPlot;
	unsigned int n=0;
	bool autoWidth = !graphWidth;
	printf("\tMaking histograms for each division...\n");
	for(unsigned int i=0; i<ndivs; i++) {
		// book histogram based on mean/sigma estimate
		float c = p->GetBinContent(i+1);
		if(autoWidth)
			graphWidth = 3*p->GetBinError(i+1);
		int x0 = int(c-graphWidth);
		int x1 = int(c+graphWidth);
		TH1F* hdiv = registeredTH1F(mon_name+"_Mon_Div_"+itos(i),mon_name+" Pedestals",x1-x0,x0,x1);
		while(n<npts && n <= float((i+1)*npts)/float(ndivs)) {
			hdiv->Fill(vdata[n]);
			n++;
		}
		hToPlot.push_back(hdiv);
		centers.push_back(hdiv->GetMean());
		sigmas.push_back(hdiv->GetRMS());
		dcenters.push_back(sigmas.back()/sqrt(hdiv->GetEntries()));
		times.push_back(pTime->GetBinContent(i+1));
		dtimes.push_back(pTime->GetBinError(i+1));
		tg->SetPoint(i,times.back(),centers.back());
		tgw->SetPoint(i,times.back(),sigmas.back());
	}
	if(ndivs==1) {
		printf("Notice: only 1 graph point found; extending to 2.");
		tg->SetPoint(1,p->GetBinCenter(1)+10.0,centers[0]);
		tgw->SetPoint(1,p->GetBinCenter(1)+10.0,sigmas[0]);
	}
	printf("\tMaking plots...\n");
	if(isPed) {
		drawSimulHistos(hToPlot);
		for(unsigned int i=0; i<ndivs; i++) {
			drawVLine(centers[i],defaultCanvas,2);
			drawVLine(centers[i]+sigmas[i],defaultCanvas,4);
			drawVLine(centers[i]-sigmas[i],defaultCanvas,4);		
		}
		printCanvas("Pedestals/"+mon_name);
	} else {
		tg->Draw("AP");
		tg->SetTitle(mon_name.c_str());
		tg->GetXaxis()->SetTitle("Time [s]");
		tg->GetYaxis()->SetTitle("Value");
		tg->Draw("ALP");
		tgw->SetLineColor(2);
		tgw->Draw("LP");
		printCanvas("Pedestals/"+mon_name);
	}
	
	// load pedestals graph into memory
	if(isPed) PCal.insertPedestal(mon_name,tg);
	
	// optionally add to Calibrations DB
	if(CDBout) {
		printf("Uploading pedestal '%s'...\n",mon_name.c_str());
		unsigned int cgid = CDBout->uploadGraph(itos(rn)+" "+mon_name+(isPed?" Pedestal":" Peak")+" Centers",times,centers,dtimes,dcenters);
		unsigned int wgid = CDBout->uploadGraph(itos(rn)+" "+mon_name+(isPed?" Pedestal":" Peak")+" Widths",times,sigmas,dtimes);
		CDBout->deleteRunMonitor(rn,mon_name,isPed?"pedestal":"GMS_peak");
		CDBout->addRunMonitor(rn,mon_name,isPed?"pedestal":"GMS_peak",cgid,wgid);
	}
	
	delete(p);
	delete(tgw);
	delete(pTime);
}
