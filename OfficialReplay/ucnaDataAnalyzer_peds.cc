#include "ucnaDataAnalyzer11b.hh"
#include "GraphicsUtils.hh"
#include <TProfile.h>
#include <TGraph.h>
#include <utility>

void ucnaDataAnalyzer11b::pedestalPrePass() {
	
	// get total run time; force pre-pass if missing
	wallTime = PCal.CDB->totalTime(rn);
	unsigned int needsPeds = !wallTime;
	if(needsPeds)
		printf("Total run time not found in DB; force pre-pass\n");
	else {
		// check if pedestal pre-pass needed
		printf("Checking for pedestals data...\n");
		for(Side s = EAST; s <= WEST; s = nextSide(s)) {
			for(unsigned int t=0; t<nBetaTubes; t++)
				needsPeds += !PCal.checkPedestals(PCal.sensorNames[s][t]);
			for(unsigned int p = X_DIRECTION; p <= Y_DIRECTION; p++)
				for(unsigned int c=0; c<cathNames[s][p].size(); c++)
					needsPeds += !PCal.checkPedestals(cathNames[s][p][c]);
			needsPeds += !PCal.checkPedestals(sideSubst("MWPC%cAnode",s));
		}
	}
	if(!needsPeds) return;
	
	printf("Pre-pass for pedestals and run time...\n");
	
	// collect pedestal data points
	std::vector< std::pair<float,float> > pmtPeds[2][nBetaTubes];
	std::vector< std::pair<float,float> > anodePeds[2];
	std::vector< std::pair<float,float> > cathPeds[2][2][kMWPCWires];
	startScan();
	while (nextPoint()) {
		calibrateTimes();
		for(Side s = EAST; s <= WEST; s = nextSide(s)) {
			if( !trig2of4(s) && qadcSum(otherSide(s))>2000 && !isLED() )
				for(unsigned int t=0; t<nBetaTubes; t++)
					pmtPeds[s][t].push_back(std::make_pair(fTimeScaler.t[BOTH],sevt[s].adc[t]));
			if(isLED() || (isPulserTrigger() && !nFiring(s)) || isUCNMon()) {
				for(unsigned int p = X_DIRECTION; p <= Y_DIRECTION; p++)
					for(unsigned int c=0; c<kMWPCWires; c++)
						cathPeds[s][p][c].push_back(std::make_pair(fTimeScaler.t[BOTH],fMWPC_caths[s][p][c]));
				anodePeds[s].push_back(std::make_pair(fTimeScaler.t[BOTH],fMWPC_anode[s].val));
			}
		}
	}
	
	// fit pedestals, save results
	for(Side s = EAST; s <= WEST; s = nextSide(s)) {
		for(unsigned int t=0; t<nBetaTubes; t++)
			monitorPedestal(pmtPeds[s][t],PCal.sensorNames[s][t],50);
		for(unsigned int p = X_DIRECTION; p <= Y_DIRECTION; p++)
			for(unsigned int c=0; c<cathNames[s][p].size(); c++)
				monitorPedestal(cathPeds[s][p][c],cathNames[s][p][c],150);
		monitorPedestal(anodePeds[s],sideSubst("MWPC%cAnode",s),100);
	}
	
	// re-set for next scan
	wallTime = totalTime.t[BOTH];
}

void ucnaDataAnalyzer11b::monitorPedestal(std::vector< std::pair<float,float> > dpts, const std::string& mon_name, double graphWidth) {
	
	printf("Monitoring data '%s'\n",mon_name.c_str());
	// collect data
	unsigned int npts = dpts.size();
	float t0 = dpts[0].first;
	float t1 = dpts.back().first;
	printf("\tfound %i points over %.2f minutes.\n",npts,(t1-t0)/60.0);
	// determine data division
	const float tmin = 60.0;
	const unsigned int cmin = 3000;
	unsigned int ndivs = (unsigned int)((t1-t0)/tmin);
	if(npts/cmin < ndivs)
		ndivs = npts/cmin;
	if(!ndivs)
		ndivs = 1;
	printf("\tdividing into %i intervals.\n",ndivs);
	
	// estimate data range using TProfile
	TProfile* p = new TProfile("pmon","PeakMonitor",ndivs,0,dpts.size()+1);
	TProfile* pTime = new TProfile("pmonT","PeakMonitor_Time",ndivs,0,dpts.size()+1);
	for(unsigned int i=0; i<dpts.size(); i++) {
		pTime->Fill(i,dpts[i].first);
		p->Fill(i,dpts[i].second);
	}
	printf("\tEstimated peak in each interval.\n");
	
	// extract pedestal in each interval
	QFile pedOut(dataPath+"/Monitors/Run_"+itos(rn)+"/"+mon_name+"_ped.txt",false);
	std::vector< std::pair<float,float> >::iterator it = dpts.begin();
	defaultCanvas->cd();
	TGraph* tg = new TGraph(ndivs>1?ndivs:2);
	std::vector<double> times;
	std::vector<double> dtimes;
	std::vector<double> centers;
	std::vector<double> dcenters;
	std::vector<double> sigmas;
	std::vector<double> dwidths;
	std::vector<TH1*> hToPlot;
	printf("\tMaking histograms for each division...\n");
	for(unsigned int i=0; i<ndivs; i++) {
		// book histogram based on mean/sigma estimate
		float c = p->GetBinContent(i+1);
		int x0 = int(c-graphWidth);
		int x1 = int(c+graphWidth);
		TH1F* hdiv = registeredTH1F(mon_name+"_Mon_Div_"+itos(i),mon_name+" Pedestals",x1-x0,x0,x1);
		while(it != dpts.end() && int(it-dpts.begin()) <= float((i+1)*dpts.size())/float(ndivs)) {
			hdiv->Fill(it->second);
			it++;
		}
		hToPlot.push_back(hdiv);
		centers.push_back(hdiv->GetMean());
		sigmas.push_back(hdiv->GetRMS());
		dcenters.push_back(sigmas.back()/sqrt(hdiv->GetEntries()));
		times.push_back(pTime->GetBinContent(i+1));
		dtimes.push_back(pTime->GetBinError(i+1));
		Stringmap m;
		m.insert("n",i);
		m.insert("t",times.back());
		m.insert("dt",dtimes.back());
		m.insert("ped",centers.back());
		m.insert("pedw",sigmas.back());
		m.insert("dped",dcenters.back());
		pedOut.insert("pedestal",m);
		tg->SetPoint(i,times.back(),centers.back());
	}
	if(ndivs==1) {
		printf("Notice: only 1 graph point found; extending to 2.");
		tg->SetPoint(1,p->GetBinCenter(1)+10.0,centers[0]);
	}
	printf("\tMaking plots...\n");
	drawSimulHistos(hToPlot);
	for(unsigned int i=0; i<ndivs; i++) {
		drawVLine(centers[i],defaultCanvas,2);
		drawVLine(centers[i]+sigmas[i],defaultCanvas,4);
		drawVLine(centers[i]-sigmas[i],defaultCanvas,4);		
	}
	printCanvas("Pedestals/"+mon_name);
	
	// load pedestals graph into memory
	PCal.insertPedestal(mon_name,tg);
	pedOut.commit();
	
	// optionally add to Calibrations DB
	if(CDBout) {
		printf("Uploading pedestal '%s'...\n",mon_name.c_str());
		unsigned int cgid = CDBout->uploadGraph(itos(rn)+" "+mon_name+" Pedestal Centers",times,centers,dtimes,dcenters);
		unsigned int wgid = CDBout->uploadGraph(itos(rn)+" "+mon_name+" Pedestal Widths",times,sigmas,dtimes);
		CDBout->deleteRunMonitor(rn,mon_name,"pedestal");
		CDBout->addRunMonitor(rn,mon_name,"pedestal",cgid,wgid);
	}
	
	delete(p);
	delete(pTime);
}
