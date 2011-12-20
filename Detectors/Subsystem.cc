#include "Subsystem.hh"
#include "RunManager.hh"
#include "PathUtils.hh"
#include "GraphicsUtils.hh"

Subsystem::Subsystem(RunManager* T, std::string nm, Side s):
OutputManager(nm, T), mySide(s), PC(T->RI.runNum,T->CDB), thisRun(T), nEvents(T->nEvents) {
	assert(false); // TODO removed registerToParent for objects; should eventually fully deprecate this code
	triggers = new bool[nEvents];
	T->attachSubsystem(this);
}

float* Subsystem::getData(std::string dataName, int offset) {
	std::map<std::string,unsigned int>::const_iterator cit = nameMap.find(dataName);
	assert(cit != nameMap.end());
	return edata[cit->second+offset]; 
}
const float* Subsystem::getData(std::string dataName, int offset) const { 
	std::map<std::string,unsigned int>::const_iterator cit = nameMap.find(dataName);
	assert(cit != nameMap.end());
	return edata[cit->second+offset]; 
}

void Subsystem::loadData(std::string branchName, std::string sensorName, std::string dataName) {
	assert(thisRun->checkBranch(branchName));
	std::map<std::string,unsigned int>::const_iterator cit = nameMap.find(dataName);
	if(dataName != "" && cit == nameMap.end())
		nameMap[dataName] = edata.size();
	edata.push_back(thisRun->readBranch(branchName));
	sensorNames.push_back(sensorName);
}

#define FORCE_REMAKE_PED true
bool Subsystem::verifyPedestal(int n) {
	bool needsPeds = !PC.checkPedestals(sensorNames[n]);
	if(FORCE_REMAKE_PED || needsPeds) {
		if(needsPeds)
			warn(BENIGN_WARNING,std::string("Missing_")+sensorNames[n]+"_Pedestals");
		else
			printf("Pedestals found, but forcing remake anyway...\n");
		return false;
	}
	return true;
}

void Subsystem::bitsHisto(UInt_t n, bool* selector) const {
	TH1F* foo0 = new TH1F("zeros","Zeros (black), Ones (red), Not 11111...1 (blue)",13,-0.5,12.5);
	TH1F* foo1 = new TH1F("ones","Bits Distribution",13,-0.5,12.5);
	TH1F* foo11 = new TH1F("ones","Bits Distribution",13,-0.5,12.5);
	foo1->SetLineColor(2);
	foo11->SetLineColor(4);
	
	for(UInt_t e = 0; e<nEvents; e++) {
		int z = int(edata[n][e]);
		if( selector && !selector[e] )
			continue;
		for(UInt_t b = 0; b<13; b++) {
			if((z>>b)%2) {
				foo1->Fill(b);
				if(z < (1<<12)-1)
					foo11->Fill(b);
			}
			else
				foo0->Fill(b);
		}
	}
	
	foo0->SetMinimum(0);
	foo0->Draw();
	foo1->Draw("Same");
	foo11->Draw("Same");
	printCanvas(std::string("Bits_")+itos(n));
	delete(foo0);
	delete(foo1);
}

void Subsystem::monitorPeak(std::vector< std::pair<float,float> > (*selector)(Subsystem*,void*),
								 float tmin, unsigned int cmin, void* selector_params, std::string mon_name, bool savePed) {
	printf("Monitoring data '%s'\n",mon_name.c_str());
	// collect data
	std::vector< std::pair<float,float> > dpts = selector(this,selector_params);
	unsigned int npts = dpts.size();
	float t0 = dpts[0].first;
	float t1 = dpts.back().first;
	printf("\tfound %i points over %.2f minutes.\n",npts,(t1-t0)/60.0);
	// determine data division scheme
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
	
	// make plots for each division
	std::string peakpath = thisRun->basePath+"/RunData/PeakMonitors/Run_"+itos(thisRun->RI.runNum)+"/";
	makePath(peakpath);
	FILE* f = fopen( (peakpath + mon_name + "_ped.txt").c_str(), "w");
	std::vector< std::pair<float,float> >::iterator it = dpts.begin();
	defaultCanvas->cd();
	TGraph* tg = new TGraph(ndivs>1?ndivs:2);
	std::vector<float> centers;
	std::vector<float> sigmas;
	std::vector<TH1*> hToPlot;
	printf("\tMaking histograms for each division...\n");
	for(unsigned int i=0; i<ndivs; i++) {
		// book histogram based on mean/sigma estimate
		float c = p->GetBinContent(i+1);
		float dc = p->GetBinError(i+1);
		float sigma = dc*sqrt(p->GetBinEntries(i+1));
		int x0 = int(c-4*sigma)-1;
		int x1 = int(c+4*sigma)+1;
		TH1F* hdiv = registeredTH1F(mon_name+"_Mon_Div_"+itos(i),"Peak Monitor",x1-x0,x0,x1);
		while(it != dpts.end() && int(it-dpts.begin()) <= float((i+1)*dpts.size())/float(ndivs)) {
			hdiv->Fill(it->second);
			it++;
		}
		hToPlot.push_back(hdiv);
		centers.push_back(hdiv->GetMean());
		sigmas.push_back(hdiv->GetRMS());
		float tcenter = pTime->GetBinContent(i+1);
		fprintf(f,"%.2f\t%.2f\t%.4f\t%.2f\n",tcenter,centers.back(),sigmas.back()/sqrt(hdiv->GetEntries()),sigmas.back());
		printf("\t\t%i\t%.2f\t%.2f\t%.4f\t%.2f\n",i,tcenter,centers.back(),sigmas.back()/sqrt(hdiv->GetEntries()),sigmas.back());
		tg->SetPoint(i,t0,centers.back());
	}
	if(ndivs==1) {
		printf("Notice: only 1 graph point found; extending to 2.");
		tg->SetPoint(1,p->GetBinCenter(1)+10.0,centers[0]);
	}
	printf("\tMaking plots...\n");
	defaultCanvas->cd();
	drawSimulHistos(hToPlot);
	for(unsigned int i=0; i<ndivs; i++) {
		drawVLine(centers[i],defaultCanvas,2);
		drawVLine(centers[i]+sigmas[i],defaultCanvas,4);
		drawVLine(centers[i]-sigmas[i],defaultCanvas,4);		
	}
	printCanvas(mon_name+"_Tracker");
	fclose(f);
	
	
	// load pedestals graph into memory
	if(savePed) {
		PC.insertPedestal(mon_name,tg);
	} else {
		delete(tg);
	}
	
	delete(p);
	delete(pTime);
}
