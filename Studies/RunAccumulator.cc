#include "RunAccumulator.hh"
#include "GraphUtils.hh"
#include "SMExcept.hh"

TRandom3 RunAccumulator::rnd_source;

fgbgPair::fgbgPair(const std::string& nm, const std::string& ttl, AFPState a, Side s):
baseName(nm), baseTitle(ttl), afp(a), mySide(s), doSubtraction(true), doTimeScale(true), isSubtracted(false) { }

void fgbgPair::bgSubtract(BlindTime tFG, BlindTime tBG) {
	assert(!isSubtracted); // don't BG subtract twice!
	double bgScale = (doTimeScale && tBG[BOTH])?tFG[BOTH]/tBG[BOTH]:1.0;
	if(doSubtraction)
		h[GV_OPEN]->Add(h[GV_CLOSED],-bgScale);
	else
		h[GV_CLOSED]->Scale(bgScale);
	isSubtracted = true;
}

void fgbgPair::operator+=(const fgbgPair& p) {
	for(GVState gv=GV_CLOSED; gv<=GV_OPEN; ++gv) {
		assert(h[gv] && p.h[gv]);
		h[gv]->Add(p.h[gv]);
	}
}

void fgbgPair::operator*=(double c) {
	for(GVState gv=GV_CLOSED; gv<=GV_OPEN; ++gv) {
		assert(h[gv]);
		h[gv]->Scale(c);
	}
}

void fgbgPair::setAxisTitle(AxisDirection d, const std::string& ttl) {
	if(d > Z_DIRECTION) return;
	for(GVState gv=GV_CLOSED; gv<=GV_OPEN; ++gv)
		(d==X_DIRECTION?h[gv]->GetXaxis():d==Y_DIRECTION?h[gv]->GetYaxis():h[gv]->GetZaxis())->SetTitle(ttl.c_str());
}

//------------------------------------------------------------------------

fgbgPair* RunAccumulator::registerFGBGPair(const std::string& hname, const std::string& title,
										  unsigned int nbins, float xmin, float xmax, AFPState a, Side s) {
	fgbgPair* p = new fgbgPair(hname,title,a,s);
	if(fgbgHists.find(p->getName())!=fgbgHists.end()) {
		SMExcept e("duplicateHistoName");
		e.insert("name",p->getName());
		throw e;
	}
	for(GVState gv=GV_CLOSED; gv<=GV_OPEN; ++gv) {
		p->h[gv] = registerSavedHist(p->getHistoName(gv),p->getHistoTitle(gv),nbins,xmin,xmax);
		if(!p->h[gv]->GetSumw2N())
			p->h[gv]->Sumw2();
	}
	fgbgHists.insert(std::make_pair(p->getName(),p));
	return p;
}

fgbgPair* RunAccumulator::registerFGBGPair(const TH1& hTemplate, AFPState a, Side s) {
	std::string hname = hTemplate.GetName();
	std::string htitle = hTemplate.GetTitle();
	fgbgPair* p = new fgbgPair(hname,htitle,a,s);
	if(fgbgHists.find(p->getName())!=fgbgHists.end()) {
		SMExcept e("duplicateHistoName");
		e.insert("name",p->getName());
		throw e;
	}
	for(int fg = 1; fg >=0; fg--) {
		p->h[fg] = registerSavedHist(p->getHistoName(fg),hTemplate);
		if(!p->h[fg]->GetSumw2N())
			p->h[fg]->Sumw2();
		p->h[fg]->SetTitle(p->getHistoTitle(fg).c_str());
	}
	fgbgHists.insert(std::make_pair(p->getName(),p));	
	return p;
}

fgbgPair* RunAccumulator::cloneFGBGPair(const fgbgPair& p, const std::string& newName, const std::string& newTitle) {
	fgbgPair* pnew = new fgbgPair(newName,newTitle,p.afp,p.mySide);
	for(GVState gv=GV_CLOSED; gv<=GV_OPEN; ++gv) {
		std::string qname = pnew->getHistoName(gv);
		pnew->h[gv] = (TH1*)addObject(p.h[gv]->Clone(qname.c_str()));
		pnew->h[gv]->SetTitle(pnew->getHistoTitle(gv).c_str());
	}
	return pnew;
}

RunAccumulator::RunAccumulator(OutputManager* pnt, const std::string& nm, const std::string& inflName):
SegmentSaver(pnt,nm,inflName), needsSubtraction(false), isSimulated(false) {
	
	// initialize blind time to 0
	zeroCounters();
	
	// load existing data (if any)
	if(fIn) {
		QFile qOld(inflname+".txt");
		// transfer octet data to new output file
		qOut.transfer(qOld, "Octet");
		// transfer run calibration data
		qOut.transfer(qOld, "runcal");
		// fetch total time
		std::vector<Stringmap> times = qOld.retrieve("totalTime");
		for(std::vector<Stringmap>::iterator it = times.begin(); it != times.end(); it++) {
			AFPState afp = strToAfp(it->getDefault("afp","Other"));
			unsigned int fg = (unsigned int)(it->getDefault("fg",3));
			assert(afp<=AFP_OTHER && fg <= GV_OPEN);
			totalTime[afp][fg] = BlindTime(*it);
		}
		// fetch total counts
		std::vector<Stringmap> counts = qOld.retrieve("totalCounts");
		for(std::vector<Stringmap>::iterator it = counts.begin(); it != counts.end(); it++) {
			AFPState afp = strToAfp(it->getDefault("afp","Other"));
			unsigned int fg = (unsigned int)(it->getDefault("fg",3));
			assert(fg <= GV_OPEN);
			totalCounts[afp][fg] = it->getDefault("counts",0);
		}
		// fetch run counts, run times
		runCounts += TagCounter<RunNum>(qOld.getFirst("runCounts"));
		runTimes += TagCounter<RunNum>(qOld.getFirst("runTimes"));
	}
}

RunAccumulator::~RunAccumulator() {
	for(std::map<std::string,fgbgPair*>::iterator it = fgbgHists.begin(); it != fgbgHists.end(); it++)
		delete it->second;
	for(std::map<std::string,AnalyzerPlugin*>::iterator it = myPlugins.begin(); it != myPlugins.end(); it++)
		delete it->second;
}

void RunAccumulator::setCurrentState(AFPState afp, GVState gv) {
	currentAFP = afp;
	currentGV = gv;
	for(std::map<std::string,AnalyzerPlugin*>::iterator it = myPlugins.begin(); it != myPlugins.end(); it++) {
		it->second->currentGV = currentGV;
		it->second->currentAFP = currentAFP;
	}
}

BlindTime RunAccumulator::getTotalTime(AFPState afp, bool fg) const {
	if(afp==AFP_OTHER)
		return totalTime[AFP_ON][fg]+totalTime[AFP_OFF][fg]+totalTime[AFP_OTHER][fg];
	return totalTime[afp][fg];
}


void RunAccumulator::addPlugin(AnalyzerPlugin* AP) {
	if(myPlugins.find(AP->name) != myPlugins.end()) {
		SMExcept e("DuplicatePluginName");
		e.insert("myName",name);
		e.insert("pluginName",AP->name);
		throw(e);
	}
	myPlugins.insert(std::make_pair(AP->name,AP));
}


AnalyzerPlugin* RunAccumulator::getPlugin(const std::string& nm) {
	std::map<std::string,AnalyzerPlugin*>::iterator it = myPlugins.find(nm);
	return it != myPlugins.end() ? it->second : NULL;

}

void RunAccumulator::fillCoreHists(ProcessedDataScanner& PDS, double weight) {
	for(std::map<std::string,AnalyzerPlugin*>::iterator it = myPlugins.begin(); it != myPlugins.end(); it++)
		it->second->fillCoreHists(PDS,weight);
}

void RunAccumulator::calculateResults() {
	for(std::map<std::string,AnalyzerPlugin*>::iterator it = myPlugins.begin(); it != myPlugins.end(); it++)
		it->second->calculateResults();
}

void RunAccumulator::uploadAnaResults() {
	for(std::map<std::string,AnalyzerPlugin*>::iterator it = myPlugins.begin(); it != myPlugins.end(); it++)
		it->second->uploadAnaResults();
}

void RunAccumulator::makePlots() {
	defaultCanvas->cd();
	for(std::map<std::string,AnalyzerPlugin*>::iterator it = myPlugins.begin(); it != myPlugins.end(); it++)
		it->second->makePlots();
}

void RunAccumulator::compareMCtoData(RunAccumulator& OAdata) {
	defaultCanvas->cd();
	for(std::map<std::string,AnalyzerPlugin*>::iterator it = myPlugins.begin(); it != myPlugins.end(); it++) {
		AnalyzerPlugin* AP = OAdata.getPlugin(it->second->name);
		if(AP) it->second->compareMCtoData(AP);
	}
}

void RunAccumulator::zeroCounters() {
	for(AFPState afp = AFP_OFF; afp <= AFP_OTHER; ++afp) {
		for(GVState gv=GV_CLOSED; gv<=GV_OPEN; ++gv) {
			totalTime[afp][gv] = BlindTime(0.0);
			totalCounts[afp][gv] = 0;
		}
	}
	runCounts = TagCounter<RunNum>();
	runTimes = TagCounter<RunNum>();
}

void RunAccumulator::addSegment(const SegmentSaver& S) {
	// histograms add
	SegmentSaver::addSegment(S);
	// recast
	const RunAccumulator& RA = (const RunAccumulator&)S;
	if(RA.isSimulated) isSimulated = true;
	// add times, counts
	for(AFPState afp = AFP_OFF; afp <= AFP_OTHER; ++afp) {
		for(GVState gv=GV_CLOSED; gv<=GV_OPEN; ++gv) {
			totalTime[afp][gv] += RA.totalTime[afp][gv];
			totalCounts[afp][gv] += RA.totalCounts[afp][gv];
		}
	}
	// add run counts, times
	runCounts += RA.runCounts;
	runTimes += RA.runTimes;
	// transfer run calibration data
	qOut.transfer(S.qOut,"runcal");
}

void RunAccumulator::scaleData(double s) {
	SegmentSaver::scaleData(s);
	for(AFPState afp = AFP_OFF; afp <= AFP_OTHER; ++afp)
		for(GVState fg = GV_CLOSED; fg <= GV_OPEN; ++fg)
			totalCounts[afp][fg] *= s;
	runCounts.scale(s);
}

bool RunAccumulator::hasFGBGPair(const std::string& qname) const {
	return fgbgHists.count(qname);
}

fgbgPair& RunAccumulator::getFGBGPair(const std::string& qname) {
	std::map<std::string,fgbgPair*>::iterator it = fgbgHists.find(qname);
	assert(it != fgbgHists.end());
	return *(it->second);
}

const fgbgPair& RunAccumulator::getFGBGPair(const std::string& qname) const {
	std::map<std::string,fgbgPair*>::const_iterator it = fgbgHists.find(qname);
	assert(it != fgbgHists.end());
	return *(it->second);
}

RunAccumulator* RunAccumulator::getErrorEstimator() {
	std::string epath = estimatorHistoLocation();
	if(!epath.size()) return NULL;
	static std::map<std::string,RunAccumulator*> estimators;
	std::map<std::string,RunAccumulator*>::const_iterator it = estimators.find(epath);
	if(it != estimators.end()) return it->second;
	RunAccumulator* OA = NULL;
	if(fileExists(epath+".root") && fileExists(epath+".txt"))
		OA = (RunAccumulator*)makeAnalyzer("MasterRates",epath);
	else
		printf("*** Unable to locate master histograms at '%s'\n",epath.c_str());
	estimators.insert(std::make_pair(epath,OA));
	return getErrorEstimator();
}

/// estimate errorbars for low-rate histogram from high-rate total events; assume both histograms are unscaled counts
void errorbarsFromMasterHisto(TH1* lowrate, const TH1* master) {
	assert(lowrate && master);
	if(!master->Integral()) return;
	float errscaling = sqrt(lowrate->Integral()/master->Integral());
	assert(errscaling==errscaling);
	for(unsigned int i=0; i<totalBins(lowrate); i++)
		if(lowrate->GetBinContent(i)<25)
			lowrate->SetBinError(i,errscaling*master->GetBinError(i));
}

void RunAccumulator::bgSubtractAll() {
	for(std::map<std::string,fgbgPair*>::iterator it = fgbgHists.begin(); it != fgbgHists.end(); it++) {
		if(getErrorEstimator())
			errorbarsFromMasterHisto(it->second->h[GV_CLOSED],getErrorEstimator()->getFGBGPair(it->second->getName()).h[GV_CLOSED]);
			it->second->bgSubtract(getTotalTime(it->second->afp,GV_OPEN),getTotalTime(it->second->afp,GV_CLOSED));
	}
	needsSubtraction = false;
}

void RunAccumulator::simBgFlucts(const RunAccumulator& RefOA, double simfactor, bool addFluctCounts) {
	printf("Adding background fluctuations to simulation...\n");
	for(std::map<std::string,fgbgPair*>::iterator it = fgbgHists.begin(); it != fgbgHists.end(); it++) {
		if(!RefOA.hasFGBGPair(it->first)) continue;
		const fgbgPair& qhRef = RefOA.getFGBGPair(it->first);
		double bgRatio = RefOA.getTotalTime(it->second->afp,true)[it->second->mySide]/RefOA.getTotalTime(it->second->afp,false)[it->second->mySide];
		for(unsigned int i=0; i<totalBins(it->second->h[0]); i++) {
			double rootn = qhRef.h[0]->GetBinError(i)*sqrt(simfactor);		// root(bg counts) from ref histogram errorbars
			double n = rootn*rootn;											// background counts from reference histogram
			double bgObsCounts = rnd_source.PoissonD(n);					// simulated background counts
			double fgBgCounts = rnd_source.PoissonD(n*bgRatio);				// simulated foreground counts due to background
			if(!addFluctCounts) bgObsCounts = fgBgCounts = 0.;
			it->second->h[0]->AddBinContent(i,bgObsCounts);					// fill fake background counts
			it->second->h[0]->SetBinError(i,rootn);							// set root-n statitics for background
			it->second->h[1]->AddBinContent(i,fgBgCounts);					// add simulated background to foreground histogram				
		}
		printf("\t%i counts for %s [%i bins]\n",int(it->second->h[0]->Integral()),it->second->getName().c_str(),totalBins(it->second->h[0]));
		it->second->h[1]->Add(it->second->h[0],-bgRatio);						// subtract back off simulated background
		it->second->isSubtracted = true;
	}
}

void RunAccumulator::makeRatesSummary() {
	for(std::map<std::string,fgbgPair*>::const_iterator it = fgbgHists.begin(); it != fgbgHists.end(); it++) {
		for(GVState gv=GV_CLOSED; gv<=GV_OPEN; ++gv) {
			Stringmap rt;
			rt.insert("side",ctos(sideNames(it->second->mySide)));
			rt.insert("afp",afpWords(it->second->afp));
			rt.insert("name",it->second->getName());
			rt.insert("fg",itos(gv));
			Double_t d_counts;
			double counts = it->second->h[gv]->IntegralAndError(-1, -1, d_counts);
			double rtime = getTotalTime(it->second->afp,gv)[it->second->mySide];
			rt.insert("counts",counts);
			rt.insert("d_counts",d_counts);
			rt.insert("rate",counts?counts/rtime:0);
			rt.insert("d_rate",d_counts?d_counts/rtime:0);
			qOut.insert("rate",rt);
		}
	}
}

void RunAccumulator::write(std::string outName) {
	// record total times, counts
	for(AFPState afp = AFP_OFF; afp <= AFP_OTHER; ++afp) {
		for(GVState gv=GV_CLOSED; gv<=GV_OPEN; ++gv) {
			Stringmap tm = totalTime[afp][gv].toStringmap();
			tm.insert("afp",afpWords(afp));
			tm.insert("fg",itos(gv));
			qOut.insert("totalTime",tm);
			Stringmap ct;
			ct.insert("afp",afpWords(afp));
			ct.insert("fg",itos(gv));
			ct.insert("counts",totalCounts[afp][gv]);
			qOut.insert("totalCounts",ct);
		}
	}
	// record run counts, times
	qOut.insert("runCounts",runCounts.toStringmap());
	qOut.insert("runTimes",runTimes.toStringmap());
	// base class write
	SegmentSaver::write(outName);
}

void RunAccumulator::loadProcessedData(AFPState afp, GVState gv, ProcessedDataScanner& PDS) {
	printf("Loading AFP=%i, fg=%i processed data...\n",afp,gv);
	assert(afp <= AFP_OTHER);
	assert(gv==GV_CLOSED || gv==GV_OPEN);
	setCurrentState(afp,gv);
	if(!PDS.getnFiles())
		return;
	PDS.startScan();
	unsigned int nScanned = 0;
	while(PDS.nextPoint()) {
		nScanned++;
		if(PDS.withCals)
			PDS.recalibrateEnergy();
		if(PDS.fPID==PID_BETA && PDS.fType==TYPE_0_EVENT) {
			runCounts.add(PDS.getRun(),1.0);
			totalCounts[afp][gv]++;
		}
		fillCoreHists(PDS,PDS.physicsWeight);
	}
	printf("\tFG=%i: scanned %i points\n",gv,nScanned);
	if(gv==GV_CLOSED)
		needsSubtraction = true;
	runTimes += PDS.runTimes;
	totalTime[afp][gv] += PDS.totalTime;
	PDS.writeCalInfo(qOut,"runcal");
}

void RunAccumulator::loadSimData(Sim2PMT& simData, unsigned int nToSim, bool countAll) {
	isSimulated = true;
	setCurrentState(simData.getAFP(),GV_OPEN);
	printf("Loading %i events of simulated data (AFP=%i)...\n",nToSim,currentAFP);
	simData.resetSimCounters();
	simData.startScan(nToSim);
	while(!nToSim || (countAll?simData.nSimmed:simData.nCounted)<=nToSim) {
		bool np = simData.nextPoint();
		loadSimPoint(simData);
		if(nToSim && !(int(simData.nSimmed)%(nToSim/20))) {
			if(nToSim>1e6) {
				printf("* %s\n",simData.evtInfo().toString().c_str());
			} else {
				printf("*");
				fflush(stdout);
			}
		}
		if(!nToSim && !np) break;
	}
	printf("\n--Scan complete.--\n");
}

void RunAccumulator::loadSimPoint(Sim2PMT& simData) {
	fillCoreHists(simData,simData.physicsWeight);
	if(double evtc = simData.simEvtCounts()) {
		runCounts.add(simData.getRun(),evtc);
		totalCounts[currentAFP][1] += evtc;
	}	
}

void RunAccumulator::simForRun(Sim2PMT& simData, RunNum rn, unsigned int nToSim, bool countAll) {
	RunInfo RI = CalDBSQL::getCDB()->getRunInfo(rn);
	if(RI.gvState != GV_OPEN) return;
	assert(RI.afpState <= AFP_OTHER);
	
	PMTCalibrator PCal(rn);
	simData.setCalibrator(PCal);
	simData.setAFP(RI.afpState);
	loadSimData(simData,nToSim,countAll);
	
	double rntime = CalDBSQL::getCDB()->fiducialTime(rn)[BOTH];
	runTimes.add(rn,rntime);
	totalTime[RI.afpState][GV_OPEN] += rntime;
}

unsigned int RunAccumulator::simMultiRuns(Sim2PMT& simData, const TagCounter<RunNum>& runReqs, unsigned int nCounts) {
	unsigned int nSimmed = 0;
	// if nCounts==0, simulate all at requested levels
	if(!nCounts) {
		for(std::map<RunNum,double>::const_iterator it = runReqs.counts.begin(); it != runReqs.counts.end(); it++) {
			simForRun(simData, it->first, it->second, false);
			nSimmed += simData.nSimmed;
		}
	} else {
		double nRequested = runReqs.total();
		assert(nRequested);
		double nGranted = 0;
		printf("Dividing %i simulation events between %i runs requesting %i events...\n",nCounts,runReqs.nTags(),(int)nRequested);
		for(std::map<RunNum,double>::const_iterator it = runReqs.counts.begin(); it != runReqs.counts.end(); it++) {
			// calculate alloted number of events for this run
			nGranted += it->second;
			int nToSim = int((nGranted/nRequested)*nCounts)-nSimmed;
			// simulate alloted requests and re-scale to requested counts
			RunAccumulator* subRA = (RunAccumulator*)makeAnalyzer("nameUnused","");
			subRA->simForRun(simData, it->first, nToSim, true);
			nSimmed += simData.nSimmed;
			printf("From %i input points, simulated %i/%i requested events for Run %i\n",simData.nSimmed,(int)simData.nCounted,(int)it->second,it->first);
			subRA->scaleData(it->second/simData.nCounted);
			addSegment(*subRA);
			delete(subRA);
		}
	}
	return nSimmed;
}
