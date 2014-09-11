#include "RunAccumulator.hh"
#include "GraphUtils.hh"
#include "SMExcept.hh"
#include "PostOfficialAnalyzer.hh"
#include <time.h>

TRandom3 RunAccumulator::rnd_source;

fgbgPair::fgbgPair(const std::string& nm, const std::string& ttl, AFPState a, Side s):
baseName(nm), baseTitle(ttl), afp(a), mySide(s), doSubtraction(true), doTimeScale(true), isSubtracted(false) { }

void fgbgPair::bgSubtract(BlindTime tFG, BlindTime tBG) {
	smassert(!isSubtracted); // don't BG subtract twice!
	double bgScale = (doTimeScale && tBG[BOTH])?tFG[BOTH]/tBG[BOTH]:1.0;
	if(doSubtraction)
		h[GV_OPEN]->Add(h[GV_CLOSED],-bgScale);
	else
		h[GV_CLOSED]->Scale(bgScale);
	isSubtracted = true;
}

void fgbgPair::operator+=(const fgbgPair& p) {
	for(GVState gv=GV_CLOSED; gv<=GV_OPEN; ++gv) {
		smassert(h[gv] && p.h[gv]);
		h[gv]->Add(p.h[gv]);
	}
}

void fgbgPair::operator*=(double c) {
	for(GVState gv=GV_CLOSED; gv<=GV_OPEN; ++gv) {
		smassert(h[gv]);
		h[gv]->Scale(c);
	}
}

void fgbgPair::setAxisTitle(AxisDirection d, const std::string& ttl) {
	if(d > Z_DIRECTION) return;
	for(GVState gv=GV_CLOSED; gv<=GV_OPEN; ++gv)
		(d==X_DIRECTION?h[gv]->GetXaxis():d==Y_DIRECTION?h[gv]->GetYaxis():h[gv]->GetZaxis())->SetTitle(ttl.c_str());
}

//------------------------------------------------------------------------

std::string RunAccumulator::processedLocation = "";
std::string RunAccumulator::AnaDB_xtag = "";

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
SegmentSaver(pnt,nm,inflName), needsSubtraction(false), isSimulated(false), grouping(GROUP_RANGE), simPerfectAsym(false) {
	
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
			smassert(afp<=AFP_OTHER && fg <= GV_OPEN);
			totalTime[afp][fg] = BlindTime(*it);
		}
		// fetch total counts
		std::vector<Stringmap> counts = qOld.retrieve("totalCounts");
		for(std::vector<Stringmap>::iterator it = counts.begin(); it != counts.end(); it++) {
			AFPState afp = strToAfp(it->getDefault("afp","Other"));
			unsigned int fg = (unsigned int)(it->getDefault("fg",3));
			smassert(fg <= GV_OPEN);
			totalCounts[afp][fg] = it->getDefault("counts",0);
		}
		// fetch run counts, run times
		runCounts += TagCounter<RunNum>(qOld.getFirst("runCounts"));
		runTimes += TagCounter<RunNum>(qOld.getFirst("runTimes"));
	}
	
	// initialize Analysis DB variables
	for(GVState g = GV_CLOSED; g <= GV_OTHER; ++g)
		for(AFPState a = AFP_OFF; a <= AFP_OTHER; ++a)
			AR_IDs[g][a] = 0;
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

BlindTime RunAccumulator::getTotalTime(AFPState afp, GVState gv) const {
	smassert(gv==GV_OPEN || gv==GV_CLOSED);
	if(afp==AFP_OTHER)
		return totalTime[AFP_ON][gv]+totalTime[AFP_OFF][gv]+totalTime[AFP_OTHER][gv];
	return totalTime[afp][gv];
}

BlindTime RunAccumulator::getTotalTime(GVState gv) const {
	return getTotalTime(AFP_OFF,gv) + getTotalTime(AFP_ON,gv) + getTotalTime(AFP_OTHER,gv);
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
	printf("Calculating results for %s...\n",name.c_str());
	if(isCalculated) printf("*** Warning: repeat calculation!\n");
	for(std::map<std::string,AnalyzerPlugin*>::iterator it = myPlugins.begin(); it != myPlugins.end(); it++) {
		printf("... results in '%s' ...\n",it->first.c_str());
		it->second->calculateResults();
	}
	printf("Done calculating results for %s.\n",name.c_str());
	isCalculated = true;
}

void RunAccumulator::uploadAnaNumber(AnaNumber& AN, GVState g, AFPState a) {
	assert(g<=GV_OTHER && a <= AFP_OTHER);
	AN.date = time(NULL);
	AN.source = getEnvSafe("UCNA_ANA_AUTHOR") + (isSimulated?"_Sim":"_Data") + AnaDB_xtag;
	anaResults[g][a].push_back(AN);
}

void RunAccumulator::uploadAnaResults() {
	for(GVState gv = GV_CLOSED; gv <= GV_OTHER; ++gv) {
		for(AFPState afp = AFP_OFF; afp <= AFP_OTHER; ++afp) {
			if(!anaResults[gv][afp].size()) continue;
			unsigned int rsid = get_ADB_Runset_ID(gv,afp);
			if(!rsid) continue;
			for(std::vector<AnaNumber>::iterator it = anaResults[gv][afp].begin(); it != anaResults[gv][afp].end(); it++) {
				it->rsid = rsid;
				it->anid = 0;
				printf("Analysis result '%s' [%i]:\t%g\t~ %g\n", it->name.c_str(), it->n, it->value, it->err);
				AnalysisDB::getADB()->uploadAnaNumber(*it);
			}
		}
	}
}

unsigned int RunAccumulator::get_ADB_Runset_ID(GVState g, AFPState a) {
	AnalysisDB* ADB = AnalysisDB::getADB();
	if(!ADB) return 0;
	
	smassert(g<=GV_OTHER && a<=AFP_OTHER);
	if(!AR_IDs[g][a]) {
		AnaRunset AR;
		if(!runCounts.counts.empty()) {
			AR.start_run = runCounts.counts.begin()->first;
			AR.end_run = runCounts.counts.rbegin()->first;
		}
		AR.gate_valve = g;
		AR.afp = a;
		AR.grouping = grouping;
		ADB->locateRunset(AR);
		AR_IDs[g][a] = AR.rsid;
	}
	assert(AR_IDs[g][a]);
	return AR_IDs[g][a];
}

void RunAccumulator::makePlots() {
	defaultCanvas->cd();
	if(!isCalculated) calculateResults();
	printf("Generating plots for %s...\n",name.c_str());
	for(std::map<std::string,AnalyzerPlugin*>::iterator it = myPlugins.begin(); it != myPlugins.end(); it++) {
		printf("... plots in '%s' ...\n",it->first.c_str());
		it->second->makePlots();
	}
	printf("Done generating plots for %s...\n",name.c_str());
}

void RunAccumulator::compareMCtoData(RunAccumulator& OAdata) {
	defaultCanvas->cd();
	if(!isCalculated) calculateResults();
	if(!OAdata.isCalculated) OAdata.calculateResults();
	printf("Comparing MC %s and data %s...\n",name.c_str(),OAdata.name.c_str());
	for(std::map<std::string,AnalyzerPlugin*>::iterator it = myPlugins.begin(); it != myPlugins.end(); it++) {
		AnalyzerPlugin* AP = OAdata.getPlugin(it->second->name);
		if(AP) {
			printf("... comparison in '%s' ...\n",it->first.c_str());
			it->second->compareMCtoData(AP);
		}
	}
	printf("Done comparing MC %s and data %s.\n",name.c_str(),OAdata.name.c_str());
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
	smassert(it != fgbgHists.end());
	return *(it->second);
}

const fgbgPair& RunAccumulator::getFGBGPair(const std::string& qname) const {
	std::map<std::string,fgbgPair*>::const_iterator it = fgbgHists.find(qname);
	smassert(it != fgbgHists.end());
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
	smassert(lowrate && master);
	if(!master->Integral()) return;
	float errscaling = sqrt(lowrate->Integral()/master->Integral());
	smassert(errscaling==errscaling);
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
		double bgRatio = RefOA.getTotalTime(it->second->afp,GV_OPEN)[it->second->mySide]/RefOA.getTotalTime(it->second->afp,GV_CLOSED)[it->second->mySide];
		for(unsigned int i=0; i<totalBins(it->second->h[GV_CLOSED]); i++) {
			double rootn = qhRef.h[0]->GetBinError(i)*sqrt(simfactor);		// root(bg counts) from ref histogram errorbars
			double n = rootn*rootn;											// background counts from reference histogram
			double bgObsCounts = rnd_source.PoissonD(n);					// simulated background counts
			double fgBgCounts = rnd_source.PoissonD(n*bgRatio);				// simulated foreground counts due to background
			if(!addFluctCounts) bgObsCounts = fgBgCounts = 0.;
			it->second->h[GV_CLOSED]->AddBinContent(i,bgObsCounts);			// fill fake background counts
			it->second->h[GV_CLOSED]->SetBinError(i,rootn);					// set root-n statitics for background
			it->second->h[GV_OPEN]->AddBinContent(i,fgBgCounts);			// add simulated background to foreground histogram
		}
		printf("\t%i counts for %s [%i bins]\n",int(it->second->h[0]->Integral()),it->second->getName().c_str(),totalBins(it->second->h[0]));
		it->second->h[GV_OPEN]->Add(it->second->h[GV_CLOSED],-bgRatio);		// subtract back off simulated background
		it->second->isSubtracted = true;
	}
}

TH1* RunAccumulator::rateHisto(const fgbgPair* p, GVState gv) const {
	smassert(p);
	smassert(gv == GV_OPEN || gv == GV_CLOSED);
	TH1* h = (TH1*)p->h[gv]->Clone();
	h->SetTitle(p->getTitle().c_str());
	Side s = p->mySide;
	h->Scale(1.0/h->GetXaxis()->GetBinWidth(1)/getTotalTime(gv)[s]);
	return h;
}

void RunAccumulator::write(std::string outName) {
	printf("Writing data to file '%s'...\n",outName.c_str());
	
	// clear previous tallies
	qOut.erase("totalTime");
	qOut.erase("totalCounts");
	qOut.erase("runCounts");
	qOut.erase("runTimes");
	
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
	smassert(afp <= AFP_OTHER);
	smassert(gv==GV_CLOSED || gv==GV_OPEN);
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

void RunAccumulator::copyTimes(const RunAccumulator& RA) {
	runTimes = RA.runTimes;
	for(AFPState afp = AFP_OFF; afp<=AFP_OTHER; ++afp)
		for(GVState gv = GV_CLOSED; gv <= GV_OPEN; ++gv)
			totalTime[afp][gv] = RA.totalTime[afp][gv];
}

void RunAccumulator::loadSimPoint(Sim2PMT& simData) {
	if(!simData.physicsWeight) return;
	fillCoreHists(simData,simData.physicsWeight);
	if(double evtc = simData.simEvtCounts()) {
		runCounts.add(simData.getRun(),evtc);
		totalCounts[currentAFP][GV_OPEN] += evtc;
	}
}

void RunAccumulator::loadSimData(Sim2PMT& simData, unsigned int nToSim, bool countAll) {
  //printf("Got here!\n");
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

void RunAccumulator::simForRun(Sim2PMT& simData, RunNum rn, unsigned int nToSim, bool countAll) {
	RunInfo RI = CalDBSQL::getCDB()->getRunInfo(rn);
	if(RI.gvState != GV_OPEN) { printf("Skipping simulation for background run "); RI.display(); return; }
	smassert(RI.afpState <= AFP_OTHER);
	
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
		smassert(nRequested);
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

void RunAccumulator::makeOutput(bool doPlots) {
	if(needsSubtraction)
		bgSubtractAll();
	calculateResults();
	uploadAnaResults();
	if(doPlots)
		makePlots();
	//getPlugin("mySimAsym")->printSimTree();
	write();
	setWriteRoot(true);
}

unsigned int RunAccumulator::mergeDir() {
	std::vector<std::string> fnames = listdir(basePath);
	unsigned int nMerged = 0;
	for(std::vector<std::string>::iterator it = fnames.begin(); it != fnames.end(); it++) {
		// check whether data directory contains cloneable subdirectories
		std::string datinfl = basePath+"/"+(*it)+"/"+(*it);
		if(!inflExists(datinfl)) continue;
		SegmentSaver* subRA = makeAnalyzer(*it,datinfl);
		addSegment(*subRA);
		delete(subRA);
		nMerged++;
	}
	makeOutput();
	return nMerged;
}


void RunAccumulator::mergeOcts(const std::vector<Octet>& Octs) {
	for(std::vector<Octet>::const_iterator octit = Octs.begin(); octit != Octs.end(); octit++) {
		std::string inflname = basePath+"/"+octit->octName()+"/"+octit->octName();
		if(!inflExists(inflname)) {
			printf("Octet '%s' missing!\n",inflname.c_str());
			continue;
		}
		SegmentSaver* subRA = makeAnalyzer(octit->octName(),inflname);
		addSegment(*subRA);
		delete(subRA);
		qOut.insert("Octet",octit->toStringmap());
	}
	makeOutput();
}


void RunAccumulator::mergeSims(const std::string& basedata, RunAccumulator* origRA) {
	std::vector<std::string> fnames = listdir(basedata);
	for(std::vector<std::string>::iterator it = fnames.begin(); it != fnames.end(); it++) {
		// check whether data directory contains cloneable subdirectories
		std::string datinfl = basedata+"/"+(*it)+"/"+(*it);
		if(!inflExists(datinfl)) continue;
		std::string siminfl = basePath+"/"+(*it)+"/"+(*it);
		if(!inflExists(siminfl)) { printf("*** Missing simulation for '%s'!\n",siminfl.c_str()); continue; }
		// load cloned sub-data
		SegmentSaver* subRA = makeAnalyzer(*it,siminfl);
		addSegment(*subRA);
		delete(subRA);
	}
	makeOutput();
	if(origRA) {
		origRA->calculateResults();
		compareMCtoData(*origRA);
	}
	write();
	setWriteRoot(true);
}


unsigned int RunAccumulator::simuClone(const std::string& basedata, Sim2PMT& simData, double simfactor, double replaceIfOlder, bool doPlots, bool doCompare) {
	
	printf("\n------ Cloning data in '%s'\n------              to '%s'...\n",basedata.c_str(),basePath.c_str());
	int nCloned = 0;
	
	// check if simulation is already up-to-date
	if(getInflAge() && getInflAge() < replaceIfOlder) {
		printf("\tSimulations in '%s' already recently generated, update skipped...\n",basePath.c_str());
		return nCloned;
	}
	
	// clear any existing data from out-of-date input
	zeroSavedHists();
	zeroCounters();
	
	// load original data for comparison
	std::vector<std::string> datpath = split(strip(basedata,"/"),"/");
	smassert(datpath.size()>0);
	isSimulated = false;
	RunAccumulator* origRA = (RunAccumulator*)makeAnalyzer("nameUnused",basedata+"/"+datpath.back());
	isSimulated = true;
	// copy over octet information
	qOut.erase("Octet");
	qOut.transfer(origRA->qOut,"Octet");
	
	// load/clone data in all subdirectories, if they exist
	int nClonable = 0;
	std::vector<std::string> fnames = listdir(basedata);
	for(std::vector<std::string>::iterator it = fnames.begin(); it != fnames.end(); it++) {
		// check whether data directory contains cloneable subdirectories
		std::string datinfl = basedata+"/"+(*it)+"/"+(*it);
		if(!RunAccumulator::inflExists(datinfl)) continue;
		std::string siminfl = basePath+"/"+(*it)+"/"+(*it);
		nClonable++;
		
		// load cloned sub-data
		RunAccumulator* subRA = (RunAccumulator*)makeAnalyzer(*it,RunAccumulator::inflExists(siminfl)?siminfl:"");
		subRA->grouping = subdivide(grouping);
		subRA->simPerfectAsym = simPerfectAsym;
		nCloned += subRA->simuClone(basedata+"/"+(*it),simData,simfactor,replaceIfOlder,doPlots,doCompare);
		addSegment(*subRA);
		delete(subRA);
	}
	
	// if there were no clonable subdrectories, we need to actually simulate the data
	if(!nClonable) {
		printf("\tNo data subdirectories found; assume data here needs cloning...\n");
		printf("\tProcessing simulation data for each pulse-pair run...\n");
		
		// collate requested simulation counts
		TagCounter<RunNum> countRequests[2];
		for(std::map<RunNum,double>::iterator it = origRA->runCounts.counts.begin(); it != origRA->runCounts.counts.end(); it++) {
			nCloned++;
			if(!it->first || !it->second) continue;
			RunInfo RI = CalDBSQL::getCDB()->getRunInfo(it->first);
			// no simulation for background runs
			if(RI.gvState != GV_OPEN) {
				runCounts.add(it->first,0);
				runTimes.add(it->first,0);
				continue;
			}
			smassert(RI.afpState <= AFP_ON);	// are you really trying to clone non-beta-octet runs here??
			// estimate background count share for this run (and reduce simulation by this amount)
			double bgEst = origRA->getTotalCounts(RI.afpState,GV_CLOSED)*origRA->getRunTime(it->first)/origRA->getTotalTime(RI.afpState,GV_CLOSED)[BOTH];
			if(it->second <= bgEst) continue;
			double nToSim = simfactor*(it->second-bgEst);
			printf("\t---Simulation cloning for run %i (%i+%i counts)---\n",it->first,int(nToSim),int(simfactor*bgEst));
			countRequests[RI.afpState].add(it->first,nToSim);
		}
		
		// see which kind of run got more AFP counts; clone first
		AFPState bigafp = countRequests[AFP_OFF].total()>=countRequests[AFP_ON].total()?AFP_OFF:AFP_ON;
		AFPState smallafp = bigafp?AFP_OFF:AFP_ON;
		smassert(countRequests[bigafp].nTags());
		// determine starting point in simulation data to which we can return later
		unsigned int startEvt = 0;
		while(!startEvt) {
			simData.startScan(true);
			startEvt = simData.getCurrentEvent();
		}
		unsigned int rseed = countRequests[bigafp].counts.begin()->first;
		PMTGenerator::sim_rnd_source.SetSeed(rseed);
		RunAccumulator::rnd_source.SetSeed(rseed);
		unsigned int nSimmed = simMultiRuns(simData, countRequests[bigafp]);
		// return to starting point or continue on with different data
		if(simPerfectAsym) {
			PMTGenerator::sim_rnd_source.SetSeed(rseed);
			RunAccumulator::rnd_source.SetSeed(rseed);
			simData.gotoEvent(startEvt);
			simMultiRuns(simData, countRequests[smallafp], nSimmed);
		} else {
			simMultiRuns(simData, countRequests[smallafp]);
		}
		
		// clone background counts in original data
		simBgFlucts(*origRA,simfactor,!simPerfectAsym);
		// scale back out simulation factor
		scaleData(1.0/simfactor);
	}
	
	// generate output
	makeOutput(doPlots);
	if(doCompare) {
		origRA->calculateResults();
		compareMCtoData(*origRA);
	}
	write();
	
	delete(origRA);
	return nCloned;
}

/* --------------------------------------------------- */

fgbgPair* AnalyzerPlugin::registerFGBGPair(const std::string& hname, const std::string& title,
										   unsigned int nbins, float xmin, float xmax,
										   AFPState a, Side s) {
	myHists.push_back(myA->registerFGBGPair(hname,title,nbins,xmin,xmax,a,s));
	return myHists.back();
}

fgbgPair* AnalyzerPlugin::registerFGBGPair(const TH1& hTemplate, AFPState a, Side s) {
	myHists.push_back(myA->registerFGBGPair(hTemplate,a,s));
	return myHists.back();
}

void AnalyzerPlugin::makeRatesSummary() {
	for(std::vector<fgbgPair*>::const_iterator it = myHists.begin(); it != myHists.end(); it++) {
		for(GVState gv=GV_CLOSED; gv<=GV_OPEN; ++gv) {
			
			Double_t d_counts;
			double counts = (*it)->h[gv]->IntegralAndError(-1, -1, d_counts);
			double rtime = myA->getTotalTime((*it)->afp,gv)[(*it)->mySide];
			
			AnaNumber AN("counts_"+(*it)->baseName);
			AN.s = (*it)->mySide;
			
			AN.value = counts;
			AN.err = d_counts;
			myA->uploadAnaNumber(AN,gv,(*it)->afp);
			
			AN.name = "rate_"+(*it)->baseName;
			AN.value = counts?counts/rtime:0;
			AN.err = d_counts?d_counts/rtime:0;
			myA->uploadAnaNumber(AN,gv,(*it)->afp);
		}
	}
}

/* --------------------------------------------------- */

unsigned int recalcOctets(RunAccumulator& RA, const std::vector<Octet>& Octs, bool doPlots) {
	unsigned int nproc = 0;
	
	for(std::vector<Octet>::const_iterator octit = Octs.begin(); octit != Octs.end(); octit++) {
		
		printf("Reprocessing octet '%s'...\n",octit->octName().c_str());
		if(!octit->getNRuns()) {
			printf("\tThat was too easy (octet contained zero runs).\n");
			continue;
		}
		RA.qOut.insert("Octet",octit->toStringmap());
		
		std::string inflname = RA.basePath+"/"+octit->octName()+"/"+octit->octName();
		if(SegmentSaver::inflExists(inflname)) {
			RunAccumulator* subRA = (RunAccumulator*)RA.makeAnalyzer(octit->octName(),inflname);
			subRA->grouping = octit->grouping;
			nproc += recalcOctets(*subRA,octit->getSubdivs(subdivide(octit->grouping),false), doPlots);
			RA.addSegment(*subRA);
			delete subRA;
		}
	}
	
	RA.makeOutput(doPlots);
	RA.setWriteRoot(false);
	return nproc;
}

unsigned int processOctets(RunAccumulator& RA, const std::vector<Octet>& Octs, double replaceIfOlder, bool doPlots) {
	
	unsigned int nproc = 0;
	if(Octs.size()==1) RA.grouping = Octs[0].grouping;
	
	for(std::vector<Octet>::const_iterator octit = Octs.begin(); octit != Octs.end(); octit++) {
		
		// check if there are any runs to process
		printf("Processing beta runs '%s'...\n",octit->octName().c_str());
		if(!octit->getNRuns()) {
			printf("\tThat was too easy (octet contained zero runs).\n");
			continue;
		}
		RA.qOut.insert("Octet",octit->toStringmap());
		
		if(octit->grouping > GROUP_FGBG) {
			// make sub-Analyzer for this octet, to load data if already available, otherwise re-process
			RunAccumulator* subRA;
			std::string inflname = RA.basePath+"/"+octit->octName()+"/"+octit->octName();
			double fAge = fileAge(inflname+".root");
			if(SegmentSaver::inflExists(inflname) && fAge < replaceIfOlder) {
				printf("Octet '%s' already scanned %.1fh ago; skipping\n",octit->octName().c_str(),fAge/3600);
				subRA = (RunAccumulator*)RA.makeAnalyzer(octit->octName(),inflname);
			} else {
				subRA = (RunAccumulator*)RA.makeAnalyzer(octit->octName(),"");
				subRA->grouping = octit->grouping;
				nproc += processOctets(*subRA,octit->getSubdivs(subdivide(octit->grouping),false),replaceIfOlder, doPlots);
			}
			RA.addSegment(*subRA);
			delete(subRA);
		} else {
			// add data from FG, BG runs pair
			for(GVState gv=GV_CLOSED; gv<=GV_OPEN; ++gv) {
				PostOfficialAnalyzer PDS(true);
				PDS.addRuns(octit->getAsymRuns(gv));
				if(PDS.getnFiles())
					RA.loadProcessedData(octit->octAFPState(), gv, PDS);
			}
			nproc++;
		}
	}
	
	RA.makeOutput(doPlots);
	return nproc;
}

