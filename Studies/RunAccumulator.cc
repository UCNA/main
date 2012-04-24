#include "RunAccumulator.hh"
#include <TRandom3.h>

void fgbgPair::bgSubtract(BlindTime tFG, BlindTime tBG) {
	assert(!isSubtracted); // don't BG subtract twice!
	double bgScale = tBG.t[mySide]?tFG.t[mySide]/tBG.t[mySide]:1.0;
	h[1]->Add(h[0],-bgScale);
	isSubtracted = true;
}

void fgbgPair::operator+=(const fgbgPair& p) {
	for(unsigned int fg=0; fg<=1; fg++) {
		assert(h[fg] && p.h[fg]);
		h[fg]->Add(p.h[fg]);
	}
}

void fgbgPair::operator*=(double c) {
	for(unsigned int fg=0; fg<=1; fg++) {
		assert(h[fg]);
		h[fg]->Scale(c);
	}
}



fgbgPair RunAccumulator::registerFGBGPair(const std::string& hname, const std::string& title,
										  unsigned int nbins, float xmin, float xmax, AFPState a, Side s) {
	fgbgPair p(hname,a,s);
	assert(fgbgHists.find(p.getName())==fgbgHists.end()); // don't duplicate names!
	for(unsigned int fg = 0; fg <= 1; fg++) {
		p.h[fg] = registerSavedHist(p.getHistoName(fg),title+(fg?"":" Background"),nbins,xmin,xmax);
		if(!p.h[fg]->GetSumw2N())
			p.h[fg]->Sumw2();
	}
	fgbgHists.insert(std::make_pair(p.getName(),p));
	return p;
}

fgbgPair RunAccumulator::registerFGBGPair(const TH1& hTemplate, AFPState a, Side s) {
	std::string hname = hTemplate.GetName();
	fgbgPair p(hname,a,s);
	assert(fgbgHists.find(p.getName())==fgbgHists.end()); // don't duplicate names!
	for(int fg = 1; fg >=0; fg--) {
		p.h[fg] = registerSavedHist(p.getHistoName(fg),hTemplate);
		if(!p.h[fg]->GetSumw2N())
			p.h[fg]->Sumw2();
		if(!fg)
			p.h[0]->SetTitle((std::string(p.h[1]->GetTitle())+" Background").c_str());
	}
	fgbgHists.insert(std::make_pair(p.getName(),p));	
	return p;
}

fgbgPair RunAccumulator::cloneFGBGPair(const fgbgPair& p, const std::string& newName) {
	fgbgPair pnew(newName,p.afp,p.mySide);
	for(unsigned int fg = 0; fg <= 1; fg++) {
		std::string qname = pnew.getHistoName(fg);
		pnew.h[fg] = (TH1*)addObject(p.h[fg]->Clone(qname.c_str()));
	}
	return pnew;
}

RunAccumulator::RunAccumulator(OutputManager* pnt, const std::string& nm, const std::string& inflName):
SegmentSaver(pnt,nm,inflName), needsSubtraction(false) {
	
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
			int afp = int(it->getDefault("afp",3));
			int fg = int(it->getDefault("fg",3));
			assert(afp<=AFP_OTHER && fg <= 1);
			totalTime[afp][fg] = BlindTime(*it);
		}
		// fetch total counts
		std::vector<Stringmap> counts = qOld.retrieve("totalCounts");
		for(std::vector<Stringmap>::iterator it = counts.begin(); it != counts.end(); it++) {
			int afp = int(it->getDefault("afp",3));
			int fg = int(it->getDefault("fg",3));
			assert(afp<=AFP_OTHER && fg <= 1);
			totalCounts[afp][fg] = it->getDefault("counts",0);
		}
		// fetch run counts, run times
		runCounts += TagCounter<RunNum>(qOld.getFirst("runCounts"));
		runTimes += TagCounter<RunNum>(qOld.getFirst("runTimes"));
	}
}

void RunAccumulator::zeroCounters() {
	for(unsigned int afp = AFP_OFF; afp <= AFP_OTHER; afp++) {
		for(unsigned int fg = 0; fg <= 1; fg++) {
			totalTime[afp][fg] = BlindTime(0.0);
			totalCounts[afp][fg] = 0;
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
	// add times, counts
	for(unsigned int afp = AFP_OFF; afp <= AFP_OTHER; afp++) {
		for(unsigned int fg = 0; fg <= 1; fg++) {
			totalTime[afp][fg] += RA.totalTime[afp][fg];
			totalCounts[afp][fg] += RA.totalCounts[afp][fg];
		}
	}
	// add run counts, times
	runCounts += RA.runCounts;
	runTimes += RA.runTimes;
	// transfer run calibration data
	qOut.transfer(S.qOut,"runcal");
}

fgbgPair& RunAccumulator::getFGBGPair(const std::string& qname) {
	std::map<std::string,fgbgPair>::iterator it = fgbgHists.find(qname);
	assert(it != fgbgHists.end());
	return it->second;
}

const fgbgPair& RunAccumulator::getFGBGPair(const std::string& qname) const {
	std::map<std::string,fgbgPair>::const_iterator it = fgbgHists.find(qname);
	assert(it != fgbgHists.end());
	return it->second;
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
	for(std::map<std::string,fgbgPair>::iterator it = fgbgHists.begin(); it != fgbgHists.end(); it++) {
		if(getErrorEstimator())
			errorbarsFromMasterHisto(it->second.h[0],getErrorEstimator()->getFGBGPair(it->second.getName()).h[0]);
		it->second.bgSubtract(totalTime[it->second.afp][1],totalTime[it->second.afp][0]);
	}
	needsSubtraction = false;
}

void RunAccumulator::simBgFlucts(const RunAccumulator& RefOA, double simfactor) {
	assert(isEquivalent(RefOA));
	printf("Adding background fluctuations to simulation...\n");
	static TRandom3* OA_sim_rnd_source = NULL;
	if(!OA_sim_rnd_source) {
		OA_sim_rnd_source = new TRandom3();
		OA_sim_rnd_source->SetSeed(0);
	}
	for(std::map<std::string,fgbgPair>::iterator it = fgbgHists.begin(); it != fgbgHists.end(); it++) {
		fgbgPair qhRef = RefOA.getFGBGPair(it->first);
		double bgRatio = RefOA.getTotalTime(it->second.afp,true).t[it->second.mySide]/RefOA.getTotalTime(it->second.afp,false).t[it->second.mySide];
		for(unsigned int i=0; i<totalBins(it->second.h[0]); i++) {
			double rootn = qhRef.h[0]->GetBinError(i)*sqrt(simfactor);		// root(bg counts) from ref histogram errorbars
			double n = rootn*rootn;											// background counts from reference histogram
			double bgObsCounts = OA_sim_rnd_source->PoissonD(n);			// simulated background counts
			double fgBgCounts = OA_sim_rnd_source->PoissonD(n*bgRatio);		// simulated foreground counts due to background
			it->second.h[0]->AddBinContent(i,bgObsCounts);					// fill fake background counts
			it->second.h[0]->SetBinError(i,rootn);							// set root-n statitics for background
			it->second.h[1]->AddBinContent(i,fgBgCounts);					// add simulated background to foreground histogram				
		}
		printf("\t%i counts for %s [%i bins]\n",int(it->second.h[0]->Integral()),it->second.getName().c_str(),totalBins(it->second.h[0]));
		it->second.h[1]->Add(it->second.h[0],-bgRatio);						// subtract back off simulated background
		it->second.isSubtracted = true;
	}
}

void RunAccumulator::makeRatesSummary() {
	for(std::map<std::string,fgbgPair>::const_iterator it = fgbgHists.begin(); it != fgbgHists.end(); it++) {
		for(unsigned int fg = 0; fg <= 1; fg++) {
			Stringmap rt;
			rt.insert("side",ctos(sideNames(it->second.mySide)));
			rt.insert("afp",itos(it->second.afp));
			rt.insert("name",it->second.getName());
			rt.insert("fg",itos(fg));
			double counts = it->second.h[fg]->Integral();
			rt.insert("counts",counts);
			rt.insert("rate",counts/totalTime[it->second.afp][fg].t[it->second.mySide]);
			qOut.insert("rate",rt);
		}
	}
}

void RunAccumulator::write(std::string outName) {
	// record total times, counts
	for(unsigned int afp = AFP_OFF; afp <= AFP_OTHER; afp++) {
		for(unsigned int fg = 0; fg <= 1; fg++) {
			Stringmap tm = totalTime[afp][fg].toStringmap();
			tm.insert("afp",itos(afp));
			tm.insert("fg",itos(fg));
			qOut.insert("totalTime",tm);
			Stringmap ct;
			ct.insert("afp",itos(afp));
			ct.insert("fg",itos(fg));
			ct.insert("counts",totalCounts[afp][fg]);
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
	currentGV = gv;
	currentAFP = afp;
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

void RunAccumulator::loadSimData(Sim2PMT& simData, unsigned int nToSim) {
	currentGV = GV_OPEN;
	currentAFP = simData.getAFP();
	printf("Loading %i events of simulated data (AFP=%i)...\n",nToSim,currentAFP);
	simData.startScan(nToSim);
	while(simData.nSimmed<=nToSim) {
		simData.nextPoint();
		loadSimPoint(simData);
		if(!(int(simData.nSimmed)%(nToSim/20))) {
			if(nToSim>1e6) {
				printf("* %s\n",simData.evtInfo().toString().c_str());
			} else {
				printf("*");
				fflush(stdout);
			}
		}
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
