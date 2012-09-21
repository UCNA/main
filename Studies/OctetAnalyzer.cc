#include "OctetAnalyzer.hh"
#include "PathUtils.hh"
#include "strutils.hh"
#include "GraphicsUtils.hh"
#include "EnergyCalibrator.hh"
#include "PostOfficialAnalyzer.hh"
#include "CalDBSQL.hh"
#include "Types.hh"

void quadHists::setFillPoint(AFPState afp, GVState gv) {
	assert(gv==GV_CLOSED || gv==GV_OPEN);
	assert(afp==AFP_OFF||afp==AFP_ON);
	fillPoint = fgbg[afp]->h[gv];
	assert(fillPoint);
}

bool quadHists::isEquivalent(const quadHists& qh) const {
	if(name != qh.name) return false;
	return true;
	//return (h[AFP_OFF][1]->GetNbinsX() == qh.h[AFP_OFF][1]->GetNbinsX()
	//		&& h[AFP_OFF][1]->GetNbinsY() == qh.h[AFP_OFF][1]->GetNbinsY()
	//		&& h[AFP_OFF][1]->GetNbinsZ() == qh.h[AFP_OFF][1]->GetNbinsZ());	//< TODO more rigorous check
}

void quadHists::operator+=(const quadHists& qh) {
	for(AFPState afp = AFP_OFF; afp <= AFP_ON; ++afp)
		*fgbg[afp] += *qh.fgbg[afp];
}

void quadHists::operator*=(double c) {
	for(AFPState afp = AFP_OFF; afp <= AFP_ON; ++afp)
		*fgbg[afp] *= c;
}

void quadHists::setDrawRange(double y, bool maximum) {
	for(AFPState afp = AFP_OFF; afp <= AFP_ON; ++afp)
		for(GVState gv=GV_CLOSED; gv<=GV_OPEN; ++gv)
			if(maximum)
				fgbg[afp]->h[gv]->SetMaximum(y);
			else
				fgbg[afp]->h[gv]->SetMinimum(y);
}

void quadHists::setRangeUser(double mn, double mx, AxisDirection d) {
	for(AFPState afp = AFP_OFF; afp <= AFP_ON; ++afp)
		for(GVState gv=GV_CLOSED; gv<=GV_OPEN; ++gv)
			(d==X_DIRECTION?fgbg[afp]->h[gv]->GetXaxis():fgbg[afp]->h[gv]->GetYaxis())->SetRangeUser(mn,mx);
}

void quadHists::setAxisTitle(AxisDirection d, const std::string& ttl) {
	for(AFPState afp = AFP_OFF; afp <= AFP_ON; ++afp)
		fgbg[afp]->setAxisTitle(d,ttl);
}

void quadHists::setSubtraction(bool b) {
	for(AFPState afp = AFP_OFF; afp <= AFP_ON; ++afp)
		fgbg[afp]->doSubtraction = b;
}

void quadHists::setTimeScaling(bool b) {
	for(AFPState afp = AFP_OFF; afp <= AFP_ON; ++afp)
		fgbg[afp]->doTimeScale = b;
}

/* --------------------------------------------------- */

OctetAnalyzer::~OctetAnalyzer() {
	for(std::map<std::string,quadHists*>::iterator it = coreHists.begin(); it != coreHists.end(); it++)
		delete it->second;
}

quadHists* OctetAnalyzer::registerCoreHist(const std::string& hname, const std::string& title,
										  unsigned int nbins, float xmin, float xmax, Side s) {
	quadHists* qh = new quadHists(hname,title,s);
	assert(coreHists.find(qh->getName())==coreHists.end());	// don't duplicate names!
	for(AFPState afp = AFP_OFF; afp <= AFP_ON; ++afp)
		qh->fgbg[afp] = registerFGBGPair(hname,title,nbins,xmin,xmax,afp,s);
	coreHists.insert(std::make_pair(qh->getName(),qh));
	return qh;
}

quadHists* OctetAnalyzer::registerCoreHist(const TH1& hTemplate, Side s) {
	quadHists* qh = new quadHists(hTemplate.GetName(),hTemplate.GetTitle(),s);
	assert(coreHists.find(qh->getName())==coreHists.end());	// don't duplicate names!
	for(AFPState afp = AFP_OFF; afp <= AFP_ON; ++afp)
		qh->fgbg[afp] = registerFGBGPair(hTemplate,AFPState(afp),s);
	coreHists.insert(std::make_pair(qh->getName(),qh));
	return qh;
}

quadHists* OctetAnalyzer::cloneQuadHist(const quadHists* qh, const std::string& newName, const std::string& newTitle) {
	assert(qh);
	quadHists* qnew = new quadHists(newName,newTitle,qh->mySide);
	for(AFPState afp = AFP_OFF; afp <= AFP_ON; ++afp)
		qnew->fgbg[afp] = cloneFGBGPair(*qh->fgbg[afp],newName,newTitle);
	coreHists.insert(std::make_pair(qnew->getName(),qnew));
	return qnew;
}

OctetAnalyzer::OctetAnalyzer(OutputManager* pnt, const std::string& nm, const std::string& inflName):
RunAccumulator(pnt,nm,inflName), depth(-1), simPerfectAsym(false) { }

void OctetAnalyzer::setFillPoints(AFPState afp, GVState gv) {
	for(std::map<std::string,quadHists*>::iterator it = coreHists.begin(); it != coreHists.end(); it++)
		it->second->setFillPoint(afp,gv);
}

quadHists* OctetAnalyzer::getCoreHist(const std::string& qname) {
	std::map<std::string,quadHists*>::iterator it = coreHists.find(qname);
	assert(it != coreHists.end());
	return it->second;
}

const quadHists* OctetAnalyzer::getCoreHist(const std::string& qname) const {
	std::map<std::string,quadHists*>::const_iterator it = coreHists.find(qname);
	assert(it != coreHists.end());
	return it->second;
}

void OctetAnalyzer::loadProcessedData(AFPState afp, ProcessedDataScanner& FG, ProcessedDataScanner& BG) {
	assert(afp == AFP_OFF || afp == AFP_ON);
	for(GVState gv = GV_CLOSED; gv <= GV_OPEN; ++gv) {
		ProcessedDataScanner& PDS = gv?FG:BG;
		setFillPoints(afp,gv);
		RunAccumulator::loadProcessedData(afp,gv,PDS);
	}
}

void OctetAnalyzer::loadSimData(Sim2PMT& simData, unsigned int nToSim, bool countAll) {
	setFillPoints(simData.getAFP(),GV_OPEN);
	RunAccumulator::loadSimData(simData,nToSim,countAll);
}

unsigned int OctetAnalyzer::simuClone(const std::string& basedata, Sim2PMT& simData, double simfactor, double replaceIfOlder, bool doPlots) {
	
	printf("\n------ Cloning data in '%s'\n------                        to '%s'...\n",basedata.c_str(),basePath.c_str());
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
	assert(datpath.size()>0);
	isSimulated = false;
	OctetAnalyzer* origOA = (OctetAnalyzer*)makeAnalyzer("nameUnused",basedata+"/"+datpath.back());
	isSimulated = true;
	// copy over octet information
	qOut.erase("Octet");
	qOut.transfer(origOA->qOut,"Octet");
	
	// load/clone data in all subdirectories, if they exist
	int nClonable = 0;
	std::vector<std::string> fnames = listdir(basedata);
	for(std::vector<std::string>::iterator it = fnames.begin(); it != fnames.end(); it++) {
		// check whether data directory contains cloneable subdirectories
		std::string datinfl = basedata+"/"+(*it)+"/"+(*it);
		if(!OctetAnalyzer::inflExists(datinfl)) continue;
		std::string siminfl = basePath+"/"+(*it)+"/"+(*it);		
		nClonable++;
		
		// load cloned sub-data
		OctetAnalyzer* subOA = (OctetAnalyzer*)makeAnalyzer(*it,OctetAnalyzer::inflExists(siminfl)?siminfl:"");
		subOA->depth = depth+1;
		subOA->simPerfectAsym = simPerfectAsym;
		nCloned += subOA->simuClone(basedata+"/"+(*it),simData,simfactor,replaceIfOlder,doPlots);
		addSegment(*subOA);
		delete(subOA);
	}
	
	// if there were no clonable subdrectories, we need to actually simulate the data
	if(!nClonable) {
		printf("\tNo data subdirectories found; assume data here needs cloning...\n");
		printf("\tProcessing simulation data for each pulse-pair run...\n");
		
		// collate requested simulation counts
		TagCounter<RunNum> countRequests[2];
		for(std::map<RunNum,double>::iterator it = origOA->runCounts.counts.begin(); it != origOA->runCounts.counts.end(); it++) {
			nCloned++;
			if(!it->first || !it->second) continue;
			RunInfo RI = CalDBSQL::getCDB()->getRunInfo(it->first);
			// no simulation for background runs
			if(RI.gvState != GV_OPEN) {
				runCounts.add(it->first,0);
				runTimes.add(it->first,0);
				continue;
			}
			assert(RI.afpState <= AFP_ON);	// are you really trying to clone non-beta-octet runs here??
			// estimate background count share for this run (and reduce simulation by this amount)
			double bgEst = origOA->getTotalCounts(RI.afpState,GV_CLOSED)*origOA->getRunTime(it->first)/origOA->getTotalTime(RI.afpState,GV_CLOSED)[BOTH];
			if(it->second <= bgEst) continue;
			double nToSim = simfactor*(it->second-bgEst);
			printf("\t---Simulation cloning for run %i (%i+%i counts)---\n",it->first,int(nToSim),int(simfactor*bgEst));
			countRequests[RI.afpState].add(it->first,nToSim);
		}
		
		// see which kind of run got more AFP counts; clone first
		AFPState bigafp = countRequests[AFP_OFF].total()>=countRequests[AFP_ON].total()?AFP_OFF:AFP_ON;
		AFPState smallafp = bigafp?AFP_OFF:AFP_ON;
		assert(countRequests[bigafp].nTags());
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
		simBgFlucts(*origOA,simfactor,!simPerfectAsym);
		// scale back out simulation factor
		scaleData(1.0/simfactor);
	}
	
	// generate output
	calculateResults();
	uploadAnaResults();
	if(doPlots) {
		makePlots();
		origOA->calculateResults();
		compareMCtoData(*origOA);
	}
	write();
	setWriteRoot(true);
	
	delete(origOA);
	return nCloned;
}

/* --------------------------------------------------- */

TH1* OctetAnalyzerPlugin::calculateSR(const std::string& hname, const quadHists* qEast, const quadHists* qWest, bool fg, bool instr) {
	assert(qEast && qWest);
	// first calculate R = E+W-/E-W+
	TH1* hR = (TH1*)qEast->fgbg[AFP_ON]->h[fg]->Clone("SR_intermediate_R");
	if(instr) {
		hR->Multiply(qEast->fgbg[AFP_OFF]->h[fg]);
	} else {
		hR->Multiply(qWest->fgbg[AFP_OFF]->h[fg]);
		hR->Scale(1.0/(myA->totalTime[AFP_ON][fg][EAST]*myA->totalTime[AFP_OFF][fg][WEST]));
	}
	
	TH1* hAsym = (TH1*)qWest->fgbg[AFP_ON]->h[fg]->Clone(hname.c_str());
	if(instr) {
		hAsym->Multiply(qWest->fgbg[AFP_OFF]->h[fg]);
	} else {
		hAsym->Multiply(qEast->fgbg[AFP_OFF]->h[fg]);
		hAsym->Scale(1.0/(myA->totalTime[AFP_OFF][fg][EAST]*myA->totalTime[AFP_ON][fg][WEST]));
	}
	
	hR->Divide(hAsym);
	
	// super-ratio asymmetry (1-sqrt(R))/(1+sqrt(R))
	for(unsigned int i=0; i<totalBins(hR); i++) {
		double r = hR->GetBinContent(i);
		double dr = hR->GetBinError(i);
		if(r>0 && r==r) {
			hAsym->SetBinContent(i,(1.0-sqrt(r))/(1.0+sqrt(r)));
			hAsym->SetBinError(i, dr/(sqrt(r)*(1+sqrt(r))*(1+sqrt(r))) * (myA->simPerfectAsym && !instr?0.33:1.0));
		} else {
			hAsym->SetBinContent(i,0);
			hAsym->SetBinError(i,1e6);
		}
	}
	
	hAsym->SetLineColor(1);
	hAsym->SetLineStyle(1);
	hAsym->SetTitle((qEast->title+" Asymmetry").c_str());
	hAsym->GetXaxis()->SetTitle(qEast->fgbg[AFP_ON]->h[fg]->GetXaxis()->GetTitle());
	
	delete(hR);
	return (TH1*)myA->addObject(hAsym);
}

/// square root of a histogram
void sqrtHist(TH1* h) {
	for(unsigned int i=0; i<totalBins(h); i++) {
		double y = h->GetBinContent(i);
		double dy = h->GetBinError(i);
		if(y>0 && y==y) {
			h->SetBinContent(i,sqrt(y));
			h->SetBinError(i,dy/(2*sqrt(y)));
		} else {
			h->SetBinContent(i,0);
			h->SetBinError(i,0);
		}
	}
}

TH1* OctetAnalyzerPlugin::calculateSuperSum(const std::string& hname, const quadHists* qEast, const quadHists* qWest, bool fg) {
	assert(qEast && qWest);
	TH1* hR = (TH1*)qEast->fgbg[AFP_ON]->h[fg]->Clone("SuperSum_intermediate");
	hR->Multiply(qWest->fgbg[AFP_OFF]->h[fg]);
	hR->Scale(1.0/(myA->totalTime[AFP_ON][fg][EAST]*myA->totalTime[AFP_OFF][fg][WEST]));
	
	TH1* hSS = (TH1*)qEast->fgbg[AFP_OFF]->h[fg]->Clone(hname.c_str());
	hSS->Multiply(qWest->fgbg[AFP_ON]->h[fg]);
	hSS->Scale(1.0/(myA->totalTime[AFP_OFF][fg][EAST]*myA->totalTime[AFP_ON][fg][WEST]));
	
	sqrtHist(hR);
	sqrtHist(hSS);
	hSS->Add(hR);
	hSS->Scale(0.5);
	
	hSS->SetTitle((qEast->title+" SuperSum").c_str());
	hSS->SetLineColor(1);
	hSS->SetLineStyle(1);
	hSS->GetXaxis()->SetTitle(qEast->fgbg[AFP_ON]->h[fg]->GetXaxis()->GetTitle());
	
	delete(hR);
	return (TH1*)myA->addObject(hSS);
}

void OctetAnalyzerPlugin::drawQuad(quadHists* qh, const std::string& subfolder, const char* opt) {
	assert(qh);
	for(AFPState afp = AFP_OFF; afp <= AFP_ON; ++afp) {
		for(GVState gv=GV_CLOSED; gv<=GV_OPEN; ++gv) {
			if(!qh->fgbg[afp]->h[gv]->Integral()) continue; // automatically skip empty histograms
			qh->fgbg[afp]->h[gv]->Draw(opt);
			printCanvas(subfolder+"/"+qh->getHistoName(AFPState(afp),gv));
		}
	}
}

void OctetAnalyzerPlugin::drawQuadSides(quadHists* qhE, quadHists* qhW, bool combineAFP, const std::string& subfolder, const std::string& opt) {
	assert(qhE && qhW);
	std::vector<TH1*> hToPlot;
	for(GVState gv=GV_CLOSED; gv<=GV_OPEN; ++gv) {
		for(AFPState afp = AFP_OFF; afp <= AFP_ON; ++afp) {
			qhE->fgbg[afp]->h[gv]->SetLineColor(2);
			qhW->fgbg[afp]->h[gv]->SetLineColor(4);
			if(!qhE->fgbg[afp]->h[gv]->Integral() && !qhW->fgbg[afp]->h[gv]->Integral()) continue;
			hToPlot.push_back(qhE->fgbg[afp]->h[gv]);
			hToPlot.push_back(qhW->fgbg[afp]->h[gv]);
			if(!combineAFP) {
				drawSimulHistos(hToPlot,opt);
				printCanvas(subfolder+"/"+qhE->name+(afp?"_On":"_Off")+(gv?"":"_BG"));
				hToPlot.clear();
			} else {
				qhE->fgbg[afp]->h[gv]->SetLineStyle(1+2*afp);
				qhW->fgbg[afp]->h[gv]->SetLineStyle(1+2*afp);
			}
		}
		if(combineAFP && hToPlot.size()) {
			drawSimulHistos(hToPlot,opt,qhE->title+(gv?"":" Background"));
			printCanvas(subfolder+"/"+qhE->name+(gv?"":"_BG"));
			hToPlot.clear();
		}
	}
}


/* --------------------------------------------------- */


unsigned int processPulsePair(OctetAnalyzer& OA, const Octet& PP) {
	unsigned int nproc = 0;
	std::vector<Octet> triads = PP.getSubdivs(3,false);
	printf("Processing pulse pair for %s containing %i triads...\n",PP.octName().c_str(),int(triads.size()));
	for(std::vector<Octet>::iterator sd = triads.begin(); sd != triads.end(); sd++) {
		nproc++;
		ProcessedDataScanner* PDSs[2] = {NULL,NULL};
		for(GVState gv=GV_CLOSED; gv<=GV_OPEN; ++gv) {
			PDSs[gv] = new PostOfficialAnalyzer(true);
			PDSs[gv]->addRuns(sd->getAsymRuns(gv));
		}
		if(PDSs[0]->getnFiles()+PDSs[1]->getnFiles())
			OA.loadProcessedData(sd->octAFPState(),*PDSs[1],*PDSs[0]);
		delete(PDSs[0]);
		delete(PDSs[1]);
	}
	return nproc;
}

unsigned int processOctets(OctetAnalyzer& OA, const std::vector<Octet>& Octs, double replaceIfOlder, bool doPlots) {
	
	unsigned int nproc = 0;
	
	unsigned int octn = 0;
	
	for(std::vector<Octet>::const_iterator octit = Octs.begin(); octit != Octs.end(); octit++) {
		octn++;
		//if(octn == 2+1 || octn == 17+1 || octn == 30+1) continue; // anomalous gamma backgrounds
		//if(octn>29) continue;
		
		// check if there are any runs to process
		printf("Processing octet '%s' at division %i...\n",octit->octName().c_str(),octit->divlevel);
		if(!octit->getNRuns()) {
			printf("\tThat was too easy (octet contained zero runs).\n");
			continue;
		}
		OA.qOut.insert("Octet",octit->toStringmap());
		
		if(octit->divlevel<=2) {
			// make sub-Analyzer for this octet, to load data if already available, otherwise re-process
			OctetAnalyzer* subOA;
			std::string inflname = OA.basePath+"/"+octit->octName()+"/"+octit->octName();
			double fAge = fileAge(inflname+".root");
			if(fileExists(inflname+".root") && fileExists(inflname+".txt") && fAge < replaceIfOlder) {
				printf("Octet '%s' already scanned %.1fh ago; skipping\n",octit->octName().c_str(),fAge/3600);
				subOA = (OctetAnalyzer*)OA.makeAnalyzer(octit->octName(),inflname);
			} else {
				subOA = (OctetAnalyzer*)OA.makeAnalyzer(octit->octName(),"");
				subOA->depth = octit->divlevel;
				nproc += processOctets(*subOA,octit->getSubdivs(octit->divlevel+1,false),replaceIfOlder, doPlots);
			}
			OA.addSegment(*subOA);
			delete(subOA);
		} else {
			nproc += processPulsePair(OA,*octit);
		}
	}
	
	// make output for octet
	if(OA.needsSubtraction)
		OA.bgSubtractAll();
	OA.calculateResults();
	OA.uploadAnaResults();
	if(doPlots)
		OA.makePlots();
	OA.write();
	OA.setWriteRoot(true);
	
	return nproc;
}


