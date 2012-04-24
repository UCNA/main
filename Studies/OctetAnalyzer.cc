#include "OctetAnalyzer.hh"
#include "PathUtils.hh"
#include "strutils.hh"
#include "GraphicsUtils.hh"
#include "EnergyCalibrator.hh"
#include "PostOfficialAnalyzer.hh"
#include "CalDBSQL.hh"
#include "Types.hh"

void quadHists::setFillPoint(AFPState afp, GVState gv) {
	if(!fillPoint) return;
	assert(gv==GV_CLOSED || gv==GV_OPEN);
	assert(afp==AFP_OFF||afp==AFP_ON);
	*fillPoint = fgbg[afp].h[gv];
}

bool quadHists::isEquivalent(const quadHists& qh) const {
	if(name != qh.name) return false;
	return true;
	//return (h[AFP_OFF][1]->GetNbinsX() == qh.h[AFP_OFF][1]->GetNbinsX()
	//		&& h[AFP_OFF][1]->GetNbinsY() == qh.h[AFP_OFF][1]->GetNbinsY()
	//		&& h[AFP_OFF][1]->GetNbinsZ() == qh.h[AFP_OFF][1]->GetNbinsZ());	//< TODO more rigorous check
}

void quadHists::operator+=(const quadHists& qh) {
	for(unsigned int afp=AFP_OFF; afp<=AFP_ON; afp++)
		fgbg[afp] += qh.fgbg[afp];
}

void quadHists::operator*=(double c) {
	for(unsigned int afp=AFP_OFF; afp<=AFP_ON; afp++)
		fgbg[afp] *= c;
}

void quadHists::setDrawMinimum(double y) {
	for(unsigned int afp = AFP_OFF; afp <= AFP_ON; afp++)
		for(unsigned int fg = 0; fg <= 1; fg++)
			fgbg[afp].h[fg]->SetMinimum(y);
}


/* --------------------------------------------------- */



quadHists OctetAnalyzer::registerCoreHist(const std::string& hname, const std::string& title,
										  unsigned int nbins, float xmin, float xmax,
										  Side s, TH1F** fillPoint) {
	quadHists qh(hname,s);
	assert(coreHists.find(qh.getName())==coreHists.end());	// don't duplicate names!
	qh.fillPoint = (TH1**)fillPoint;
	for(unsigned int afp = AFP_OFF; afp <= AFP_ON; afp++)
		qh.fgbg[afp] = registerFGBGPair(hname,title,nbins,xmin,xmax,AFPState(afp),s);
	coreHists.insert(std::make_pair(qh.getName(),qh));
	return qh;
}

quadHists OctetAnalyzer::registerCoreHist(const TH1& hTemplate, Side s, TH1** fillPoint) {
	std::string hname = hTemplate.GetName();
	quadHists qh(hname,s);
	assert(coreHists.find(qh.getName())==coreHists.end());	// don't duplicate names!
	qh.fillPoint = fillPoint;
	for(unsigned int afp = AFP_OFF; afp <= AFP_ON; afp++)
		qh.fgbg[afp] = registerFGBGPair(hTemplate,AFPState(afp),s);
	coreHists.insert(std::make_pair(qh.getName(),qh));
	return qh;
}

quadHists OctetAnalyzer::cloneQuadHist(const quadHists& qh, const std::string& newName) {
	quadHists qnew(newName,qh.mySide);
	for(unsigned int afp = AFP_OFF; afp <= AFP_ON; afp++)
		qnew.fgbg[afp] = cloneFGBGPair(qh.fgbg[afp],newName);
	return qnew;
}

OctetAnalyzer::OctetAnalyzer(OutputManager* pnt, const std::string& nm, const std::string& inflName):
RunAccumulator(pnt,nm,inflName), depth(-1) { }

void OctetAnalyzer::setFillPoints(AFPState afp, GVState gv) {
	for(std::map<std::string,quadHists>::iterator it = coreHists.begin(); it != coreHists.end(); it++)
		it->second.setFillPoint(afp,gv);
}

quadHists& OctetAnalyzer::getCoreHist(const std::string& qname) {
	std::map<std::string,quadHists>::iterator it = coreHists.find(qname);
	assert(it != coreHists.end());
	return it->second;
}

const quadHists& OctetAnalyzer::getCoreHist(const std::string& qname) const {
	std::map<std::string,quadHists>::const_iterator it = coreHists.find(qname);
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

void OctetAnalyzer::loadSimData(Sim2PMT& simData, unsigned int nToSim) {
	for(AFPState afp=AFP_OFF; afp<=AFP_ON; ++afp) {
		if(afp != simData.getAFP() && simData.getAFP() != AFP_OTHER) continue;
		setFillPoints(afp,GV_OPEN);
		RunAccumulator::loadSimData(simData,nToSim);
	}
}

TH1* OctetAnalyzer::calculateSR(const std::string& hname, const quadHists& qEast, const quadHists& qWest, bool fg) {
	// first calculate R = E+W-/E-W+
	TH1* hR = (TH1*)qEast.fgbg[AFP_ON].h[fg]->Clone("SR_intermediate_R");
	hR->Multiply(qWest.fgbg[AFP_OFF].h[fg]);
	hR->Scale(1.0/(totalTime[AFP_ON][fg].t[EAST]*totalTime[AFP_OFF][fg].t[WEST]));
	
	TH1* hAsym = (TH1*)qEast.fgbg[AFP_OFF].h[fg]->Clone(hname.c_str());
	hAsym->Multiply(qWest.fgbg[AFP_ON].h[fg]);
	hAsym->Scale(1.0/(totalTime[AFP_OFF][fg].t[EAST]*totalTime[AFP_ON][fg].t[WEST]));
	
	hR->Divide(hAsym);

	// super-ratio asymmetry (1-sqrt(R))/(1+sqrt(R))
	for(unsigned int i=0; i<totalBins(hR); i++) {
		double r = hR->GetBinContent(i);
		double dr = hR->GetBinError(i);
		if(r>0 && r==r) {
			hAsym->SetBinContent(i,(1.0-sqrt(r))/(1.0+sqrt(r)));
			hAsym->SetBinError( i, dr/(sqrt(r)*(1+sqrt(r))*(1+sqrt(r))) );
		} else {
			hAsym->SetBinContent(i,0);
			hAsym->SetBinError(i,1e6);
		}
	}
	
	hAsym->SetLineColor(1);
	hAsym->SetLineStyle(1);
	hAsym->SetTitle((std::string(hAsym->GetTitle())+" Asymmetry").c_str());
	
	delete(hR);
	return (TH1*)addObject(hAsym);
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

TH1* OctetAnalyzer::calculateSuperSum(const std::string& hname, const quadHists& qEast, const quadHists& qWest, bool fg) {
	TH1* hR = (TH1*)qEast.fgbg[AFP_ON].h[fg]->Clone("SuperSum_intermediate");
	hR->Multiply(qWest.fgbg[AFP_OFF].h[fg]);
	hR->Scale(1.0/(totalTime[AFP_ON][fg].t[EAST]*totalTime[AFP_OFF][fg].t[WEST]));
	
	TH1* hSS = (TH1*)qEast.fgbg[AFP_OFF].h[fg]->Clone(hname.c_str());
	hSS->Multiply(qWest.fgbg[AFP_ON].h[fg]);
	hSS->Scale(1.0/(totalTime[AFP_OFF][fg].t[EAST]*totalTime[AFP_ON][fg].t[WEST]));
	
	sqrtHist(hR);
	sqrtHist(hSS);
	hSS->Add(hR);
	hSS->Scale(0.5);
	
	hSS->SetTitle((std::string(hSS->GetTitle())+" SuperSum").c_str());
	hSS->SetLineColor(1);
	hSS->SetLineStyle(1);
	
	delete(hR);
	return (TH1*)addObject(hSS);
}

void OctetAnalyzer::drawQuad(quadHists& qh, const std::string& subfolder, const char* opt) {
	defaultCanvas->cd();
	for(unsigned int afp = AFP_OFF; afp <= AFP_ON; afp++) {
		for(unsigned int fg = 0; fg <= 1; fg++) {
			if(!qh.fgbg[afp].h[fg]->Integral()) continue; // automatically skip empty histograms
			qh.fgbg[afp].h[fg]->Draw(opt);
			printCanvas(subfolder+"/"+qh.getHistoName(AFPState(afp),fg));
		}
	}
}

void OctetAnalyzer::drawQuadSides(quadHists& qhE, quadHists& qhW, bool combineAFP, const std::string& subfolder, const std::string& opt) {
	defaultCanvas->cd();
	std::vector<TH1*> hToPlot;
	for(unsigned int fg = 0; fg <= 1; fg++) {
		for(unsigned int afp = AFP_OFF; afp <= AFP_ON; afp++) {
			qhE.fgbg[afp].h[fg]->SetLineColor(2);
			qhW.fgbg[afp].h[fg]->SetLineColor(4);
			hToPlot.push_back(qhE.fgbg[afp].h[fg]);
			hToPlot.push_back(qhW.fgbg[afp].h[fg]);
			if(!combineAFP) {
				drawSimulHistos(hToPlot,opt);
				printCanvas(subfolder+"/"+qhE.name+(afp?"_On":"_Off")+(fg?"":"_BG"));
				hToPlot.clear();
			} else {
				qhE.fgbg[afp].h[fg]->SetLineStyle(1+2*afp);
				qhW.fgbg[afp].h[fg]->SetLineStyle(1+2*afp);
			}
		}
		if(combineAFP) {
			drawSimulHistos(hToPlot,opt);
			printCanvas(subfolder+"/"+qhE.name+(fg?"":"_BG"));
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
		for(unsigned int fg = 0; fg <= 1; fg++) {
			PDSs[fg] = new PostOfficialAnalyzer(true);
			PDSs[fg]->addRuns(sd->getAsymRuns(fg));
		}
		if(PDSs[0]->getnFiles()+PDSs[1]->getnFiles())
			OA.loadProcessedData(sd->octAFPState(),*PDSs[1],*PDSs[0]);
		delete(PDSs[0]);
		delete(PDSs[1]);
	}
	return nproc;
}

unsigned int processOctets(OctetAnalyzer& OA, const std::vector<Octet>& Octs, double replaceIfOlder) {
	
	unsigned int nproc = 0;
	
	for(std::vector<Octet>::const_iterator octit = Octs.begin(); octit != Octs.end(); octit++) {
		
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
				nproc += processOctets(*subOA,octit->getSubdivs(octit->divlevel+1,false),replaceIfOlder);
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
	OA.makePlots();
	OA.write();
	OA.setWriteRoot(true);
	
	return nproc;
}

void simForRun(OctetAnalyzer& OA, Sim2PMT& simData, RunNum rn, AFPState afp, unsigned int nToSim, double rntime) {
	assert(afp == AFP_ON || afp == AFP_OFF || afp == AFP_OTHER);
	PMTCalibrator PCal(rn,CalDBSQL::getCDB());
	simData.setCalibrator(PCal);
	simData.setAFP(afp);
	OA.loadSimData(simData,nToSim);
	OA.runTimes.add(rn,rntime);
	OA.totalTime[afp][GV_OPEN] += rntime;
}

unsigned int simuClone(const std::string& basedata, OctetAnalyzer& OA, Sim2PMT& simData, double simfactor, double replaceIfOlder) {
	
	printf("\n------ Cloning asymmetry data in '%s'\n------                        to '%s'...\n",basedata.c_str(),OA.basePath.c_str());
	int nCloned = 0;
	
	// check if simulation is already up-to-date
	if(OA.getInflAge() && OA.getInflAge() < replaceIfOlder) {
		printf("\tSimulations in '%s' already recently generated, update skipped...\n",OA.basePath.c_str());
		return nCloned;
	}
	
	// clear any existing data from out-of-date input
	OA.zeroSavedHists();
	OA.zeroCounters();
	
	// load original data for comparison
	std::vector<std::string> datpath = split(strip(basedata,"/"),"/");
	assert(datpath.size()>0);
	OctetAnalyzer* origOA = (OctetAnalyzer*)OA.makeAnalyzer("nameUnused",basedata+"/"+datpath.back());
	
	// load/clone data in all subdirectories, if they exist
	int nClonable = 0;
	std::vector<std::string> fnames = listdir(basedata);
	for(std::vector<std::string>::iterator it = fnames.begin(); it != fnames.end(); it++) {
		// check whether data directory contains cloneable subdirectories
		std::string datinfl = basedata+"/"+(*it)+"/"+(*it);
		if(!OctetAnalyzer::inflExists(datinfl)) continue;
		std::string siminfl = OA.basePath+"/"+(*it)+"/"+(*it);		
		nClonable++;
		
		// load cloned sub-data
		OctetAnalyzer* subOA = (OctetAnalyzer*)OA.makeAnalyzer(*it,OctetAnalyzer::inflExists(siminfl)?siminfl:"");
		subOA->depth = OA.depth+1;
		nCloned += simuClone(basedata+"/"+(*it),*subOA,simData,simfactor,replaceIfOlder);
		OA.addSegment(*subOA);
		delete(subOA);
	}
	
	// if there were no clonable subdrectories, we need to actually simulate the data
	if(!nClonable) {
		printf("\tNo data subdirectories found; assume data here needs cloning...\n");
		printf("\tProcessing simulation data for each pulse-pair run...\n");
		for(std::map<RunNum,double>::iterator it = origOA->runCounts.counts.begin(); it != origOA->runCounts.counts.end(); it++) {
			nCloned++;
			if(!it->first || !it->second) continue;
			RunInfo RI = CalDBSQL::getCDB()->getRunInfo(it->first);
			if(RI.gvState != GV_OPEN) continue;	// no simulation for background runs
			// estimate background count share for this run (and reduce simulation by this amount)
			double bgEst = origOA->getTotalCounts(RI.afpState,GV_CLOSED)*origOA->getRunTime(it->first)/origOA->getTotalTime(RI.afpState,0).t[BOTH];
			if(it->second <= bgEst) continue;
			int nToSim = (int)it->second-bgEst;
			printf("\n\t---Simulation cloning for run %i (%i+%i counts)---\n",it->first,nToSim,(int)bgEst);
			simForRun(OA, simData, it->first,  RI.afpState, (unsigned int)nToSim,origOA->getRunTime(it->first));
		}
		// and clone background counts in original data
		OA.simBgFlucts(*origOA,simfactor);
	}
	
	// generate output
	OA.calculateResults();
	OA.makePlots();
	origOA->calculateResults();
	OA.compareMCtoData(*origOA,simfactor);
	OA.write();
	OA.setWriteRoot(true);
	
	delete(origOA);
	return nCloned;
}
