#include "OctetAnalyzer.hh"
#include "PathUtils.hh"
#include "strutils.hh"
#include "GraphicsUtils.hh"
#include "EnergyCalibrator.hh"
#include "PostOfficialAnalyzer.hh"
#include "CalDBSQL.hh"
#include "Types.hh"

void quadHists::setFillPoint(AFPState afp, GVState gv) {
	smassert(gv==GV_CLOSED || gv==GV_OPEN);
	smassert(afp==AFP_OFF||afp==AFP_ON);
	fillPoint = fgbg[afp]->h[gv];
	smassert(fillPoint);
}

bool quadHists::isEquivalent(const quadHists& qh) const {
	if(name != qh.name) return false;
	return true;
	//return (h[AFP_OFF][1]->GetNbinsX() == qh.h[AFP_OFF][1]->GetNbinsX()
	//		&& h[AFP_OFF][1]->GetNbinsY() == qh.h[AFP_OFF][1]->GetNbinsY()
	//		&& h[AFP_OFF][1]->GetNbinsZ() == qh.h[AFP_OFF][1]->GetNbinsZ());	///< TODO more rigorous check
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
	smassert(coreHists.find(qh->getName())==coreHists.end());	// don't duplicate names!
	for(AFPState afp = AFP_OFF; afp <= AFP_ON; ++afp)
		qh->fgbg[afp] = registerFGBGPair(hname,title,nbins,xmin,xmax,afp,s);
	coreHists.insert(std::make_pair(qh->getName(),qh));
	return qh;
}

quadHists* OctetAnalyzer::registerCoreHist(const TH1& hTemplate, Side s) {
	quadHists* qh = new quadHists(hTemplate.GetName(),hTemplate.GetTitle(),s);
	smassert(coreHists.find(qh->getName())==coreHists.end());	// don't duplicate names!
	for(AFPState afp = AFP_OFF; afp <= AFP_ON; ++afp)
		qh->fgbg[afp] = registerFGBGPair(hTemplate,AFPState(afp),s);
	coreHists.insert(std::make_pair(qh->getName(),qh));
	return qh;
}

quadHists* OctetAnalyzer::cloneQuadHist(const quadHists* qh, const std::string& newName, const std::string& newTitle) {
	smassert(qh);
	quadHists* qnew = new quadHists(newName,newTitle,qh->mySide);
	for(AFPState afp = AFP_OFF; afp <= AFP_ON; ++afp)
		qnew->fgbg[afp] = cloneFGBGPair(*qh->fgbg[afp],newName,newTitle);
	coreHists.insert(std::make_pair(qnew->getName(),qnew));
	return qnew;
}

OctetAnalyzer::OctetAnalyzer(OutputManager* pnt, const std::string& nm, const std::string& inflName):
RunAccumulator(pnt,nm,inflName) { }

void OctetAnalyzer::setFillPoints(AFPState afp, GVState gv) {
	for(std::map<std::string,quadHists*>::iterator it = coreHists.begin(); it != coreHists.end(); it++)
		it->second->setFillPoint(afp,gv);
}

quadHists* OctetAnalyzer::getCoreHist(const std::string& qname) {
	std::map<std::string,quadHists*>::iterator it = coreHists.find(qname);
	smassert(it != coreHists.end());
	return it->second;
}

const quadHists* OctetAnalyzer::getCoreHist(const std::string& qname) const {
	std::map<std::string,quadHists*>::const_iterator it = coreHists.find(qname);
	smassert(it != coreHists.end());
	return it->second;
}

void OctetAnalyzer::loadSimData(Sim2PMT& simData, unsigned int nToSim, bool countAll) {
	setFillPoints(simData.getAFP(),GV_OPEN);
	RunAccumulator::loadSimData(simData,nToSim,countAll);
}

void OctetAnalyzer::loadBothAFP(Sim2PMT& simData, unsigned int nToSim) {
	printf("Loading %i events of simulated data doubled over AFP states...\n",nToSim);
	simData.resetSimCounters();
	while(!nToSim || simData.nSimmed<=nToSim) {
		bool np = simData.nextPoint();
		for(AFPState afp = AFP_OFF; afp <= AFP_ON; ++afp) {
			setCurrentState(afp,GV_OPEN);
			setFillPoints(afp,GV_OPEN);
			simData.setAFP(afp);
			simData.calcReweight();
			loadSimPoint(simData);
		}
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

TH1* OctetAnalyzer::flipperSummedRate(const quadHists* qh, GVState gv, bool doNorm) const {
	smassert(qh);
	smassert(gv == GV_OPEN || gv == GV_CLOSED);
	TH1* h = (TH1*)qh->fgbg[AFP_OFF]->h[gv]->Clone();
	h->SetTitle(qh->title.c_str());
	h->Add(qh->fgbg[AFP_ON]->h[gv]);
	Side s = qh->mySide;
	if(doNorm)
		h->Scale(1.0/h->GetXaxis()->GetBinWidth(1)/(totalTime[AFP_ON][gv][s]+totalTime[AFP_OFF][gv][s]));
	return h;
}

/* --------------------------------------------------- */

TH1* OctetAnalyzerPlugin::calculateSR(const std::string& hname, const quadHists* qEast, const quadHists* qWest, bool fg, bool instr) {
	smassert(qEast && qWest);
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

TH1* OctetAnalyzerPlugin::calculateSuperSum(const std::string& hname, const quadHists* qEast, const quadHists* qWest, GVState gv, bool toRate) {
	smassert(qEast && qWest);
	smassert(gv <= GV_OPEN);
	TH1* hR = (TH1*)qEast->fgbg[AFP_ON]->h[gv]->Clone("SuperSum_intermediate");
	hR->Multiply(qWest->fgbg[AFP_OFF]->h[gv]);
	if(toRate) hR->Scale(1.0/(myA->totalTime[AFP_ON][gv][EAST]*myA->totalTime[AFP_OFF][gv][WEST]));
	
	TH1* hSS = (TH1*)qEast->fgbg[AFP_OFF]->h[gv]->Clone(hname.c_str());
	hSS->Multiply(qWest->fgbg[AFP_ON]->h[gv]);
	if(toRate) hSS->Scale(1.0/(myA->totalTime[AFP_OFF][gv][EAST]*myA->totalTime[AFP_ON][gv][WEST]));
	
	sqrtHist(hR);
	sqrtHist(hSS);
	hSS->Add(hR);
	hSS->Scale(1.0/hSS->GetBinWidth(1));
	
	hSS->SetTitle((qEast->title+" SuperSum").c_str());
	hSS->SetLineColor(1);
	hSS->SetLineStyle(1);
	hSS->GetXaxis()->SetTitle(qEast->fgbg[AFP_ON]->h[gv]->GetXaxis()->GetTitle());
	
	delete(hR);
	return (TH1*)myA->addObject(hSS);
}

void OctetAnalyzerPlugin::drawQuad(quadHists* qh, const std::string& subfolder, const char* opt) {
	smassert(qh);
	for(AFPState afp = AFP_OFF; afp <= AFP_ON; ++afp) {
		for(GVState gv=GV_CLOSED; gv<=GV_OPEN; ++gv) {
			if(!qh->fgbg[afp]->h[gv]->Integral()) continue; // automatically skip empty histograms
			qh->fgbg[afp]->h[gv]->Draw(opt);
			printCanvas(subfolder+"/"+qh->getHistoName(AFPState(afp),gv));
		}
	}
}

void OctetAnalyzerPlugin::drawQuadSides(quadHists* qhE, quadHists* qhW, bool combineAFP, const std::string& subfolder, const std::string& opt) {
	smassert(qhE && qhW);
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

