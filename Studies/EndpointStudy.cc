#include "EndpointStudy.hh"
#include "CalDBSQL.hh"
#include "MultiGaus.hh"
#include "LinHistCombo.hh"
#include "G4toPMT.hh"
#include "KurieFitter.hh"
#include "PostOfficialAnalyzer.hh"
#include "GraphicsUtils.hh"
#include "GraphUtils.hh"
#include <TStyle.h>
#include <time.h>
#include <TH3F.h>
#include <time.h>
#include <stdlib.h>
#include <utility>
#include <algorithm>

/// convert cathode segment to Stringmap
Stringmap cathseg2sm(const CathodeSeg& c) {
	Stringmap m;
	m.insert("side",sideSubst("%c",c.s));
	m.insert("plane",c.d==X_DIRECTION?"x":"y");
	m.insert("i",c.i);
	m.insert("max",c.max.x);
	m.insert("d_max",c.max.err);
	m.insert("height",c.height.x);
	m.insert("d_height",c.height.err);
	m.insert("width",c.width.x);
	m.insert("d_width",c.width.err);
	m.insert("center",c.center.x);
	m.insert("d_center",c.center.x);
	m.insert("position",c.pos);
	m.insert("fill_frac",c.fill_frac);
	return m;
}

CathodeSeg sm2cathseg(const Stringmap& m) {
	CathodeSeg c;
	c.s=(m.getDefault("side","E")=="E")?EAST:WEST;
	c.d=m.getDefault("plane","x")=="x"?X_DIRECTION:Y_DIRECTION;
	c.i=m.getDefault("i",0);
	c.pos=m.getDefault("position",0);
	c.height=float_err(m.getDefault("height",0),m.getDefault("d_height",0));
	c.width=float_err(m.getDefault("width",0),m.getDefault("d_width",0));
	c.center=float_err(m.getDefault("center",0),m.getDefault("d_center",0));
	return c;
}

/// convert Stringmap to SectorDat
SectorDat sm2sd(const Stringmap& m) {
	SectorDat sd;
	sd.s = m.getDefault("side","")=="E"?EAST:WEST;
	sd.t = (unsigned int)(m.getDefault("tube",0));
	sd.m = (unsigned int)(m.getDefault("m",0));
	sd.eta = m.getDefault("eta",0);
	sd.low_peak = float_err(m.getDefault("xe_lo",0),m.getDefault("d_xe_lo",0));
	sd.low_peak_width = float_err(m.getDefault("xe_lo_w",0),m.getDefault("d_xe_lo_w",0));
	sd.xe_ep = float_err(m.getDefault("xe_hi",0),m.getDefault("d_xe_hi",0));
	return sd;
}

/// convert SectorDat to Stringmap
Stringmap sd2sm(const SectorDat& sd) {
	Stringmap m;
	m.insert("side",sd.s==EAST?"E":"W");
	m.insert("tube",sd.t);
	m.insert("m",sd.m);
	m.insert("eta",sd.eta);
	m.insert("xe_lo",sd.low_peak.x);
	m.insert("d_xe_lo",sd.low_peak.err);
	m.insert("xe_lo_w",sd.low_peak_width.x);
	m.insert("d_xe_lo_w",sd.low_peak_width.err);
	m.insert("xe_hi",sd.xe_ep.x);
	m.insert("d_xe_hi",sd.xe_ep.err);
	return m;
}


//-----------------------------------------------


WirechamberAnalyzer::WirechamberAnalyzer(OutputManager* pnt, const std::string& nm, const std::string& infl):
RunAccumulator(pnt,nm,infl) {
	TH2F hPositionsTemplate("hPositions","Event Positions",200,-60,60,200,-60,60);
	hPositionsTemplate.GetXaxis()->SetTitle("x position [mm]");
	hPositionsTemplate.GetYaxis()->SetTitle("y position [mm]");
	TH2F hPositionsRawTemplate("hPositionsRaw","Event Raw Positions",200,-60,60,200,-60,60);
	hPositionsRawTemplate.GetXaxis()->SetTitle("x position [mm]");
	hPositionsRawTemplate.GetYaxis()->SetTitle("y position [mm]");
	for(Side s = EAST; s <= WEST; ++s) {
		hitPos[s] = registerFGBGPair(hPositionsTemplate,AFP_OTHER,s);
		hitPosRaw[s] = registerFGBGPair(hPositionsRawTemplate,AFP_OTHER,s);
		for(AxisDirection d = X_DIRECTION; d <= Y_DIRECTION; ++d) {
			for(unsigned int c=0; c < kMaxCathodes; c++) {
				TH2F hCathHitposTemplate((std::string("hCathHitpos_")+(d==X_DIRECTION?"x_":"y_")+itos(c)).c_str(),
										 "Cathode Hit Positions",64,-0.5,0.5,1000,0,1000);
				TH2F hCathNormTemplate((std::string("hCathNorm_")+(d==X_DIRECTION?"x_":"y_")+itos(c)).c_str(),
									   "Normalized Cathode",51,-0.5,0.5,100,0,5);
				
				cathHitpos[s][d][c] = registerFGBGPair(hCathHitposTemplate,AFP_OTHER,s);
				cathHitpos[s][d][c]->setAxisTitle(X_DIRECTION,"Normalized Raw Position");
				cathHitpos[s][d][c]->setAxisTitle(Y_DIRECTION,"Energy [keV]");
				cathNorm[s][d][c] = registerFGBGPair(hCathNormTemplate,AFP_OTHER,s);
				cathNorm[s][d][c]->setAxisTitle(X_DIRECTION,"Normalized Position");
				cathNorm[s][d][c]->setAxisTitle(Y_DIRECTION,"Normalized Cathode");
			}
		}
	}
}

void WirechamberAnalyzer::fillCoreHists(ProcessedDataScanner& PDS, double weight) {
	const Side s = PDS.fSide;
	if(!(PDS.fType == TYPE_0_EVENT && PDS.fPID == PID_BETA && (s==EAST||s==WEST))) return;
	for(AxisDirection d = X_DIRECTION; d <= Y_DIRECTION; ++d) {
		unsigned int n;
		float c;
		PDS.ActiveCal->toLocal(s,d,PDS.wires[s][d].rawCenter,n,c);	
		assert(n<kMaxCathodes);
		((TH2F*)cathHitpos[s][d][n]->h[currentGV])->Fill(c,PDS.scints[s].energy.x,weight);
		PDS.ActiveCal->toLocal(s,d,PDS.wires[s][d].center,n,c);	
		if(PDS.mwpcs[s].anode > 1000) continue;	// avoid most massively clipped events
		assert(n<kMaxCathodes);
		((TH2F*)cathNorm[s][d][n]->h[currentGV])->Fill(c,PDS.cathodes[s][d][n]/PDS.mwpcs[s].anode,weight);
	}
	((TH2F*)(hitPos[s]->h[currentGV]))->Fill(PDS.wires[s][X_DIRECTION].center,PDS.wires[s][Y_DIRECTION].center,weight);
	((TH2F*)(hitPosRaw[s]->h[currentGV]))->Fill(PDS.wires[s][X_DIRECTION].rawCenter,PDS.wires[s][Y_DIRECTION].rawCenter,weight);
}

void WirechamberAnalyzer::makePlots() {
	for(Side s = EAST; s <= WEST; ++s) {
		TH2F* hPos = (TH2F*)hitPos[s]->h[GV_OPEN]->Clone();
		hPos->Scale(10*M_PI*52*52/(hPos->GetXaxis()->GetBinWidth(1)*hPos->GetYaxis()->GetBinWidth(1)*hPos->Integral()));
		hPos->SetMaximum(14);
		hPos->Draw("COL Z");
		printCanvas(sideSubst("hPos_%c",s));
		delete hPos;
		
		TH2F* hRaw = (TH2F*)hitPosRaw[s]->h[GV_OPEN]->Clone();
		hRaw->Scale(10*M_PI*52*52/(hRaw->GetXaxis()->GetBinWidth(1)*hRaw->GetYaxis()->GetBinWidth(1)*hRaw->Integral()));
		hRaw->SetMaximum(14);
		hRaw->Draw("COL Z");
		printCanvas(sideSubst("hPosRaw_%c",s));
		delete hRaw;
	}
}


//-----------------------------------------------

double PositionBinner::fidRadius = 50.;

PositionBinner::PositionBinner(OutputManager* pnt, const std::string& nm, unsigned int nr, const std::string& infl):
WirechamberAnalyzer(pnt,nm,infl), sects(nr,PositionBinner::fidRadius), sectorPlots(false) {
	
	// load sector cutter
	if(fIn) {
		QFile qOld(inflname+".txt");
		Stringmap sct = qOld.getFirst("SectorCutter");
		sects = SectorCutter(int(sct.getDefault("nRings",0)),sct.getDefault("radius",0));
		assert(sects.n && sects.r);
	}
	
	// save sector cutter
	Stringmap ms;
	ms.insert("nRings",sects.n);
	ms.insert("radius",sects.r);
	ms.insert("nSectors",sects.nSectors());
	qOut.insert("SectorCutter",ms);
	
	// set up histograms, data
	energySpectrum = registerFGBGPair("hEnergy","Combined Energy",200,-100,1200,AFP_OTHER);
	energySpectrum->h[GV_OPEN]->SetLineColor(2);
	for(Side s = EAST; s <= WEST; ++s) {
		for(unsigned int m=0; m<sects.nSectors(); m++) {
			sectAnode[s].push_back(registerFGBGPair("hAnode_Sector_"+itos(m),"Sector "+itos(m)+" Anode",200,0,4000,AFP_OTHER,s));
			for(unsigned int t=0; t<nBetaTubes; t++) {
				sectEnergy[s][t].push_back(registerFGBGPair("hTuben_"+itos(t)+"_Sector_"+itos(m),
															"Sector "+itos(m)+" Tube "+itos(t)+" Energy",200,-100,2000,AFP_OTHER,s));
				SectorDat sd;
				sd.s = s;
				sd.t = t;
				sd.m = m;
				sectDat[s][t].push_back(sd);
			}
		}
	}
	
	// load sector data
	if(fIn) {
		QFile qOld(inflname+".txt");
		std::vector<Stringmap> sds = qOld.retrieve("sectDat");
		for(std::vector<Stringmap>::iterator it = sds.begin(); it != sds.end(); it++) {
			SectorDat sd = sm2sd(*it);
			sectDat[sd.s][sd.t][sd.m] = sd;
		}
	}
	
}

void PositionBinner::fillCoreHists(ProcessedDataScanner& PDS, double weight) {
	WirechamberAnalyzer::fillCoreHists(PDS,weight);
	const Side s = PDS.fSide;
	if(!(PDS.fType == TYPE_0_EVENT && PDS.fPID == PID_BETA && (s==EAST||s==WEST))) return;
	unsigned int m = sects.sector(PDS.wires[s][X_DIRECTION].center,PDS.wires[s][Y_DIRECTION].center);
	if(m>=sects.nSectors()) return;
	for(unsigned int t=0; t<nBetaTubes; t++)
		sectEnergy[s][t][m]->h[currentGV]->Fill(PDS.scints[s].tuben[t].x,weight);
	if(PDS.radius2(s) <= 45*45)
		energySpectrum->h[currentGV]->Fill(PDS.getEtrue(),weight);
	if(PDS.getEtrue() > 225)
		sectAnode[s][m]->h[currentGV]->Fill(PDS.mwpcs[s].anode,weight);
}

void PositionBinner::fitSpectrum(TH1* hSpec,SectorDat& sd) {
	
	hSpec->SetLineColor(2+sd.t);
	
	//----------------------
	// Low peak fit
	//----------------------
	double epGuess = 915.;
	TF1 gausFit("gasufit","gaus",0,500);
	gausFit.SetLineColor(2+sd.t);
	if(!iterGaus(hSpec,&gausFit,3,hSpec->GetBinCenter(hSpec->GetMaximumBin()),100,1.0)) {
		sd.low_peak = float_err(gausFit.GetParameter(1),gausFit.GetParError(1));
		sd.low_peak_width = float_err(gausFit.GetParameter(2),gausFit.GetParError(2));
		epGuess = 6.6*sd.low_peak.x;
	} else {
		sd.low_peak = sd.low_peak_width = 0;
	}
	
	//----------------------
	// 915keV endpoint fit
	//----------------------
	sd.xe_ep = kurieIterator(hSpec,epGuess,NULL,915.,450,750);
}

void PositionBinner::fitSectors() {
	assert(runCounts.counts.size());
	PMTCalibrator PCal(runCounts.counts.begin()->first);
	printf("\n\n---- Using Calibrator: ----\n");
	PCal.printSummary();
	for(Side s = EAST; s <= WEST; ++s) {
		for(unsigned int t=0; t<nBetaTubes; t++) {
			for(unsigned int m=0; m<sects.nSectors(); m++) {
				float x,y;
				sects.sectorCenter(m,x,y);
				sectDat[s][t][m].eta = PCal.eta(s,t,x,y);
				TH1* hSpec = sectEnergy[s][t][m]->h[GV_OPEN];
				fitSpectrum(hSpec,sectDat[s][t][m]);
				qOut.insert("sectDat",sd2sm(sectDat[s][t][m]));
			}
		}
	}	
}

void PositionBinner::calculateResults() {
	WirechamberAnalyzer::calculateResults();
	for(Side s = EAST; s <= WEST; ++s) {
		for(unsigned int t=0; t<nBetaTubes; t++) {
			hTuben[s][t] = (TH1F*)addObject(sectEnergy[s][t][0]->h[GV_OPEN]->Clone((sideSubst("hTuben_%c",s)+itos(t)).c_str()));
			hTuben[s][t]->Reset();
			for(unsigned int m=0; m<sects.nSectors(); m++) {
				if(sects.sectorCenterRadius(m) < 45.)
					hTuben[s][t]->Add(sectEnergy[s][t][m]->h[GV_OPEN]);
			}
			SectorDat sd;
			sd.s = s;
			sd.t = t;
			sd.m = 0;
			sd.eta = 1.;
			fitSpectrum(hTuben[s][t],sd);
			qOut.insert("tuben",sd2sm(sd));
		}
	}
	
}

void PositionBinner::compareMCtoData(RunAccumulator& OAdata) {
	WirechamberAnalyzer::compareMCtoData(OAdata);
	PositionBinner* PB = (PositionBinner*)&OAdata;
	defaultCanvas->cd();
	
	// overall energy spectrum
	int b0 = PB->energySpectrum->h[GV_OPEN]->FindBin(400);
	int b1 = PB->energySpectrum->h[GV_OPEN]->FindBin(800);
	energySpectrum->h[GV_OPEN]->Scale(PB->energySpectrum->h[GV_OPEN]->Integral(b0,b1)/energySpectrum->h[GV_OPEN]->Integral(b0,b1));
	drawHistoPair(PB->energySpectrum->h[GV_OPEN],energySpectrum->h[GV_OPEN]);
	printCanvas("Comparison/hEnergy");		
	
	// upload posmap
	std::string pmapname = itos(PB->runCounts.counts.begin()->first)+"-"+itos(PB->runCounts.counts.rbegin()->first)+"/"+itos(time(NULL));
	CalDBSQL* CDBout = CalDBSQL::getCDB(false);
	unsigned int pmid_ep = CDBout->newPosmap("Xe Endpoint "+pmapname,sects.n,sects.r);
	unsigned int pmid_lp = CDBout->newPosmap("Xe Low Peak "+pmapname,sects.n,sects.r);
	float x,y;
	for(unsigned int m=0; m<sects.nSectors(); m++) {
		sects.sectorCenter(m,x,y);
		for(Side s=EAST; s<=WEST; ++s) {
			for(unsigned int t=0; t<nBetaTubes; t++) {
				CDBout->addPosmapPoint(pmid_ep,s,t,m,
									   PB->sectDat[s][t][m].xe_ep.x*PB->sectDat[s][t][m].eta,
									   sectDat[s][t][m].xe_ep.x,x,y);
				CDBout->addPosmapPoint(pmid_lp,s,t,m,
									   PB->sectDat[s][t][m].low_peak.x*PB->sectDat[s][t][m].eta,
									   sectDat[s][t][m].low_peak.x,x,y);
			}
		}
	}
	Stringmap pmsm;
	pmsm.insert("pmid_ep",pmid_ep);
	pmsm.insert("pmid_lp",pmid_lp);
	pmsm.insert("name",pmapname);
	qOut.insert("posmap",pmsm);
}

void PositionBinner::makePlots() {
	WirechamberAnalyzer::makePlots();
	std::vector<TH1*> hToPlot;
	
	// overall energy spectrum
	defaultCanvas->cd();
	energySpectrum->h[GV_OPEN]->Draw();
	printCanvas("hEnergy");
	
	for(Side s = EAST; s <= WEST; ++s) {
		// tube energy
		hToPlot.clear();
		for(unsigned int t=0; t<nBetaTubes; t++)
			hToPlot.push_back(hTuben[s][t]);
		drawSimulHistos(hToPlot);
		printCanvas(sideSubst("hTuben_%c",s));
		
		// energy in each sector
		if(sectorPlots) {
			for(unsigned int m=0; m<sects.nSectors(); m++) {
				hToPlot.clear();
				for(unsigned int t=0; t<nBetaTubes; t++)
					hToPlot.push_back(sectEnergy[s][t][m]->h[GV_OPEN]);
				drawSimulHistos(hToPlot);
				for(unsigned int t=0; t<nBetaTubes; t++)
					drawVLine(sectDat[s][t][m].xe_ep.x, defaultCanvas, 2+t);
				printCanvas(sideSubst("SectorEnergy/h_%c_",s)+itos(m));
			}
		}
	}
}


//-----------------------------------------------


void process_xenon(RunNum r0, RunNum r1, unsigned int nrings) {	
		
	// scan data from each run
	std::vector<std::string> snames;
	OutputManager OM1("NameUnused",getEnvSafe("UCNA_ANA_PLOTS")+"/PositionMaps/SingleRuns/");
	std::vector<RunNum> rlist = CalDBSQL::getCDB()->findRuns("run_type = 'Xenon'",r0,r1);
	for(std::vector<RunNum>::iterator rit = rlist.begin(); rit != rlist.end(); rit++) {
		std::string singleName = "Xenon_"+itos(*rit)+"_"+itos(nrings)+"_"+dtos(PositionBinner::fidRadius);
		std::string prevFile = OM1.basePath+"/"+singleName+"/"+singleName;
		snames.push_back(singleName);
		if(r0==r1 || !fileExists(prevFile+".root")) {
			PositionBinner PB1(&OM1, singleName, nrings);
			PostOfficialAnalyzer POA(true);
			POA.redoPositions = true;
			POA.addRun(*rit);
			PB1.loadProcessedData(AFP_OTHER, GV_OPEN, POA);
			PB1.calculateResults();
			POA.writeCalInfo(PB1.qOut,"runcal");
			PB1.write();
			PB1.setWriteRoot(true);
		}
	}
	
	if(r0==r1) return;
	
	// reload data
	OutputManager OM("NameUnused",getEnvSafe("UCNA_ANA_PLOTS")+"/PositionMaps/");
	PositionBinner PB(&OM, "Xenon_"+itos(r0)+"-"+itos(r1)+"_"+itos(nrings), nrings);	
	for(std::vector<std::string>::iterator it = snames.begin(); it != snames.end(); it++) {
		std::string prevFile = OM1.basePath+"/"+*it+"/"+*it;
		PositionBinner PB1(&OM1, *it, 0, prevFile);
		PB.addSegment(PB1);
	}
	
	// finish and output
	PB.calculateResults();
	PB.fitSectors();
	PB.makePlots();
	PB.write();
	PB.setWriteRoot(true);	
}

std::vector<unsigned int> randomPermutation(unsigned int n) {
	std::vector< std::pair<double,unsigned int> > v;
	for(unsigned int i=0; i<n; i++)
		v.push_back(std::make_pair(rand(),i));
	std::sort(v.begin(),v.end());
	std::vector<unsigned int> p;
	for(unsigned int i=0; i<n; i++)
		p.push_back(v[i].second);
	return p;
}

std::string simulate_one_xenon(RunNum r, OutputManager& OM1, PositionBinner& PB, float simFactor, bool forceResim=false) {
	std::string singleName = "Xenon_"+itos(r)+"_"+itos(PB.sects.n)+"_"+dtos(PB.sects.r);
	std::string prevFile = OM1.basePath+"/"+singleName+"/"+singleName;
	if(forceResim || !fileExists(prevFile+".root")) {
		PMTCalibrator PCal(r);
		PositionBinner PBM(&OM1,singleName,PB.sects.n);
		PBM.totalTime[AFP_OTHER][GV_OPEN] += PB.totalTime[AFP_OTHER][GV_OPEN];
		PBM.runTimes += PB.runTimes;
		
		unsigned int nToSim=simFactor*PB.runCounts[r];
		printf("Data counts West: %f; to sim = %i\n",PB.energySpectrum->h[GV_OPEN]->Integral(),nToSim);
		
		// simulate for each isotope
		std::vector<PositionBinner*> PBMi;
		std::vector<std::string> isots;
		LinHistCombo LHC;
		isots.push_back("Xe125_1-2+");
		isots.push_back("Xe129_11-2-");
		isots.push_back("Xe131_11-2-");
		isots.push_back("Xe133_3-2+");
		isots.push_back("Xe133_11-2-");
		isots.push_back("Xe135_3-2+");
		if((14264 <= r && r <= 14273) || (15991 <= r && r <= 16010))
			isots.push_back("Xe135_11-2-");
		int b1 = PB.energySpectrum->h[GV_OPEN]->FindBin(1075);
		int b2 = PB.energySpectrum->h[GV_OPEN]->FindBin(1175);
		if(PB.energySpectrum->h[GV_OPEN]->Integral(b1,b2) > 100)
			isots.push_back("Xe137_7-2-");
		
		srand(time(NULL));
		std::vector<unsigned int> p = randomPermutation(isots.size());
		
		for(unsigned int n=0; n<isots.size(); n++) {
			printf("Simulating for component %s...\n",isots[p[n]].c_str());
			PBMi.push_back(new PositionBinner(&OM1,singleName+"_"+isots[p[n]],PB.sects.n));
			//G4SegmentMultiplier GSM(PB.sects);
			G4SegmentMultiplier GSM(SectorCutter(4,52.));
			GSM.setCalibrator(PCal);
			std::string simFile = "/home/mmendenhall/geant4/output/WideKev_"+isots[p[n]]+"/analyzed_";
			unsigned int nTot = 18;
			unsigned int stride = 7;
			for(unsigned int i=0; i<stride; i++)
				GSM.addFile(simFile+itos((stride*r+i)%nTot)+".root");
			PBMi.back()->loadSimData(GSM, nToSim*(isots[p[n]]=="Xe135_3-2+"?1.5:0.25));
			LHC.addTerm(PBMi.back()->energySpectrum->h[GV_OPEN]);
			printf("Done.\n");
		}
		
		// determine spectrum composition and accumulate segments
		LHC.Fit(PB.energySpectrum->h[GV_OPEN],50,1000);
		std::vector<double> counts;
		for(unsigned int i=0; i<LHC.coeffs.size(); i++) {
			PBMi[i]->scaleData(LHC.coeffs[i]);
			counts.push_back(PBMi[i]->energySpectrum->h[GV_OPEN]->Integral());
			PBM.addSegment(*PBMi[i]);
			delete(PBMi[i]);
		}
		Stringmap m;
		m.insert("nTerms",LHC.coeffs.size());
		m.insert("terms",vtos(LHC.coeffs));
		m.insert("errs",vtos(LHC.dcoeffs));
		m.insert("isots",join(isots,","));
		m.insert("counts",vtos(counts));
		m.display();
		PBM.qOut.insert("spectrumComp",m);
		PBM.qOut.insert("runcal",PCal.calSummary());
		PBM.calculateResults();
		PBM.write();
		PBM.setWriteRoot(true);
	}
	return singleName;
}

void xenon_posmap(RunNum r0, RunNum r1, unsigned int nRings) {
	std::string basePath = getEnvSafe("UCNA_ANA_PLOTS")+"/PositionMaps/";
	OutputManager OM("NameUnused",basePath);
	std::string readname = itos(r0)+"-"+itos(r1)+"_"+itos(nRings);
	
	// read in data
	PositionBinner PBdat(&OM, "Xenon_"+readname, 0,
						 basePath+"/Xenon_"+readname+"/Xenon_"+readname);
	
	// read in simulation
	PositionBinner PBsim(&OM, "SimXe_"+readname, 0,
						 basePath+"/SimXe_"+readname+"/SimXe_"+readname);
	
	// data comparison / posmap generation
	PBsim.compareMCtoData(PBdat);
}

void simulate_xenon(RunNum r0, RunNum r1, RunNum rsingle, unsigned int nRings) {
	
	std::string basePath = getEnvSafe("UCNA_ANA_PLOTS")+"/PositionMaps/";
	OutputManager OM("NameUnused",basePath);
	std::string readname = itos(r0)+"-"+itos(r1)+"_"+itos(nRings);
	printf("Simulating for position map %i-%i\n",r0,r1);
	
	// MC data for each run
	OutputManager OM1("NameUnused",basePath+"/SingleRunsSim/");
	float simFactor = 4.0;
	if(rsingle) {
		if(!CalDBSQL::getCDB()->findRuns("run_type = 'Xenon'",rsingle,rsingle).size()) return;
		std::string singleName = "Xenon_"+itos(rsingle)+"_"+itos(nRings)+"_"+dtos(PositionBinner::fidRadius);
		PositionBinner PB1(&OM, singleName, 0, basePath+"/SingleRuns/"+singleName+"/"+singleName);
		simulate_one_xenon(rsingle,OM1,PB1,simFactor,true);
		return;
	}
	
	// read in comparison data
	PositionBinner PB(&OM, "Xenon_"+readname, 0, basePath+"/Xenon_"+readname+"/Xenon_"+readname);
	std::vector<std::string> snames;
	for(std::map<RunNum,double>::const_iterator rit = PB.runCounts.counts.begin(); rit != PB.runCounts.counts.end(); rit++)
		snames.push_back(simulate_one_xenon(rit->first,OM1,PB,simFactor));
	
	// merge simulated data
	PositionBinner PBM(&OM, "SimXe_"+readname, PB.sects.n);
	PBM.isSimulated = true;
	for(std::vector<std::string>::iterator it = snames.begin(); it != snames.end(); it++) {
		std::string prevFile = OM1.basePath+"/"+*it+"/"+*it;
		PositionBinner PBM1(&OM1, *it, 0, prevFile);
		PBM.addSegment(PBM1);
	}
	
	// finish and output
	PBM.calculateResults();
	PBM.fitSectors();
	PBM.makePlots();
	PBM.compareMCtoData(PB);
	PBM.write();
	PBM.setWriteRoot(true);
}


//-----------------------------------------------


void processWirechamberCal(WirechamberAnalyzer& WCdat, WirechamberAnalyzer& WCsim) {
	RunNum r1 = WCdat.runCounts.counts.begin()->first;
	RunNum r2 = WCdat.runCounts.counts.rbegin()->first;
	OutputManager OM("MWPCCal",getEnvSafe("UCNA_ANA_PLOTS")+"/WirechamberCal/"+itos(r1)+"-"+itos(r2)+"/");
	OM.qOut.transfer(WCdat.qOut,"runcal");
	TF1 fCathCenter("fCathCenter","gaus",-0.5,0.5);
	fCathCenter.SetLineColor(2);
	
	const unsigned int nterms = 3;
	std::string fser = "[0]";
	for(unsigned int n=1; n<nterms; n++)
		fser += " + ["+itos(2*n-1)+"]*sin(2*pi*"+itos(n)+"*x) + ["+itos(2*n)+"]*cos(2*pi*"+itos(n)+"*x)";
	printf("Fitter '%s'\n",fser.c_str());
	TF1 fFourier("fFourier",fser.c_str(),-0.5,0.5);
	
	for(Side s = EAST; s <= WEST; ++s) {
		for(AxisDirection d = X_DIRECTION; d <= Y_DIRECTION; ++d) {
			for(unsigned int c=0; c<kMaxCathodes; c++) {
				////////////////////
				// cathode normalization plots
				////////////////////
				
				TH2* hCathNorm =  (TH2*)WCdat.cathNorm[s][d][c]->h[GV_OPEN];
				std::vector<TH1D*> slicefits = replaceFitSlicesY(hCathNorm);
				for(std::vector<TH1D*>::iterator it = slicefits.begin(); it != slicefits.end(); it++) {
					(*it)->SetDirectory(0);
					OM.addObject(*it);
				}
				slicefits[1]->Fit(&fCathCenter,"QR");
				CathodeSeg cs;
				cs.s = s;
				cs.d = d;
				cs.i = c;
				cs.height = float_err(fCathCenter.GetParameter(0),fCathCenter.GetParError(0));
				cs.center = float_err(fCathCenter.GetParameter(1),fCathCenter.GetParError(1));
				cs.width = float_err(fCathCenter.GetParameter(2),fCathCenter.GetParError(2));
				int bm = slicefits[1]->GetMaximumBin();
				cs.max = float_err(slicefits[1]->GetBinContent(bm),slicefits[1]->GetBinError(bm));
				hCathNorm->Draw("Col");
				slicefits[1]->Draw("Same");
				OM.printCanvas(sideSubst("Cathodes/Center_%c",s)+(d==X_DIRECTION?"x":"y")+itos(c));
				
				
				////////////////////
				// cathode positioning distributions by energy
				////////////////////
				
				TH2* hCathHitposDat =  (TH2*)WCdat.cathHitpos[s][d][c]->h[GV_OPEN];
				TH2* hCathHitposSim =  (TH2*)WCsim.cathHitpos[s][d][c]->h[GV_OPEN];
				
				// cathode event fractions
				double nDat = hCathHitposDat->Integral();
				double nSim = hCathHitposSim->Integral();
				cs.fill_frac = nDat/nSim;
				OM.qOut.insert("cathseg",cathseg2sm(cs));
				
				// determine energy re-binning
				std::vector<TH1F*> hEnDats = sliceTH2(*hCathHitposDat,Y_DIRECTION);
				std::vector<TH1F*> hEnSims = sliceTH2(*hCathHitposSim,Y_DIRECTION);
				std::vector<float> enCounts;
				for(unsigned int i=0; i<hEnDats.size(); i++)
					enCounts.push_back(hEnDats[i]->Integral());
				std::vector<unsigned int> part = equipartition(enCounts, 80);
				
				// fit each energy range
				std::vector<TH1*> hToPlot;
				for(unsigned int i = 0; i < part.size()-1; i++) {
					// accumulate over energy range
					TH1F* hEnDat = hEnDats[part[i]];
					TH1F* hEnSim = hEnSims[part[i]];
					float eavg = hEnDat->Integral()*hCathHitposDat->GetYaxis()->GetBinCenter(part[i]+1);
					for(unsigned int j=part[i]+1; j<part[i+1]; j++) {
						hEnDat->Add(hEnDats[j]);
						hEnSim->Add(hEnSims[j]);
						eavg += hEnDats[j]->Integral()*hCathHitposDat->GetYaxis()->GetBinCenter(j+1);
						delete(hEnDats[j]);
						delete(hEnSims[j]);
					}
					eavg /= hEnDat->Integral();
					//printf("Accumulating %u-%u (%.1fkeV) of %i\n",part[i],part[i+1],eavg,(int)hEnDats.size());
					//hEnDat->Rebin();
					//hEnSim->Rebin();
					hEnSim->Scale(hEnDat->Integral()/hEnSim->Integral());
					hEnDat->Divide(hEnSim);
					fixNaNs(hEnDat);
					hEnDat->SetLineColor(i+1);
					OM.addObject(hEnDat);
					OM.addObject(hEnSim);
					
					// fourier components fit
					if(c==1)
						fFourier.SetRange(-0.1,0.5);
					else if(c==kMaxCathodes-2)
						fFourier.SetRange(-0.5,0.1);
					else
						fFourier.SetRange(-0.5,0.5);
					fFourier.SetParameter(0,1.0);
					for(unsigned int n=1; n<=2*nterms; n++)
						fFourier.SetParameter(n,0.);
					fFourier.FixParameter(3,0.);
					fFourier.SetLineColor(i+1);
					hEnDat->Fit(&fFourier,"QR");
					hEnDat->SetMinimum(0);
					hEnDat->SetMaximum(1.5);
					hToPlot.push_back(hEnDat);
					
					// record results
					std::vector<double> terms;
					std::vector<double> dterms;
					double c0 = fFourier.GetParameter(0);
					for(unsigned int n=0; n<2*nterms-1; n++) {
						terms.push_back(fFourier.GetParameter(n)/c0);
						dterms.push_back(fFourier.GetParError(n)/c0);
					}
					Stringmap ff;
					ff.insert("side",sideSubst("%c",s));
					ff.insert("plane",d==X_DIRECTION?"x":"y");
					ff.insert("cathode",c);
					ff.insert("elo",hCathHitposDat->GetYaxis()->GetBinLowEdge(part[i]+1));
					ff.insert("ehi",hCathHitposDat->GetYaxis()->GetBinUpEdge(part[i+1]));
					ff.insert("eavg",eavg);
					ff.insert("terms",vtos(terms));
					ff.insert("dterms",vtos(dterms));
					//ff.display();
					OM.qOut.insert("hitdist",ff);
				}
				drawSimulHistos(hToPlot);
				OM.printCanvas(sideSubst("Cathodes/Positions_%c",s)+(d==X_DIRECTION?"x":"y")+itos(c));
			}
		}
	}
	
	OM.write();
	OM.setWriteRoot(true);
}

void processWirechamberCal(RunNum r0, RunNum r1, unsigned int nrings) {
	std::string basePath = getEnvSafe("UCNA_ANA_PLOTS")+"/PositionMaps/";
	OutputManager OM("NameUnused",basePath);
	std::string readname = itos(r0)+"-"+itos(r1)+"_"+itos(nrings);
	WirechamberAnalyzer WCdat(&OM, "NameUnused", basePath+"/Xenon_"+readname+"/Xenon_"+readname);
	WirechamberAnalyzer WCsim(&OM, "NameUnused", basePath+"/SimXe_"+readname+"/SimXe_"+readname);
	processWirechamberCal(WCdat,WCsim);
}

