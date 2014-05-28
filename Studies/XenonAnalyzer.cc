#include "XenonAnalyzer.hh"
#include "GraphicsUtils.hh"
#include "MultiGaus.hh"
#include "LinHistCombo.hh"
#include "KurieFitter.hh"
#include "PostOfficialAnalyzer.hh"
#include "G4toPMT.hh"

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

//-----------------------------------------------------------

XenonSpectrumPlugin::XenonSpectrumPlugin(RunAccumulator* RA, unsigned int nr): PositionBinnedPlugin(RA,"Xe",nr) {
	// set up histograms
	energySpectrum = registerFGBGPair("hXeSpec", "Xenon energy spectrum", 200, 0, 2000);
	for(Side s = EAST; s <= WEST; ++s) {
		for(unsigned int t=0; t<=nBetaTubes; t++) {
			TH1F hTemplate(("hXe_"+itos(t)).c_str(),"PMT Energy",200,-100,2000);
			sectEnergy[s][t] = allocateSegmentHistograms(hTemplate,AFP_OTHER,s);
		}
	}
	
	// pre-fill sector data
	for(Side s = EAST; s <= WEST; ++s) {
		for(unsigned int m=0; m<=sects.nSectors(); m++) {
			for(unsigned int t=0; t<=nBetaTubes; t++) {
				SectorDat sd;
				sd.s = s;
				sd.t = t;
				sd.m = m;
				sectDat[s][t].push_back(sd);
			}
		}
	}
	
	// load sector data
	if(myA->fIn) {
		QFile qOld(myA->inflname+".txt");
		std::vector<Stringmap> sds = qOld.retrieve("sectDat");
		for(std::vector<Stringmap>::iterator it = sds.begin(); it != sds.end(); it++) {
			SectorDat sd = sm2sd(*it);
			sectDat[sd.s][sd.t][sd.m] = sd;
		}
	}
}

void XenonSpectrumPlugin::fillCoreHists(ProcessedDataScanner& PDS, double weight) {
	Side s = PDS.fSide;
	if(!(PDS.fType == TYPE_0_EVENT && PDS.fPID == PID_BETA)) return;
	if(PDS.radius2(s) < 45*45)
		energySpectrum->h[currentGV]->Fill(PDS.getErecon(),weight);
	unsigned int m = sects.sector(PDS.wires[s][X_DIRECTION].center,PDS.wires[s][Y_DIRECTION].center);
	if(m>=sects.nSectors()) return;
	for(unsigned int t=0; t<nBetaTubes; t++) {
		sectEnergy[s][t][m]->h[currentGV]->Fill(PDS.scints[s].tuben[t].x,weight);
		sectEnergy[s][t][sects.nSectors()]->h[currentGV]->Fill(PDS.scints[s].tuben[t].x,weight);
	}
	sectEnergy[s][nBetaTubes][m]->h[currentGV]->Fill(PDS.scints[s].energy.x,weight);
	sectEnergy[s][nBetaTubes][sects.nSectors()]->h[currentGV]->Fill(PDS.scints[s].energy.x,weight);
}

void XenonSpectrumPlugin::fitSpectrum(TH1* hSpec,SectorDat& sd) {
	
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
	// 2010 analysis fit range was 450-750; expanded for better statistics
	sd.xe_ep = kurieIterator(hSpec,epGuess,NULL,915.,350,850);
}

void XenonSpectrumPlugin::fitSectors() {
	smassert(myA->runCounts.counts.size());
	PMTCalibrator PCal(myA->runCounts.counts.begin()->first);
	printf("\n\n---- Using Calibrator: ----\n");
	PCal.printSummary();
	for(Side s = EAST; s <= WEST; ++s) {
		for(unsigned int t=0; t<=nBetaTubes; t++) {
			for(unsigned int m=0; m<sects.nSectors(); m++) {
				float x = 0;
				float y = 0;
				if(m<sects.nSectors())
					sects.sectorCenter(m,x,y);
				sectDat[s][t][m].eta = PCal.eta(s,t,x,y);
				TH1* hSpec = sectEnergy[s][t][m]->h[GV_OPEN];
				fitSpectrum(hSpec,sectDat[s][t][m]);
				myA->qOut.insert("sectDat",sd2sm(sectDat[s][t][m]));
			}
		}
	}
}

void XenonSpectrumPlugin::calculateResults() {
	for(Side s = EAST; s <= WEST; ++s) {
		for(unsigned int t=0; t<=nBetaTubes; t++) {
			unsigned int m = sects.nSectors();
			TH1* hSpec = sectEnergy[s][t][m]->h[GV_OPEN];
			fitSpectrum(hSpec,sectDat[s][t][m]);
			myA->qOut.insert("sectDat",sd2sm(sectDat[s][t][m]));
			
			AnaNumber AN("XeEndpt");
			AN.s = s;								// PMT side
			AN.n = t;								// PMT number
			AN.value = sectDat[s][t][m].xe_ep.x;	// Xenon spectrum endpoint
			AN.err = sectDat[s][t][m].xe_ep.err;	// endpoint uncertainty (may not be reliable)
			myA->uploadAnaNumber(AN, GV_OPEN, AFP_OTHER);
			
			AN.name = "XeLoPk";
			AN.value = sectDat[s][t][m].low_peak.x;	// Xenon low-energy peak fit center
			AN.err = sectDat[s][t][m].low_peak.err;	// parameter uncertainty
			myA->uploadAnaNumber(AN, GV_OPEN, AFP_OTHER);
		}
	}
}

void XenonSpectrumPlugin::genComparisonPosmap(XenonSpectrumPlugin* XA) {
	// overall energy spectrum
	int b0 = energySpectrum->h[GV_OPEN]->FindBin(400);
	int b1 = energySpectrum->h[GV_OPEN]->FindBin(800);
	energySpectrum->h[GV_OPEN]->Scale(XA->energySpectrum->h[GV_OPEN]->Integral(b0,b1)/energySpectrum->h[GV_OPEN]->Integral(b0,b1));
	drawHistoPair(XA->energySpectrum->h[GV_OPEN],energySpectrum->h[GV_OPEN]);
	printCanvas("DataComparison/hEnergy");
	
	// upload posmap
	std::string pmapname = itos(XA->myA->runCounts.counts.begin()->first)+"-"+itos(XA->myA->runCounts.counts.rbegin()->first)+"/"+itos(time(NULL));
	CalDBSQL* CDBout = CalDBSQL::getCDB(false);
	unsigned int pmid_ep = CDBout->newPosmap("Xe Endpoint "+pmapname,sects.n,sects.r);
	unsigned int pmid_lp = CDBout->newPosmap("Xe Low Peak "+pmapname,sects.n,sects.r);
	float x,y;
	for(unsigned int m=0; m<sects.nSectors(); m++) {
		sects.sectorCenter(m,x,y);
		for(Side s=EAST; s<=WEST; ++s) {
			for(unsigned int t=0; t<nBetaTubes; t++) {
				CDBout->addPosmapPoint(pmid_ep,s,t,m,
									   XA->sectDat[s][t][m].xe_ep.x*XA->sectDat[s][t][m].eta,
									   sectDat[s][t][m].xe_ep.x,x,y);
				CDBout->addPosmapPoint(pmid_lp,s,t,m,
									   XA->sectDat[s][t][m].low_peak.x*XA->sectDat[s][t][m].eta,
									   sectDat[s][t][m].low_peak.x,x,y);
			}
		}
	}
}


//----------------------------------------------------------------

XenonAnalyzer::XenonAnalyzer(OutputManager* pnt, const std::string& nm, const std::string& inflName, unsigned int nrE):
RunAccumulator(pnt,nm,inflName) {
	addPlugin(myXeSpec = new XenonSpectrumPlugin(this,nrE));
	addPlugin(myWG = new MWPCGainPlugin(this));
}

//----------------------------------------------------------------

void process_xenon(RunNum r0, RunNum r1, unsigned int nrings) {
	
	// scan data from each run
	std::vector<std::string> snames;
	OutputManager OM1("NameUnused",getEnvSafe("UCNA_ANA_PLOTS")+"/PositionMaps/SingleRuns/");
	std::vector<RunNum> rlist = CalDBSQL::getCDB()->findRuns("run_type = 'Xenon'",r0,r1);
	for(std::vector<RunNum>::iterator rit = rlist.begin(); rit != rlist.end(); rit++) {
		std::string singleName = "Xenon_"+itos(*rit)+"_"+itos(nrings)+"_"+dtos(PositionBinnedPlugin::fidRadius);
		std::string prevFile = OM1.basePath+"/"+singleName+"/"+singleName;
		snames.push_back(singleName);
		if(r0==r1 || !fileExists(prevFile+".root")) {
			XenonAnalyzer XA1(&OM1, singleName, "", nrings);
			XA1.grouping = GROUP_RUN;
			PostOfficialAnalyzer POA(true);
			POA.redoPositions = true;
			POA.addRun(*rit);
			XA1.loadProcessedData(AFP_OTHER, GV_OPEN, POA);
			XA1.calculateResults();
			XA1.uploadAnaResults();
			POA.writeCalInfo(XA1.qOut,"runcal");
			XA1.write();
			XA1.setWriteRoot(true);
		}
	}
	
	if(r0==r1) return;
	
	// reload data
	printf("Combining Xe runs %i -- %i\n",r0,r1);
	OutputManager OM("NameUnused",getEnvSafe("UCNA_ANA_PLOTS")+"/PositionMaps/");
	
	XenonAnalyzer XA(&OM, "Xenon_"+itos(r0)+"-"+itos(r1)+"_"+itos(nrings), "", nrings);
	XA.grouping = GROUP_RANGE;
	for(std::vector<std::string>::iterator it = snames.begin(); it != snames.end(); it++) {
		std::string prevFile = OM1.basePath+"/"+*it+"/"+*it;
		XenonAnalyzer XA1(&OM1, *it, prevFile);
		XA.addSegment(XA1);
	}
	
	// finish and output
	XA.calculateResults();
	XA.myXeSpec->fitSectors();
	XA.makePlots();
	XA.uploadAnaResults();
	XA.write();
	XA.setWriteRoot(true);
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

std::string simulate_one_xenon(RunNum r, unsigned int nRings, bool forceResim) {
	
	if(!CalDBSQL::getCDB()->findRuns("run_type = 'Xenon'", r, r).size()) return "";	// skip non-Xenon runs
	
	// canonical output paths
	std::string basePath = getEnvSafe("UCNA_ANA_PLOTS")+"/PositionMaps/";
	OutputManager OM1("NameUnused",basePath+"/SingleRunsSim/");
	
	// naming convention for files associated with this run
	std::string singleName = "Xenon_"+itos(r)+"_"+itos(nRings)+"_"+dtos(PositionBinnedPlugin::fidRadius);

	if(forceResim || !fileExists(OM1.basePath+"/"+singleName+"/"+singleName+".root")) {
		// processed data to re-simulate
		XenonAnalyzer XA(&OM1, singleName, basePath+"/SingleRuns/"+singleName+"/"+singleName);
	
		PMTCalibrator PCal(r);
		SimXenonAnalyzer XAM(&OM1,singleName,"",XA.myXeSpec->sects.n);
		XAM.grouping = GROUP_RUN;
		XAM.totalTime[AFP_OTHER][GV_OPEN] += XA.totalTime[AFP_OTHER][GV_OPEN];
		XAM.runTimes += XA.runTimes;
		
		unsigned int nToSim = XA.runCounts[r];
		if(nToSim > 1e6) nToSim = 1e6;
		printf("Data counts West: %f; to sim = %i\n",XA.myXeSpec->energySpectrum->h[GV_OPEN]->Integral(),nToSim);
		
		// simulate for each isotope, in random order to minimize hitting the same file at once from parallel processes
		std::vector<SimXenonAnalyzer*> XAMi;
		std::vector<std::string> isotsIn;
		LinHistCombo LHC;
		
		isotsIn.push_back("Xe125_1-2+");
		isotsIn.push_back("Xe129_11-2-");
		isotsIn.push_back("Xe131_11-2-");
		isotsIn.push_back("Xe133_3-2+");
		isotsIn.push_back("Xe133_11-2-");
		isotsIn.push_back("Xe135_3-2+");
		
		// runs promptly after activation with short-lived peak
		if((14264 <= r && r <= 14273) || (15991 <= r && r <= 16010))
			isotsIn.push_back("Xe135_11-2-");
		
		// shorter-lived high-energy beta spectrum component
		int b1 = XA.myXeSpec->energySpectrum->h[GV_OPEN]->FindBin(1075);
		int b2 = XA.myXeSpec->energySpectrum->h[GV_OPEN]->FindBin(1175);
		if(XA.myXeSpec->energySpectrum->h[GV_OPEN]->Integral(b1,b2) > 100)
			isotsIn.push_back("Xe137_7-2-");
		
		
		// randomize order of isotope simulation to pick different simulation files for each one
		srand((unsigned int)time(NULL)); // set new random seed from clock
		std::vector<unsigned int> p = randomPermutation((unsigned int)isotsIn.size());
		std::vector<std::string> isots;
		for(unsigned int n=0; n<isotsIn.size(); n++)
			isots.push_back(isotsIn[p[n]]);
		
		for(unsigned int n=0; n<isots.size(); n++) {
			printf("Simulating for component %s...\n",isots[n].c_str());
			XAMi.push_back(new SimXenonAnalyzer(&OM1,singleName+"_"+isots[n],"",XA.myXeSpec->sects.n));
			G4toPMT G2P;
			G2P.runCathodeSim();
			G2P.setCalibrator(PCal);
			std::string simFile = getEnvSafe("UCNA_CALSRC_SIMS")+"/"+isots[n]+"/analyzed_*";
			G2P.addFile(simFile);
			
			unsigned int nToSimI = nToSim*(isots[n]=="Xe135_3-2+"? 1.0 : isots[n]=="Xe133_3-2+"?  0.05 : 0.25);
			if(nToSimI > 0.5*G2P.nEvents) nToSimI = 0.5*G2P.nEvents;
			XAMi.back()->loadSimData(G2P, nToSimI);
			
			AnaNumber AN("XeSimCounts_"+isots[n]);
			AN.value = nToSimI;		// number of Type 0 events simulated for isotope
			AN.err = G2P.nSimmed;	// total number of events run to meet requested count
			XAM.uploadAnaNumber(AN, GV_OPEN, AFP_OTHER);
			
			LHC.addTerm(XAMi.back()->myXeSpec->energySpectrum->h[GV_OPEN]);
			printf("Done.\n");
		}
		
		// determine spectrum composition, and sum re-scaled spectra
		LHC.Fit(XA.myXeSpec->energySpectrum->h[GV_OPEN],50,1000);
		for(unsigned int i=0; i<LHC.coeffs.size(); i++) {
			XAMi[i]->scaleData(LHC.coeffs[i]);
			XAM.addSegment(*XAMi[i]);
			
			AnaNumber AN("XeComp_"+isots[i]);
			AN.value = XAMi[i]->myXeSpec->energySpectrum->h[GV_OPEN]->Integral();	// Type 0 counts from this isotope in estimated composition
			AN.err = (LHC.dcoeffs[i]/LHC.coeffs[i])*AN.value;						// isotope composition fit error on counts
			XAM.uploadAnaNumber(AN, GV_OPEN, AFP_OTHER);
			
			delete(XAMi[i]);
		}

		XAM.qOut.insert("runcal",PCal.calSummary());
		XAM.calculateResults();
		XAM.compareMCtoData(XA);
		XAM.uploadAnaResults();
		XAM.write();
		XAM.setWriteRoot(true);
	}
	
	return singleName;
}

void combine_xenon_sims(RunNum r0, RunNum r1, unsigned int nRings) {
	
	// canonical paths
	std::string basePath = getEnvSafe("UCNA_ANA_PLOTS")+"/PositionMaps/";
	OutputManager OM("NameUnused",basePath);
	OutputManager OM1("NameUnused",basePath+"/SingleRunsSim/");
	std::string readname = itos(r0)+"-"+itos(r1)+"_"+itos(nRings);
		
	// merge simulated data
	SimXenonAnalyzer XAM(&OM, "SimXe_"+readname, "", nRings);
	XAM.grouping = GROUP_RANGE;
	XAM.isSimulated = true;
	for(RunNum rn = r0; rn <= r1; rn++) {
		std::string singleName = simulate_one_xenon(rn, nRings, false);
		if(!singleName.size()) continue;
		std::string prevFile = OM1.basePath+"/"+singleName+"/"+singleName;
		SimXenonAnalyzer XAM1(&OM1, singleName, prevFile);
		XAM.addSegment(XAM1);
	}
	
	// finish and output
	XAM.calculateResults();
	XAM.myXeSpec->fitSectors();
	XAM.makePlots();
	XAM.uploadAnaResults();
	XAM.write();
	XAM.setWriteRoot(true);
}

void xenon_posmap(RunNum r0, RunNum r1, unsigned int nRings) {
	std::string basePath = getEnvSafe("UCNA_ANA_PLOTS")+"/PositionMaps/";
	OutputManager OM("NameUnused",basePath);
	std::string readname = itos(r0)+"-"+itos(r1)+"_"+itos(nRings);
	
	// read in data
	XenonAnalyzer XAdat(&OM, "Xenon_"+readname, basePath+"/Xenon_"+readname+"/Xenon_"+readname);
	
	// read in simulation
	XenonAnalyzer XAsim(&OM, "SimXe_"+readname, basePath+"/SimXe_"+readname+"/SimXe_"+readname);
	XAsim.grouping = GROUP_RANGE;
	
	// data comparison / posmap generation
	XAsim.compareMCtoData(XAdat);
	XAsim.myXeSpec->genComparisonPosmap(XAdat.myXeSpec);
	XAsim.uploadAnaResults();
}

