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

PositionBinner::PositionBinner(OutputManager* pnt, const std::string& nm, float r, unsigned int nr, const std::string& infl):
RunAccumulator(pnt,nm,infl), sects(nr,r), sectorPlots(false) {
	
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
	energySpectrum.h[GV_OPEN]->SetLineColor(2);
	TH2F hPositionsTemplate("hPostions","Event Positions",200,-60,60,200,-60,60);
	for(Side s = EAST; s <= WEST; ++s) {
		hitPos[s] = registerFGBGPair(hPositionsTemplate,AFP_OTHER,s);
		for(unsigned int m=0; m<sects.nSectors(); m++) {
			sectAnode[s].push_back(registerFGBGPair(std::string("hAnode_Sector_")+itos(m),std::string("Sector ")+itos(m)+" Anode",200,0,4000,AFP_OTHER,s));
			for(unsigned int t=0; t<nBetaTubes; t++) {
				sectEnergy[s][t].push_back(registerFGBGPair(std::string("hTuben_")+itos(t)+"_Sector_"+itos(m),
															std::string("Sector ")+itos(m)+" Tube "+itos(t)+" Energy",200,-100,2000,AFP_OTHER,s));
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
	const Side s = PDS.fSide;
	if(!(PDS.fType == TYPE_0_EVENT && PDS.fPID == PID_BETA && (s==EAST||s==WEST))) return;
	((TH2F*)(hitPos[s].h[currentGV]))->Fill(PDS.wires[s][X_DIRECTION].center,PDS.wires[s][Y_DIRECTION].center,weight);
	unsigned int m = sects.sector(PDS.wires[s][X_DIRECTION].center,PDS.wires[s][Y_DIRECTION].center);
	if(m>=sects.nSectors()) return;
	for(unsigned int t=0; t<nBetaTubes; t++)
		sectEnergy[s][t][m].h[currentGV]->Fill(PDS.scints[s].tuben[t].x,weight);
	if(PDS.radius2(s) <= 45*45)
		energySpectrum.h[currentGV]->Fill(PDS.getEtrue(),weight);
	if(PDS.getEtrue() > 225)
		sectAnode[s][m].h[currentGV]->Fill(PDS.mwpcs[s].anode,weight);
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
				TH1* hSpec = sectEnergy[s][t][m].h[GV_OPEN];
				fitSpectrum(hSpec,sectDat[s][t][m]);
				qOut.insert("sectDat",sd2sm(sectDat[s][t][m]));
			}
		}
	}	
}

void PositionBinner::calculateResults() {
	for(Side s = EAST; s <= WEST; ++s) {
		for(unsigned int t=0; t<nBetaTubes; t++) {
			hTuben[s][t] = (TH1F*)addObject(sectEnergy[s][t][0].h[GV_OPEN]->Clone((sideSubst("hTuben_%c",s)+itos(t)).c_str()));
			hTuben[s][t]->Reset();
			for(unsigned int m=0; m<sects.nSectors(); m++) {
				if(sects.sectorCenterRadius(m) < 45.)
					hTuben[s][t]->Add(sectEnergy[s][t][m].h[GV_OPEN]);
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
	
	PositionBinner* PB = (PositionBinner*)&OAdata;
	defaultCanvas->cd();
	
	// overall energy spectrum
	int b0 = PB->energySpectrum.h[GV_OPEN]->FindBin(400);
	int b1 = PB->energySpectrum.h[GV_OPEN]->FindBin(800);
	energySpectrum.h[GV_OPEN]->Scale(PB->energySpectrum.h[GV_OPEN]->Integral(b0,b1)/energySpectrum.h[GV_OPEN]->Integral(b0,b1));
	drawHistoPair(PB->energySpectrum.h[GV_OPEN],energySpectrum.h[GV_OPEN]);
	printCanvas("Comparison/hEnergy");		
	
	// upload posmap
	std::string pmapname = itos(PB->runCounts.counts.begin()->first)+"-"+itos(PB->runCounts.counts.rbegin()->first)+"/"+itos(time(NULL));
	CalDBSQL* CDBout = CalDBSQL::getCDB(false);
	unsigned int pmid_ep = CDBout->newPosmap(std::string("Xe Endpoint ")+pmapname,sects.n,sects.r);
	unsigned int pmid_lp = CDBout->newPosmap(std::string("Xe Low Peak ")+pmapname,sects.n,sects.r);
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
	
	std::vector<TH1*> hToPlot;
	
	// overall energy spectrum
	defaultCanvas->cd();
	energySpectrum.h[GV_OPEN]->Draw();
	printCanvas("hEnergy");
	
	for(Side s = EAST; s <= WEST; ++s) {
		// positions
		hitPos[s].h[GV_OPEN]->Draw("COL");
		drawSectors(sects,6);
		//labelSectors(sects,6);
		printCanvas(sideSubst("hPos_%c",s));
		
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
					hToPlot.push_back(sectEnergy[s][t][m].h[GV_OPEN]);
				drawSimulHistos(hToPlot);
				for(unsigned int t=0; t<nBetaTubes; t++)
					drawVLine(sectDat[s][t][m].xe_ep.x, defaultCanvas, 2+t);
				printCanvas(sideSubst("SectorEnergy/h_%c_",s)+itos(m));
			}
		}
	}
}


void process_xenon(RunNum r0, RunNum r1, unsigned int nrings) {	
	
	double fidRadius = 50;
	
	// scan data from each run
	std::vector<std::string> snames;
	OutputManager OM1("NameUnused",getEnvSafe("UCNA_ANA_PLOTS")+"/PositionMaps/SingleRuns/");
	for(RunNum r = r0; r <= r1; r++) {
		std::string singleName = std::string("Xenon_")+itos(r)+"_"+itos(nrings)+"_"+dtos(fidRadius);
		std::string prevFile = OM1.basePath+"/"+singleName+"/"+singleName;
		snames.push_back(singleName);
		if(r0==r1 || !fileExists(prevFile+".root")) {
			PositionBinner PB1(&OM1, singleName, fidRadius, nrings);
			PostOfficialAnalyzer POA(true);
			POA.addRun(r);
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
	PositionBinner PB(&OM, std::string("Xenon_")+itos(r0)+"-"+itos(r1)+"_"+itos(nrings), fidRadius, nrings);	
	for(std::vector<std::string>::iterator it = snames.begin(); it != snames.end(); it++) {
		std::string prevFile = OM1.basePath+"/"+*it+"/"+*it;
		PositionBinner PB1(&OM1, *it, 0, 0, prevFile);
		PB.addSegment(PB1);
	}
	
	// finish and output
	PB.calculateResults();
	PB.fitSectors();
	PB.makePlots();
	PB.write();
	PB.setWriteRoot(true);	
}

std::string simulate_one_xenon(RunNum r, OutputManager& OM1, PositionBinner& PB, float simFactor, bool forceResim=false) {
	std::string singleName = std::string("Xenon_")+itos(r)+"_"+itos(PB.sects.n)+"_"+dtos(PB.sects.r);
	std::string prevFile = OM1.basePath+"/"+singleName+"/"+singleName;
	if(forceResim || !fileExists(prevFile+".root")) {
		PMTCalibrator PCal(r);
		PositionBinner PBM(&OM1,singleName,PB.sects.r,PB.sects.n);
		PBM.totalTime[AFP_OTHER][GV_OPEN] += PB.totalTime[AFP_OTHER][GV_OPEN];
		PBM.runTimes += PB.runTimes;
		
		unsigned int nToSim=simFactor*PB.runCounts[r];
		printf("Data counts West: %f; to sim = %i\n",PB.energySpectrum.h[GV_OPEN]->Integral(),nToSim);
		
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
		if(15991 <= r && r <= 16010)
			isots.push_back("Xe135_11-2-");
		int b1 = PB.energySpectrum.h[GV_OPEN]->FindBin(1075);
		int b2 = PB.energySpectrum.h[GV_OPEN]->FindBin(1175);
		if(PB.energySpectrum.h[GV_OPEN]->Integral(b1,b2) > 100)
			isots.push_back("Xe137_7-2-");
		
		for(unsigned int n=0; n<isots.size(); n++) {
			printf("Simulating for component %s...\n",isots[n].c_str());
			PBMi.push_back(new PositionBinner(&OM1,singleName+"_"+isots[n],PB.sects.r,PB.sects.n));
			//G4SegmentMultiplier GSM(PB.sects);
			G4SegmentMultiplier GSM(SectorCutter(4,52.));
			GSM.setCalibrator(PCal);
			std::string simFile = "/home/mmendenhall/geant4/output/WideKev_"+isots[n]+"/analyzed_";
			unsigned int nTot = 18;
			unsigned int stride = 7;
			for(unsigned int i=0; i<stride; i++)
				GSM.addFile(simFile+itos((stride*r+i)%nTot)+".root");
			PBMi.back()->loadSimData(GSM, nToSim*(isots[n]=="Xe135_3-2+"?1.5:0.25));
			LHC.addTerm(PBMi.back()->energySpectrum.h[GV_OPEN]);
			printf("Done.\n");
		}
		
		// determine spectrum composition and accumulate segments
		LHC.Fit(PB.energySpectrum.h[GV_OPEN],50,1000);
		std::vector<double> counts;
		for(unsigned int i=0; i<LHC.coeffs.size(); i++) {
			PBMi[i]->scaleData(LHC.coeffs[i]);
			counts.push_back(PBMi[i]->energySpectrum.h[GV_OPEN]->Integral());
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

void simulate_xenon(RunNum r0, RunNum r1, RunNum rsingle, unsigned int nRings) {
	
	// read in comparison data
	std::string basePath = getEnvSafe("UCNA_ANA_PLOTS")+"/PositionMaps/";
	OutputManager OM("NameUnused",basePath);
	std::string readname = std::string("Xenon_")+itos(r0)+"-"+itos(r1)+"_"+itos(nRings);
	PositionBinner PB(&OM, std::string("Xenon_")+itos(r0)+"-"+itos(r1), 0, 0, basePath+"/"+readname+"/"+readname);
	
	// MC data for each run
	printf("Simulating for position map %i-%i\n",r0,r1);
	OutputManager OM1("NameUnused",basePath+"/SingleRunsSim/");
	float simFactor = 4.0;
	if(rsingle) {
		std::string singleName = std::string("Xenon_")+itos(rsingle)+"_"+itos(PB.sects.n)+"_"+dtos(PB.sects.r);
		PositionBinner PB1(&OM, singleName, 0, 0, basePath+"/SingleRuns/"+singleName+"/"+singleName);
		simulate_one_xenon(rsingle,OM1,PB1,simFactor,true);
		return;
	}
	std::vector<std::string> snames;
	for(std::map<RunNum,double>::const_iterator rit = PB.runCounts.counts.begin(); rit != PB.runCounts.counts.end(); rit++)
		snames.push_back(simulate_one_xenon(rit->first,OM1,PB,simFactor));
	
	
	// reload data
	PositionBinner PBM(&OM, std::string("SimXe_")+itos(r0)+"-"+itos(r1)+"_"+itos(nRings), PB.sects.r, PB.sects.n);
	PBM.isSimulated = true;
	for(std::vector<std::string>::iterator it = snames.begin(); it != snames.end(); it++) {
		std::string prevFile = OM1.basePath+"/"+*it+"/"+*it;
		PositionBinner PBM1(&OM1, *it, 0, 0, prevFile);
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
