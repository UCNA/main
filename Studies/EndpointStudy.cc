#include "EndpointStudy.hh"
#include "CalDBSQL.hh"
#include "MultiGaus.hh"
#include "TH1toPMT.hh"
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
RunAccumulator(pnt,nm,infl), sects(nr,r) {
	
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
	TH2F hPositionsTemplate("hPostions","Event Positions",200,-65,65,200,-65,65);
	for(Side s = EAST; s <= WEST; ++s) {
		hitPos[s] = registerFGBGPair(hPositionsTemplate,AFP_OTHER,s);
		energySpectrum[s] = registerFGBGPair("hEnergy","Combined Energy",200,-100,1200,AFP_OTHER,s);
		energySpectrum[s].h[GV_OPEN]->SetLineColor(2+2*s);
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
	if(PDS.radius2(s) <= 25*25)
		energySpectrum[s].h[currentGV]->Fill(PDS.scints[s].energy.x,weight);
	if(PDS.getEtrue() > 225)
		sectAnode[s][m].h[currentGV]->Fill(PDS.mwpcs[s].anode,weight);
}

void PositionBinner::calculateResults() {
	
	assert(runCounts.counts.size());
	PMTCalibrator PCal(runCounts.counts.begin()->first,CalDBSQL::getCDB());
	printf("\n\n---- Using Calibrator: ----\n");
	PCal.printSummary();
	
	TF1 gausFit("gasufit","gaus",0,500);
	
	for(Side s = EAST; s <= WEST; ++s) {
		for(unsigned int t=0; t<nBetaTubes; t++) {
			for(unsigned int m=0; m<sects.nSectors(); m++) {
				
				TH1* hSpec = sectEnergy[s][t][m].h[GV_OPEN];
				hSpec->SetLineColor(2+t);
				float x,y;
				sects.sectorCenter(m,x,y);
				float eta = PCal.eta(s,t,x,y);
				sectDat[s][t][m].eta = eta;
				
				//----------------------
				// Low peak fit
				//----------------------
				gausFit.SetLineColor(2+t);
				if(!iterGaus(hSpec,&gausFit,3,hSpec->GetBinCenter(hSpec->GetMaximumBin()),100,1.0)) {
					sectDat[s][t][m].low_peak = float_err(gausFit.GetParameter(1),gausFit.GetParError(1));
					sectDat[s][t][m].low_peak_width = float_err(gausFit.GetParameter(2),gausFit.GetParError(2));
				} else {
					sectDat[s][t][m].low_peak = sectDat[s][t][m].low_peak_width = 0;
				}
				
				//----------------------
				// 915keV endpoint fit
				//----------------------
				sectDat[s][t][m].xe_ep = kurieIterator(hSpec,6.5*sectDat[s][t][m].low_peak.x,NULL,915.,400,800);
				
				qOut.insert("sectDat",sd2sm(sectDat[s][t][m]));
			}
		}
	}	
}

void PositionBinner::compareMCtoData(RunAccumulator& OAdata, float simfactor) {
	
	PositionBinner* PB = (PositionBinner*)&OAdata;
	defaultCanvas->cd();
	
	// overall energy spectrum
	for(Side s=EAST; s<=WEST; ++s) {
		energySpectrum[s].h[GV_OPEN]->Scale(PB->energySpectrum[s].h[GV_OPEN]->GetMaximum()/energySpectrum[s].h[GV_OPEN]->GetMaximum());
		energySpectrum[s].h[GV_OPEN]->SetLineColor(4);
		energySpectrum[s].h[GV_OPEN]->Draw();
		PB->energySpectrum[s].h[GV_OPEN]->SetLineColor(2);
		PB->energySpectrum[s].h[GV_OPEN]->Draw("Same");
		printCanvas(sideSubst("Comparison/hEnergy_%c",s));		
	}
	
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
	hToPlot.push_back(energySpectrum[EAST].h[GV_OPEN]);
	hToPlot.push_back(energySpectrum[WEST].h[GV_OPEN]);
	drawSimulHistos(hToPlot);
	printCanvas("hEnergy");
		
	for(Side s = EAST; s <= WEST; ++s) {
		// positions
		hitPos[s].h[GV_OPEN]->Draw("COL");
		drawSectors(sects,6);
		labelSectors(sects,6);
		printCanvas(sideSubst("hPos_%c",s));
		// energy in each sector
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


void process_xenon(RunNum r0, RunNum r1, unsigned int nrings) {	
	
	double fidRadius = 52;
	
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
			PB1.write();
			PB1.setWriteRoot(true);
		}
	}
	
	if(r0==r1) return;
	
	// reload data
	OutputManager OM("NameUnused",getEnvSafe("UCNA_ANA_PLOTS")+"/PositionMaps/");
	PositionBinner PB(&OM, std::string("Xenon_")+itos(r0)+"-"+itos(r1), fidRadius, nrings);	
	for(std::vector<std::string>::iterator it = snames.begin(); it != snames.end(); it++) {
		std::string prevFile = OM1.basePath+"/"+*it+"/"+*it;
		PositionBinner PB1(&OM1, *it, 0, 0, prevFile);
		PB.addSegment(PB1);
	}
	
	// finish and output
	PB.calculateResults();
	PB.makePlots();
	PB.write();
	PB.setWriteRoot(true);	
}

double avgSmear(PMTCalibrator& PCal, Side s, SectorCutter& sects) {
	double sm = 0;
	double e0 = 500.0;
	float x,y;
	for(unsigned int m=0; m<sects.nSectors(); m++) {
		sects.sectorCenter(m,x,y);
		sm+=PCal.nPE(s,nBetaTubes,e0,x,y)/e0;
	}
	return sm/sects.nSectors();
}

std::string simulate_one_xenon(RunNum r, OutputManager& OM1, PositionBinner& PB, float simFactor, bool forceResim=false) {
	std::string singleName = std::string("Xenon_")+itos(r)+"_"+itos(PB.sects.n)+"_"+dtos(PB.sects.r);
	std::string prevFile = OM1.basePath+"/"+singleName+"/"+singleName;
	if(forceResim || !fileExists(prevFile+".root")) {
		PositionBinner PBM(&OM1,singleName,PB.sects.r,PB.sects.n);
		SectPosGen SPG(PB.sects);
		PMTCalibrator PCal(r,CalDBSQL::getCDB());
		for(Side s = EAST; s <= WEST; ++s) {
			TH1toPMT t2p(PB.energySpectrum[s].h[GV_OPEN],&SPG);
			t2p.setCalibrator(PCal);
			t2p.setAFP(AFP_OTHER);
			t2p.PGen[s].presmear = avgSmear(PCal,s,PB.sects);
			printf("Setting pre-smearing to %g PE/keV...\n",t2p.PGen[s].presmear);
			//t2p.PGen[s].larmorField = sqrt(0.6);
			t2p.genside = s;
			for(unsigned int m=0; m<PB.sects.nSectors(); m++) {
				printf("Simulating sector %c%i...\n",sideNames(s),m);
				SPG.m = m;
				t2p.nToSim=simFactor*0.5*PB.runCounts[r]*SPG.sects.sectorArea(m)/SPG.sects.totalArea();
				PBM.loadSimData(t2p,t2p.nToSim);
			}
		}
		PBM.write();
		PBM.setWriteRoot(true);
	}
	return singleName;
}

void simulate_xenon(RunNum r0, RunNum r1, RunNum rsingle) {
	
	// read in comparison data
	OutputManager OM("NameUnused",getEnvSafe("UCNA_ANA_PLOTS")+"/PositionMaps/");
	std::string readname = std::string("Xenon_")+itos(r0)+"-"+itos(r1);
	PositionBinner PB(&OM, std::string("Xenon_")+itos(r0)+"-"+itos(r1), 0, 0, OM.basePath+"/"+readname+"/"+readname);
	
	// MC data for each run
	printf("Simulating for position map %i-%i (%g,%i)\n",r0,r1,PB.sects.r, PB.sects.n);
	OutputManager OM1("NameUnused",getEnvSafe("UCNA_ANA_PLOTS")+"/PositionMaps/SingleRunsSim/");
	float simFactor = 4.0;
	if(rsingle) {
		simulate_one_xenon(rsingle,OM1,PB,simFactor,true);
		return;
	}
	std::vector<std::string> snames;
	for(std::map<RunNum,double>::const_iterator rit = PB.runCounts.counts.begin(); rit != PB.runCounts.counts.end(); rit++)
		snames.push_back(simulate_one_xenon(rit->first,OM1,PB,simFactor));
	
	// reload data
	PositionBinner PBM(&OM, std::string("SimXe_")+itos(r0)+"-"+itos(r1), PB.sects.r, PB.sects.n);
	for(std::vector<std::string>::iterator it = snames.begin(); it != snames.end(); it++) {
		std::string prevFile = OM1.basePath+"/"+*it+"/"+*it;
		PositionBinner PBM1(&OM1, *it, 0, 0, prevFile);
		PBM.addSegment(PBM1);
	}
	
	// finish and output
	PBM.calculateResults();
	PBM.makePlots();
	PBM.compareMCtoData(PB,simFactor);
	PBM.write();
	PBM.setWriteRoot(true);
}
