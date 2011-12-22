#include "EndpointStudy.hh"
#include "CalDBSQL.hh"
#include "MultiGaus.hh"
#include "TH1toPMT.hh"
#include "PostOfficialAnalyzer.hh"
#include <TStyle.h>

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
	for(Side s = EAST; s <= WEST; ++s) {
		energySpectrum[s] = registerFGBGPair("hEnergy","Combined Energy",200,-100,1200,AFP_OTHER,s);
		energySpectrum[s].h[GV_OPEN]->SetLineColor(2+2*s);
		for(unsigned int m=0; m<getNSectors(); m++) {
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
	unsigned int m = getSector(PDS.wires[s][X_DIRECTION].center,PDS.wires[s][Y_DIRECTION].center);
	if(m>=getNSectors()) return;
	for(unsigned int t=0; t<nBetaTubes; t++)
		sectEnergy[s][t][m].h[currentGV]->Fill(PDS.scints[s].tuben[t].x,weight);
	if(PDS.radius2(s) <= 25*25)
		energySpectrum[s].h[currentGV]->Fill(PDS.scints[s].energy.x,weight);
	if(PDS.getEtrue() > 225)
		sectAnode[s][m].h[currentGV]->Fill(PDS.mwpcs[s].anode,weight);
}

void PositionBinner::calculateResults() {
	
	assert(runTimes.counts.size());
	PMTCalibrator PCal(runTimes.counts.begin()->first,CalDBSQL::getCDB());
	printf("\n\n---- Using Calibrator: ----\n");
	PCal.printSummary();
	
	TF1 gausFit("gasufit","gaus",0,500);
	
	for(Side s = EAST; s <= WEST; ++s) {
		for(unsigned int t=0; t<nBetaTubes; t++) {
			for(unsigned int m=0; m<getNSectors(); m++) {
				
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
	// TODO
	// upload as posmap
	// load & draw posmap
}

void PositionBinner::makePlots() {
	
	std::vector<TH1*> hToPlot;
	
	// overall energy spectrum
	hToPlot.push_back(energySpectrum[EAST].h[GV_OPEN]);
	hToPlot.push_back(energySpectrum[WEST].h[GV_OPEN]);
	drawSimulHistos(hToPlot);
	printCanvas("hEnergy");
	
	// energy in each sector
	for(Side s = EAST; s <= WEST; ++s) {
		for(unsigned int m=0; m<getNSectors(); m++) {
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
	// set up output
	OutputManager OM("NameUnused","../PostPlots/PositionMaps/");
	PositionBinner PB(&OM, std::string("Xenon_")+itos(r0)+"-"+itos(r1), 55, nrings);
	
	// load data from each run
	for(RunNum r = r0; r <= r1; r++) {
		PostOfficialAnalyzer POA(true);
		POA.addRun(r);
		PB.loadProcessedData(AFP_OTHER, GV_OPEN, POA);
	}
	
	// finish and output
	PB.calculateResults();
	PB.makePlots();
	PB.write();
	PB.setWriteRoot(true);	
}

void simulate_xenon(RunNum r0, RunNum r1) {
	// set up output
	OutputManager OM("NameUnused","../PostPlots/PositionMaps/");
	
	// read in comparison data
	PositionBinner PB(&OM, std::string("Xenon_")+itos(r0)+"-"+itos(r1), 0, 0, std::string("Xenon_")+itos(r0)+"-"+itos(r1));
	
	// MC data
	PositionBinner PBM(&OM, std::string("SimXe_")+itos(r0)+"-"+itos(r1), 0, 0); //TODO
	
	// finish and output
	PBM.calculateResults();
	PBM.makePlots();
	PBM.write();
	PBM.setWriteRoot(true);
}
