#include "ReSource.hh"
#include "PostOfficialAnalyzer.hh"
#include "PMTGenerator.hh"
#include "TSpectrumUtils.hh"
#include "GraphicsUtils.hh"
#include "G4toPMT.hh"
#include "PathUtils.hh"
#include "SourceDBSQL.hh"
#include "CalDBSQL.hh"
#include "SMExcept.hh"
#include <TSpectrum.h>
#include <TSpectrum2.h>
#include <utility>


ReSourcer::ReSourcer(OutputManager* O, const Source& s, PMTCalibrator* P):
OM(O), mySource(s), PCal(P), dbgplots(false), simMode(false),
nBins(300), eMin(-100), eMax(2000), pkMin(0.0), nSigma(2.0) {
	
	// source-dependent ranges
	if(mySource.t == "Bi207") {
		pkMin = 200;
		nSigma = 1.5;
		addCorrFit(400,475);
		addCorrFit(900,1100);
	} else if(mySource.t == "Sn113") {
		nBins = 150;
		eMin = -50;
		eMax = 1000;
		addCorrFit(300,450);
	} else if(mySource.t == "Ce139") {
		nBins = 100;
		eMin = -20;
		eMax = 500;
		addCorrFit(70,140);
	} else if (mySource.t == "Cd109") {
		nBins = 100;
		eMin = -20;
		eMax = 300;
	}
	
	// set up histograms
	for(unsigned int t=0; t<=nBetaTubes; t++) {
		for(unsigned int tp=TYPE_0_EVENT; tp<=TYPE_III_EVENT; tp++) {
			hTubes[t][tp] = OM->registeredTH1F(s.name()+(t==nBetaTubes?"_Combined":"_Tube_"+itos(t))+"_type_"+itos(tp),mySource.t+" energy spectrum",
											   tp==TYPE_0_EVENT?nBins:nBins/4,eMin,eMax);
			if(tp==TYPE_0_EVENT) hTubes[t][tp]->Sumw2();
			hTubes[t][tp]->GetXaxis()->SetTitle("Scintillator Visible Energy [keV]");
			hTubes[t][tp]->GetYaxis()->SetTitle("Event Rate [Hz/keV]");
			hTubes[t][tp]->GetYaxis()->SetTitleOffset(1.25);
		}
		hTubesRaw[t] = NULL;
		if(t<nBetaTubes && !simMode) {
			hTubesRaw[t] = OM->registeredTH1F(s.name()+"_"+itos(t)+"_ADC",mySource.t+" ADC spectrum",nBins,
											  PCal?-PCal->invertCorrections(mySource.mySide, t, -eMin, mySource.x, mySource.y, 0.0):2*eMin,
											  PCal?PCal->invertCorrections(mySource.mySide, t, eMax, mySource.x, mySource.y, 0.0):2*eMax);
			hTubesRaw[t]->Sumw2();
			hTubesRaw[t]->GetXaxis()->SetTitle("ADC Channels");
			hTubesRaw[t]->GetYaxis()->SetTitle("Event Rate [Hz/keV]");
		}
		for(unsigned int t2=0; t2<nBetaTubes; t2++) {
			if(t==nBetaTubes || t==t2) continue;			
			pPMTCorr[t][t2] = new TProfile((itos(t)+"_vs_"+itos(t2)).c_str(),
										   "PMT Correlation",nBins,eMin,eMax);
			pPMTCorr[t][t2]->SetMinimum(eMin);
			pPMTCorr[t][t2]->SetMaximum(eMax);
			OM->addObject(pPMTCorr[t][t2]);
		}
	}	
	for(AxisDirection d = X_DIRECTION; d <= Y_DIRECTION; ++d) {
		hitPos[d] = OM->registeredTH1F(s.name()+"_hits_profile_"+itos(d),"Hit Positions",300,-10,10);
		hitPos[d]->SetLineColor(2+2*d);
	}
}

unsigned int ReSourcer::fill(const ProcessedDataScanner& P) {
	EventType tp = P.fType;
	Side s = P.fSide;
	
	if(tp > TYPE_III_EVENT || s != mySource.mySide) return 0;
	float x = P.wires[P.fSide][X_DIRECTION].center;
	float y = P.wires[P.fSide][Y_DIRECTION].center;
	if(!mySource.inSourceRegion(x,y,4.0)) return 0;
	if(P.fPID != PID_BETA) return 0;
	
	if(tp==TYPE_0_EVENT) {
		hitPos[X_DIRECTION]->Fill(x-mySource.x,P.physicsWeight);
		hitPos[Y_DIRECTION]->Fill(y-mySource.y,P.physicsWeight);
		for(unsigned int t1=0; t1<nBetaTubes; t1++)
			for(unsigned int t2=0; t2<nBetaTubes; t2++)
				if(t1!=t2)
					pPMTCorr[t1][t2]->Fill(P.scints[s].tuben[t1].x,P.scints[s].tuben[t2].x);
	}
	for(unsigned int t=0; t<nBetaTubes; t++) {
		if(P.scints[s].adc[t]>3750 || P.scints[s].tuben[t].x < 1.0)
			continue;
		hTubes[t][tp]->Fill(P.scints[s].tuben[t].x,P.physicsWeight);
		if(PCal && !simMode && tp==TYPE_0_EVENT)
			hTubesRaw[t]->Fill(P.scints[s].adc[t]*PCal->gmsFactor(mySource.mySide,t,P.runClock[BOTH]),P.physicsWeight);
	}
	if(tp==TYPE_0_EVENT)
		hTubes[nBetaTubes][tp]->Fill(P.scints[s].energy.x,P.physicsWeight);
	else
		hTubes[nBetaTubes][tp]->Fill(P.getEnergy(),P.physicsWeight);
	return 1;
}

void ReSourcer::findSourcePeaks(float runtime) {
	
	std::vector<SpectrumPeak> expectedPeaks = mySource.getPeaks();
	printf("Fitting peaks for %s source (%i hits, %.2f hours, %i peaks)...\n",mySource.name().c_str(),counts(),runtime/3600.0,(int)expectedPeaks.size());
	OM->defaultCanvas->SetLogy(true);
	
	// normalize to rate Hz/keV; record rates
	for(unsigned int t=0; t<=nBetaTubes; t++) {
		double nType0 = hTubes[t][TYPE_0_EVENT]->Integral();
		for(unsigned int tp=TYPE_0_EVENT; tp<=TYPE_III_EVENT; tp++) {
			if(t==nBetaTubes) {
				Stringmap m;
				m.insert("type",itos(tp));
				m.insert("sID",mySource.sID);
				m.insert("name",mySource.name());
				m.insert("simulated",simMode?"yes":"no");			
				m.insert("side",sideSubst("%c",mySource.mySide));
				m.insert("counts",hTubes[t][tp]->Integral());
				m.insert("rate",runtime?hTubes[t][tp]->Integral()/runtime:0);
				m.insert("type0frac",hTubes[t][tp]->Integral()/nType0);
				OM->qOut.insert("rate",m);
			}
			hTubes[t][tp]->Scale(1.0/(runtime*hTubes[t][tp]->GetBinWidth(1)));
		}
		if(hTubesRaw[t])
			hTubesRaw[t]->Scale(1.0/(runtime*hTubesRaw[t]->GetBinWidth(1)));
	}
	
	// perform fits, draw histograms
	float searchsigma;
	for(unsigned int t=0; t<=nBetaTubes; t++) {
		
		if(expectedPeaks.size()) {
			if(PCal)
				searchsigma = 0.5*PCal->energyResolution(mySource.mySide, t, expectedPeaks[0].energy(), mySource.x, mySource.y, 0);
			else
				searchsigma = 0.5*sqrt((t==nBetaTubes?2.5:10)*expectedPeaks[0].energy());
			if(!(searchsigma==searchsigma)) {
				printf("Bad search range! Aborting!\n");
				continue;
			}
			if(searchsigma > 120)
				searchsigma = 120;
			if(searchsigma < 4)
				searchsigma = 4;		
			
			OM->defaultCanvas->SetLogy(false);
			std::string fitPlotName = "";
			if(dbgplots)
				fitPlotName = OM->plotPath+"/"+mySource.name()+"/Fit_Spectrum_"+(t==nBetaTubes?"Combined":itos(t))+".pdf";
			printf("Fitting %s for %i peaks in tube %i (sigma = %f)\n", mySource.name().c_str(), (int)expectedPeaks.size(), t, searchsigma);
			tubePeaks[t] = fancyMultiFit(hTubes[t][TYPE_0_EVENT], searchsigma, expectedPeaks, false, fitPlotName, nSigma, pkMin);
			for(std::vector<SpectrumPeak>::iterator it = tubePeaks[t].begin(); it != tubePeaks[t].end(); it++) {
				it->simulated = simMode;
				it->t = t;
				if(PCal) {
					it->eta = PCal->eta(mySource.mySide,t,mySource.x,mySource.y);
					it->nPE = PCal->nPE(mySource.mySide,t,it->energyCenter.x,mySource.x,mySource.y,0);
				}
			}
			
			if(tubePeaks[t].size()<expectedPeaks.size()) {
				Stringmap m;
				m.insert("peak",expectedPeaks.back().name());
				m.insert("sID",mySource.sID);
				m.insert("tube",t);
				m.insert("simulated",simMode?"yes":"no");
				m.insert("side",sideSubst("%c",mySource.mySide));
				OM->warn(MODERATE_WARNING,"Missing_Peak",m);
				printf("Cancelling fit.\n");
				continue;
			} else if(t==nBetaTubes) {
				for(unsigned int i=0; i<tubePeaks[t].size(); i++) {
					printf("-------- %c%i --------\n",sideNames(mySource.mySide),t);
					tubePeaks[t][i].toStringmap().display();
					OM->qOut.insert("Main_"+tubePeaks[t][i].name()+"_peak_"+itos(t),tubePeaks[t][i].toStringmap());
					SourceDBSQL::getSourceDBSQL()->addPeak(tubePeaks[t][i]);
				}
			}
		}
		
		if(t==nBetaTubes)
			continue;
		
		// individual tube plots, convert peaks back to PMT ADC values
		OM->defaultCanvas->SetLogy(true);
		if(!simMode)
			hTubesRaw[t]->Draw();
		for(unsigned int i=0; i<tubePeaks[t].size(); i++) {
			if(PCal && !simMode) {
				float_err fx = tubePeaks[t][i].energyCenter;	// center in energy --- ??????
				tubePeaks[t][i].center = PCal->invertCorrections(mySource.mySide, t, fx, mySource.x, mySource.y, 0.0);
				fx.err = tubePeaks[t][i].energyWidth.x;
				tubePeaks[t][i].width.x = PCal->invertCorrections(mySource.mySide, t, fx, mySource.x, mySource.y, 0.0).err;
				tubePeaks[t][i].gms = PCal->gmsFactor(mySource.mySide, t, 0);
				drawVLine(tubePeaks[t][i].center.x*tubePeaks[t][i].gms, OM->defaultCanvas,2);
				drawVLine(PCal->invertCorrections(mySource.mySide, t, fx.x-fx.err, mySource.x, mySource.y, 0.0)*tubePeaks[t][i].gms, OM->defaultCanvas, 3);
				drawVLine(PCal->invertCorrections(mySource.mySide, t, fx.x+fx.err, mySource.x, mySource.y, 0.0)*tubePeaks[t][i].gms, OM->defaultCanvas, 3);
				drawVLine((tubePeaks[t][i].center.x-tubePeaks[t][i].width.x)*tubePeaks[t][i].gms, OM->defaultCanvas,4);
				drawVLine((tubePeaks[t][i].center.x+tubePeaks[t][i].width.x)*tubePeaks[t][i].gms, OM->defaultCanvas,4);				
			}
			std::string qoutname = "Main_"+tubePeaks[t][i].name()+"_peak_"+itos(t);
			printf("-------- %c%i %s --------\n",sideNames(mySource.mySide),t,qoutname.c_str());
			tubePeaks[t][i].toStringmap().display();
			OM->qOut.insert(qoutname,tubePeaks[t][i].toStringmap());
			SourceDBSQL::getSourceDBSQL()->addPeak(tubePeaks[t][i]);
		}
		if(dbgplots && !simMode)
			OM->printCanvas(mySource.name()+"/Raw_ADC_"+itos(t)+(simMode?"_Sim":""));
	}
	
	printf("Completed fitting procedure.\n");
	
	// combined energy plots
	if(dbgplots) {
		for(unsigned int tp=TYPE_0_EVENT; tp<=TYPE_0_EVENT; tp++) {
			std::vector<TH1*> hToPlot;
			for(unsigned int t=0; t<=nBetaTubes; t++) {
				hTubes[t][tp]->SetLineColor(t==nBetaTubes?2:3+t);
				if(hTubes[t][tp]->GetRMS() > 5)
					hToPlot.push_back(hTubes[t][tp]);
			}
			drawSimulHistos(hToPlot);
			OM->printCanvas(mySource.name()+"/Spectrum_Combined"+(tp==TYPE_0_EVENT?"":"_type_"+itos(tp))+(simMode?"_Sim":""));
		}
	}
	
	OM->defaultCanvas->SetLogy(false);
	std::vector<TH1*> hToPlot;
	for(AxisDirection d = X_DIRECTION; d <= Y_DIRECTION; ++d)
		hToPlot.push_back(hitPos[d]);
	drawSimulHistos(hToPlot);
	OM->printCanvas(mySource.name()+"/Hit_Positions"+(simMode?"_Sim":""));
}

void ReSourcer::calcPMTcorr() {
	printf("Calculating PMT correlations...\n");
	assert(corrFitE0.size()==corrFitE1.size());
	for(unsigned int i=0; i<corrFitE0.size(); i++) {
		TF1 fLine("fLine","pol1",corrFitE0[i],corrFitE1[i]);
		fLine.SetLineColor(simMode?4:2);
		for(unsigned int t1 = 0; t1 < nBetaTubes; t1++) {
			for(unsigned int t2 = 0; t2 < nBetaTubes; t2++) {
				if(t1==t2) continue;
				pPMTCorr[t1][t2]->Fit(&fLine,"QR+");
				Stringmap m;
				m.insert("t1",t1);
				m.insert("t2",t2);
				m.insert("corr",fLine.GetParameter(1));
				m.insert("dcorr",fLine.GetParError(1));
				m.insert("E0",corrFitE0[i]);
				m.insert("E1",corrFitE1[i]);
				m.insert("sID",mySource.sID);
				m.insert("name",mySource.name());
				m.insert("simulated",simMode?"yes":"no");			
				m.insert("side",sideSubst("%c",mySource.mySide));
				m.display();
				OM->qOut.insert("correlation",m);
			}
		}
	}
}

void reSource(RunNum rn) {
	
	// load data
	PostOfficialAnalyzer* P = new PostOfficialAnalyzer(true);
	P->addRun(rn);
	
	if(!P->getnFiles()) {
		SMExcept e("MissingProcessedRun");
		e.insert("runnum",rn);
		delete(P);
		throw(e);
	}
	if(P->totalTime[BOTH] < 1.0) { //TODO
		SMExcept e("RunTooShort");
		e.insert("runnum",rn);
		e.insert("runtime",P->totalTime[BOTH]);
		delete(P);
		throw(e);
	}	
	
	// set up output paths
	std::string outPath = getEnvSafe("UCNA_ANA_PLOTS")+"/LivermoreSources/";
	PMTCalibrator PCal(rn);
	RunInfo RI = CalDBSQL::getCDB()->getRunInfo(rn);
	OutputManager TM("Run_"+itos(RI.runNum), outPath+replace(RI.groupName,' ','_')+"/"+itos(rn)+"_"+RI.roleName+"/");
	
	// get sources list; set up ReSourcer for each
	std::map<unsigned int,ReSourcer> sources;
	std::vector<Source> expectedSources = SourceDBSQL::getSourceDBSQL()->runSources(rn);
	for(std::vector<Source>::iterator it = expectedSources.begin(); it != expectedSources.end(); it++) {
		sources.insert(std::make_pair(it->sID,ReSourcer(&TM,*it,&PCal)));
		Stringmap m = it->toStringmap();
		for(unsigned int t=0; t<nBetaTubes; t++)
			m.insert("eta_"+itos(t),PCal.eta(it->mySide,t,it->x,it->y));
		TM.qOut.insert("Source",m);
	}
	
	// save GMS data to output file
	for(Side s = EAST; s <= WEST; ++s) {
		Stringmap m;
		for(unsigned int t=0; t<nBetaTubes; t++) {
			m.insert("gms0_"+itos(t),PCal.getGMS0(s,t));
			m.insert("gmsRel_"+itos(t),PCal.gmsFactor(s,t,0)/PCal.getGMS0(s,t));
		}
		TM.qOut.insert("BetaSc"+sideNames(s),m);
	}
	
	printf("Expecting %i sources...\n",(int)sources.size());
	
	// all positions histogram
	TH2F* hitPos[2];
	for(Side s = EAST; s <= WEST; ++s)
		hitPos[s] = TM.registeredTH2F(sideSubst("HitPos_%c",s),sideSubst("%s Hit Positions",s),400,-65,65,400,-65,65);
	
	// collect source data points
	P->startScan();
	unsigned int nSPts = 0;
	while(P->nextPoint()) {
		Side s = P->fSide;
		if(P->fType <= TYPE_III_EVENT && (P->fPID == PID_BETA || P->fPID == PID_MUON) && (s==EAST || s==WEST)) {
			P->recalibrateEnergy();
			hitPos[s]->Fill(P->wires[P->fSide][X_DIRECTION].center,P->wires[P->fSide][Y_DIRECTION].center);
			for(std::map<unsigned int, ReSourcer>::iterator it = sources.begin(); it != sources.end(); it++)
				nSPts += it->second.fill(*P);
		}
	}
	
	printf("Located %i source hits...\n",nSPts);
	
	// plot hit positions
	TM.defaultCanvas->SetLogy(false);
	for(Side s = EAST; s <= WEST; ++s) {
		hitPos[s]->Draw("Col");
		for(std::vector<Source>::const_iterator it = expectedSources.begin(); it != expectedSources.end(); it++)
			if(it->mySide==s)
				drawEllipseCut(*it,4.0,it->name());
		TM.printCanvas(sideSubst("HitPos_%c",s));
	}
	
	// fit, graph, output
	for(std::map<unsigned int, ReSourcer>::iterator it = sources.begin(); it != sources.end(); it++) {
		
		printf("\n------------------------------\n");
		printf("Source %i:'%s' [%i events]\n",it->first,it->second.mySource.name().c_str(),(int)it->second.counts());
		it->second.mySource.display("\t");
		printf("------------------------------\n");
		if(it->second.counts()<1000) {
			printf("\tToo few source events for consideration.\n");
			continue;
		}
		
		// delete previous source fit data from DB
		SourceDBSQL::getSourceDBSQL()->clearPeaks(it->second.mySource.sID);
		
		// fit source peaks
		it->second.dbgplots = PCal.isRefRun() || it->second.mySource.t=="Bi207" || it->second.mySource.t=="Cd109";
		it->second.simMode = false;
		it->second.findSourcePeaks(P->totalTime[BOTH]);
		it->second.calcPMTcorr();
		
		// fit simulated source data with same parameters
		Source simSource = it->second.mySource;	
		Sim2PMT* g2p = NULL;
		std::string g4dat = "/home/mmendenhall/geant4/output/FixGeom_";
		printf("Loading source simulation data...\n");
		if(simSource.t=="Bi207" || simSource.t=="Ce139" || simSource.t=="Sn113")
			g4dat = "/home/mmendenhall/geant4/output/20120823_";
		if(simSource.t=="Ce139" || simSource.t=="Sn113" || simSource.t=="Bi207" ||
		   simSource.t=="Cd109" || simSource.t=="In114E" || simSource.t=="In114W") {
			g2p = new G4toPMT();
			g2p->addFile(g4dat + simSource.t + "/analyzed_*.root");
		} else if(simSource.t=="Cd113m") {
			/*
			 G4toPMT* cd109 = new G4toPMT();
			 cd109->addFile(g4dat+"Cd109_geomC/analyzed_*.root");
			 G4toPMT* cd113m = new G4toPMT();
			 cd113m->addFile(g4dat+"Cd113m_geomC/analyzed_*.root");
			 MixSim* MS = new MixSim();
			 MS->addSim(cd113m, 1.0, 14.1*365*24*3600);
			 MS->addSim(cd109, 0.25, 461.4*24*3600);
			 g2p = MS;
			 */
		}
		if(!g2p || !g2p->getnFiles()) {
			printf("Unknown source '%s'!\n",simSource.t.c_str());
			continue;
		}
		
		printf("Preparing to simulate source data...\n");
		
		SourcedropPositioner SDP(simSource.x, simSource.y,
								 (simSource.t=="Ce139" || simSource.t=="Cd109")?1.25:1.5);
		g2p->SP = &SDP;
		g2p->setCalibrator(PCal);
		g2p->fakeClip = true;
		
		ReSourcer RS(&TM,simSource,&PCal);
		RS.simMode = true;
		RS.dbgplots = false;
		
		unsigned int nSimmed = 0;
		g2p->startScan(true);
		while(nSimmed<100000) {
			g2p->nextPoint();
			nSimmed+=RS.fill(*g2p);
		}
		delete(g2p);
		RS.findSourcePeaks(1000.0);
		RS.calcPMTcorr();
		
		// plot data and MC together
		TM.defaultCanvas->SetLeftMargin(1.25);
		for(unsigned int t=0; t<=nBetaTubes; t++) {
			float simNorm = it->second.hTubes[t][TYPE_0_EVENT]->Integral()/RS.hTubes[t][TYPE_0_EVENT]->Integral();
			for(unsigned int tp=TYPE_0_EVENT; tp<=TYPE_III_EVENT; tp++) {
				if(tp>TYPE_0_EVENT && t!=nBetaTubes) continue;
				std::vector<TH1*> hToPlot;
				RS.hTubes[t][tp]->Scale(simNorm);
				
				RS.hTubes[t][tp]->SetLineColor(4);
				it->second.hTubes[t][tp]->SetLineColor(2);
				//RS.hTubes[t][tp]->GetSumw2()->Set(0);
				
				hToPlot.push_back(RS.hTubes[t][tp]);
				hToPlot.push_back(it->second.hTubes[t][tp]);
				TM.defaultCanvas->SetLogy(true);
				drawSimulHistos(hToPlot);
				TM.printCanvas(it->second.mySource.name()+"/Spectrum_Comparison_"+itos(t)+(tp==TYPE_0_EVENT?"":"_type_"+itos(tp)));
				
				// same, linear scale
				TM.defaultCanvas->SetLogy(false);
				RS.hTubes[t][tp]->SetMinimum(0);
				it->second.hTubes[t][tp]->SetMinimum(0);
				drawSimulHistos(hToPlot);
				std::string outName = it->second.mySource.name()+"/Spectrum_Comparison_Lin_"+itos(t)+(tp==TYPE_0_EVENT?"":"_type_"+itos(tp));
				TM.printCanvas(outName);
			}
			
			if(0) {
				TM.defaultCanvas->SetLogy(false);
				for(unsigned int t2=0; t2<nBetaTubes; t2++) {
					if(t==t2 || t==nBetaTubes) continue;
					RS.pPMTCorr[t][t2]->SetLineColor(4);
					RS.pPMTCorr[t][t2]->Draw();
					it->second.pPMTCorr[t][t2]->SetLineColor(2);
					it->second.pPMTCorr[t][t2]->Draw("Same");
					TM.printCanvas(it->second.mySource.name()+"/PMTCorr/"+itos(t)+"_v_"+itos(t2));
				}
			}
		}
		
		
	}
	
	TM.write();	
	delete(P);
}

void uploadRunSources() {
	
	std::vector<std::string> sources[2];
	std::string l;
	
	printf("Loading run log...\n");
	std::ifstream fin("Aux/UCNA Run Log.txt");
	
	while (fin.good()) {
		
		std::getline(fin,l);
		if(!l.size()) continue;
		std::vector<std::string> words = split(l);
		if(!words.size()) continue;
		
		if(l[0]=='*') {
			if(words.size() < 2)
				continue;
			RunNum rn;
			if(!sscanf(words[0].c_str(),"*%i",&rn) || words[1] != "SourcesCal")
				continue;
			
			bool needsUpdate = false;
			for (Side s = EAST; s <= WEST; ++s) {
				std::vector<Source> expectedSources =  SourceDBSQL::getSourceDBSQL()->runSources(rn,s);
				if(expectedSources.size() != sources[s].size()) { needsUpdate = true; break; }
				for(unsigned int n=0; n<expectedSources.size(); n++) {
					if(expectedSources[n].t != sources[s][n])
						needsUpdate = true;
				}
			}
			if(needsUpdate) {
				SourceDBSQL::getSourceDBSQL()->clearSources(rn);
				printf("Run %i: Loading %i,%i sources\n",rn,(int)sources[EAST].size(),(int)sources[WEST].size());
				for (Side s = EAST; s <= WEST; ++s) {
					unsigned int nsrc=0;
					for(std::vector<std::string>::iterator it = sources[s].begin(); it != sources[s].end(); it++) {
						Source src(*it,s);
						src.myRun = rn;
						src.x = nsrc++;
						src.display();
						SourceDBSQL::getSourceDBSQL()->addSource(src);
					}
				}
			} else {
				printf("Run %i: sources already loaded.\n",rn);
			}
			
		} else if(words[0]=="@sources") {
			
			bool notSide[2];
			notSide[EAST] = notSide[WEST] = false;
			if(words.size() >= 2) {
				if(words[1]=="E")
					notSide[WEST] = true;
				else if(words[1]=="W")
					notSide[EAST] = true;
			}
			for(Side s = EAST; s<=WEST; ++s) {
				if(notSide[s])
					continue;
				sources[s].clear();
				for(std::vector<std::string>::const_iterator it = words.begin(); it != words.end(); it++) {
					if(*it == "Sn")
						sources[s].push_back("Sn113");
					else if(*it == "Bi")
						sources[s].push_back("Bi207");
					else if(*it == "Cd")
						sources[s].push_back("Cd109");
					else if(*it == "Cd113m")
						sources[s].push_back("Cd113m");
					else if(*it == "Sr85")
						sources[s].push_back("Sr85");
					else if(*it == "Sr90")
						sources[s].push_back("Sr90");
					else if(*it == "In")
						sources[s].push_back("In114");
					else if(*it == "InE")
						sources[s].push_back("In114E");
					else if(*it == "InW")
						sources[s].push_back("In114W");
					else if(*it == "Ce")
						sources[s].push_back("Ce139");
					else if(*it == "Cs137")
						sources[s].push_back("Ce137");
				}
			}
		}
	}
	fin.close();
}


