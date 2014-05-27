#include "ReSource.hh"
#include "PostOfficialAnalyzer.hh"
#include "PMTGenerator.hh"
#include "TSpectrumUtils.hh"
#include "GraphicsUtils.hh"
#include "G4toPMT.hh"
#include "PenelopeToPMT.hh"
#include "PathUtils.hh"
#include "SourceDBSQL.hh"
#include "CalDBSQL.hh"
#include "SMExcept.hh"
#include <TSpectrum.h>
#include <TSpectrum2.h>
#include <utility>

SourceHitsPlugin::SourceHitsPlugin(RunAccumulator* RA, const Source& s, PMTCalibrator* P):
AnalyzerPlugin(RA,"src_"+s.name()), mySource(s), PCal(P), dbgplots(false),
nBins(300), eMin(-100), eMax(2000), pkMin(0.0), pkMax(FLT_MAX), nSigma(2.0)  {
	
	// source-dependent ranges
	double wRange = 10; // position plot +/-width
	unsigned int wBins = 300;
	if(mySource.t == "Bi207") {
		eMax = 1800;
		pkMin = 200;
		nSigma = 1.5;
		wBins = 100;
	} else if(mySource.t == "Sn113") {
		nBins = 150;
		eMin = -50;
		eMax = 1000;
	} else if(mySource.t == "Ce139") {
		nBins = 100;
		eMin = -20;
		eMax = 300;
		wRange = 5;
		wBins = 100;
	} else if (mySource.t == "Cd109") {
		nBins = 100;
		eMin = -20;
		eMax = 300;
		wRange = 5;
		wBins = 100;
		nSigma = 1.0;
	} else if(mySource.t == "Cs137") {
		eMax = 1200;
		pkMin = 400;
		nSigma = 1.0;
	} else if(mySource.t == "In114E" || mySource.t == "In114W") {
		pkMax = 300;
	}
	
	
	// set up histograms
	hErec = myA->registerFGBGPair(s.name()+"_Erec",mySource.lxname()+" energy spectrum",nBins,0,eMax);
	hErec->setAxisTitle(X_DIRECTION,"Energy [keV]");
	for(unsigned int t=0; t<=nBetaTubes; t++) {
		for(unsigned int tp=TYPE_0_EVENT; tp<=TYPE_III_EVENT; tp++) {
			hTubes[t][tp] = myA->registerFGBGPair(s.name()+(t==nBetaTubes?"_Combined":"_Tube_"+itos(t))+"_type_"+itos(tp),mySource.lxname()+" energy spectrum",
												  tp==TYPE_0_EVENT?nBins:nBins/4,eMin,eMax);
			hTubes[t][tp]->setAxisTitle(X_DIRECTION,"Scintillator Visible Energy [keV]");
			hTubes[t][tp]->setAxisTitle(Y_DIRECTION,"Event Rate [Hz/keV]");		}
		if(t<nBetaTubes) {
			hTubesRaw[t] = myA->registerFGBGPair(s.name()+"_"+itos(t)+"_ADC",mySource.lxname()+" ADC spectrum",nBins,
												 PCal?-PCal->invertCorrections(mySource.mySide, t, -eMin, mySource.x, mySource.y, 0.0):2*eMin,
												 PCal?PCal->invertCorrections(mySource.mySide, t, eMax, mySource.x, mySource.y, 0.0):2*eMax);
			hTubesRaw[t]->setAxisTitle(X_DIRECTION,"ADC Channels");
		}
	}
	for(AxisDirection d = X_DIRECTION; d <= Y_DIRECTION; ++d) {
		hitPos[d] = myA->registerFGBGPair(s.name()+"_hits_profile_"+itos(d),mySource.lxname()+" event positions",wBins,-wRange,wRange);
		// TODO hitPos[d]->SetLineColor(2+2*d);
		hitPos[d]->setAxisTitle(X_DIRECTION,(std::string(d==X_DIRECTION?"x":"y")+" position [mm]").c_str());
	}
	
}

void SourceHitsPlugin::fillCoreHists(ProcessedDataScanner& PDS, double weight) {
	EventType tp = PDS.fType;
	Side s = PDS.fSide;
	
	if(tp > TYPE_III_EVENT || s != mySource.mySide) return;
	float x = PDS.wires[s][X_DIRECTION].center;
	float y = PDS.wires[s][Y_DIRECTION].center;
	if(!mySource.inSourceRegion(x,y,4.0)) return;
	if(PDS.fPID != PID_BETA) return;
	
	if(tp==TYPE_0_EVENT) {
		hitPos[X_DIRECTION]->h[currentGV]->Fill(x-mySource.x,weight);
		hitPos[Y_DIRECTION]->h[currentGV]->Fill(y-mySource.y,weight);
	}
	for(unsigned int t=0; t<nBetaTubes; t++) {
		if(PDS.scints[s].adc[t]>3750 || PDS.scints[s].tuben[t].x < 1.0)
			continue;
		hTubes[t][tp]->h[currentGV]->Fill(PDS.scints[s].tuben[t].x,weight);
		if(tp == TYPE_0_EVENT)
			hTubesRaw[t]->h[currentGV]->Fill(PDS.scints[s].adc[t]*PDS.ActiveCal->gmsFactor(mySource.mySide,t,PDS.runClock[BOTH]),weight);
	}
	hErec->h[currentGV]->Fill(PDS.getErecon(),weight);
	if(tp == TYPE_0_EVENT)
		hTubes[nBetaTubes][tp]->h[currentGV]->Fill(PDS.scints[s].energy.x,weight);
	else
		hTubes[nBetaTubes][tp]->h[currentGV]->Fill(PDS.getEnergy(),weight);
}

void SourceHitsPlugin::calculateResults() {
	
	float runtime = myA->getTotalTime(GV_OPEN)[mySource.mySide];
	mySource.nCounts = counts();
	
	// normalize to rates
	for(unsigned int t=0; t<=nBetaTubes; t++) {
		for(unsigned int tp=TYPE_0_EVENT; tp<=TYPE_III_EVENT; tp++) {
			hTubesR[t][tp] = myA->rateHisto(hTubes[t][tp]);
			if(tp != TYPE_0_EVENT) hTubesR[t][tp]->Scale(1000);
			hTubesR[t][tp]->GetYaxis()->SetTitle(tp==TYPE_0_EVENT?"Event Rate [Hz/keV]":"Event Rate [mHz/keV]");
			hTubesR[t][tp]->GetYaxis()->SetTitleOffset(1.25);
			hTubesR[t][tp]->GetYaxis()->SetNdivisions(507);
		}
		if(t<nBetaTubes) {
			hTubesRawR[t] = myA->rateHisto(hTubesRaw[t]);
			hTubesRawR[t]->GetYaxis()->SetTitle("Event Rate [Hz/keV]");
		}
	}
	hErecR = myA->rateHisto(hErec);
	hErecR->GetYaxis()->SetTitleOffset(1.3);
	hErecR->GetYaxis()->SetTitle("Event Rate [Hz/keV]");
	for(AxisDirection d = X_DIRECTION; d <= Y_DIRECTION; ++d) {
		hitPosR[d] = myA->rateHisto(hitPos[d]);
		hitPosR[d]->GetYaxis()->SetTitle("Event Rate [Hz/mm]");
		hitPosR[d]->GetYaxis()->SetTitleOffset(1.25);
		hitPosR[d]->SetMinimum(0);
	}
	
	for(EventType tp = TYPE_0_EVENT; tp <= TYPE_I_EVENT; ++tp) {
		
		// get expected peaks list
		std::vector<SpectrumPeak> expectedPeaks = mySource.getPeaks(tp);
		if(!expectedPeaks.size()) continue;
		
		if(counts(tp)<1000) {
			printf("\tToo few source events for fitting.\n");
			continue;
		}
		
		printf("Fitting peaks for %s source (%i hits, %.1f minutes, %i peaks)...\n",mySource.name().c_str(),counts(tp),runtime/60.0,(int)expectedPeaks.size());
		
		for(unsigned int t=0; t<=nBetaTubes; t++) {
			
			if(tp>TYPE_0_EVENT && t<nBetaTubes) continue;
			
			// estimate for peak width
			float searchsigma;
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
			
			std::string fitPlotName = ( myA->plotPath + "/"
									   + mySource.name() + "/Fit_Spectrum_"
									   + (t==nBetaTubes?"Combined":itos(t))
									   + (tp>TYPE_0_EVENT?"_"+typeWords(tp):"") + ".pdf");
			if(!(dbgplots || tp>TYPE_0_EVENT)) fitPlotName = "";
			printf("Fitting %s for %i peaks in tube %i (sigma = %f)\n", mySource.name().c_str(), (int)expectedPeaks.size(), t, searchsigma);
			double fitsigma = tp>TYPE_0_EVENT ? 1.3 : nSigma;
			tubePeaks[t] = fancyMultiFit(hTubesR[t][tp], searchsigma, expectedPeaks, false, fitPlotName, fitsigma, pkMin, pkMax);
			
			// exit if fit finds too few peaks
			if(tubePeaks[t].size()<expectedPeaks.size()) {
				Stringmap m;
				m.insert("peak",expectedPeaks.back().name());
				m.insert("sID",mySource.sID);
				m.insert("tube",t);
				m.insert("simulated",myA->isSimulated?"yes":"no");
				m.insert("side",sideSubst("%c",mySource.mySide));
				myA->warn(MODERATE_WARNING,"Missing_Peak",m);
				
				if(expectedPeaks.size()==1 && mySource.t != "Cs137") {
					printf("Defaulting to mean/RMS...\n");
					tubePeaks[t] = expectedPeaks;
					tubePeaks[t][0].energyCenter = hTubesR[t][tp]->GetMean();
					tubePeaks[t][0].energyWidth = hTubesR[t][tp]->GetRMS();
				} else {
					printf("Cancelling fit.\n");
					continue;
				}
			}
			
			// display and upload peaks
			for(std::vector<SpectrumPeak>::iterator it = tubePeaks[t].begin(); it != tubePeaks[t].end(); it++) {
				it->simulated = myA->isSimulated;
				it->t = t;
				if(PCal) {
					it->eta = PCal->eta(mySource.mySide, t, mySource.x, mySource.y);
					it->nPE = PCal->nPE(mySource.mySide, t, it->energyCenter.x, mySource.x, mySource.y, 0);
					if(t<nBetaTubes) {
						// individual PMT calibration inverse, from energy to ADC units
						float_err fx = it->energyCenter;
						it->center = PCal->invertCorrections(mySource.mySide, t, fx, mySource.x, mySource.y, 0.0);
						fx.err = it->energyWidth.x;
						it->width.x = PCal->invertCorrections(mySource.mySide, t, fx, mySource.x, mySource.y, 0.0).err;
						it->gms = PCal->gmsFactor(mySource.mySide, t, 0);
					}
				}
				printf("-------- %c%i %s --------\n",sideNames(mySource.mySide),t,it->name().c_str());
				it->toStringmap().display();
				SourceDBSQL::getSourceDBSQL()->addPeak(*it);
			}
		}
	}
	
	printf("Completed fitting procedure.\n");
}

void SourceHitsPlugin::makePlots() {
	
	if(counts()<1000) {
		printf("\tToo few source events for plotting.\n");
		return;
	}
	
	myA->defaultCanvas->SetLeftMargin(2.0);
	myA->defaultCanvas->SetRightMargin(0.05);
	
	// Erecon spectrum
#ifdef PUBLICATION_PLOTS
	hErecR->SetLineColor(1);
	hErecR->SetLineWidth(2);
#endif
	if(hErecR->GetMaximum() < 0.51) hErecR->GetYaxis()->SetRangeUser(0,0.51);
	hErecR->Draw();
	hErecR->Draw("HIST SAME");
	printCanvas(mySource.name()+"/Erecon");
	
	if(!dbgplots) return;
	
	myA->defaultCanvas->SetLogy(true);
	
	for(unsigned int t=0; t<nBetaTubes; t++) {
		
		// raw ADC plots
		if(!myA->isSimulated) {
			hTubesRawR[t]->Draw();
			if(PCal) {
				for(std::vector<SpectrumPeak>::iterator it = tubePeaks[t].begin(); it != tubePeaks[t].end(); it++) {
					float_err fx = it->energyCenter;	// center in energy --- ??????
					drawVLine(it->center.x*it->gms, myA->defaultCanvas,2);
					drawVLine(PCal->invertCorrections(mySource.mySide, t, fx.x-fx.err, mySource.x, mySource.y, 0.0)*it->gms, myA->defaultCanvas, 3);
					drawVLine(PCal->invertCorrections(mySource.mySide, t, fx.x+fx.err, mySource.x, mySource.y, 0.0)*it->gms, myA->defaultCanvas, 3);
					drawVLine((it->center.x-it->width.x)*it->gms, myA->defaultCanvas,4);
					drawVLine((it->center.x+it->width.x)*it->gms, myA->defaultCanvas,4);
				}
			}
			myA->printCanvas(mySource.name()+"/Raw_ADC_"+itos(t)+(myA->isSimulated?"_Sim":""));
		}
		
		// combined energy plots
		for(unsigned int tp=TYPE_0_EVENT; tp<=TYPE_0_EVENT; tp++) {
			std::vector<TH1*> hToPlot;
			for(unsigned int t=0; t<=nBetaTubes; t++) {
				hTubesR[t][tp]->SetLineColor(t==nBetaTubes?2:3+t);
				if(hTubesR[t][tp]->GetRMS() > 5)
					hToPlot.push_back(hTubesR[t][tp]);
			}
			drawSimulHistos(hToPlot);
			myA->printCanvas(mySource.name()+"/Spectrum_Combined"+(tp==TYPE_0_EVENT?"":"_type_"+itos(tp))+(myA->isSimulated?"_Sim":""));
		}
	}
	
	myA->defaultCanvas->SetLogy(false);
}

void SourceHitsPlugin::compareMCtoData(AnalyzerPlugin* AP) {
	SourceHitsPlugin* dat = (SourceHitsPlugin*)AP; // re-cast to correct type
	
	myA->defaultCanvas->SetLeftMargin(2.0);
	myA->defaultCanvas->SetRightMargin(0.05);
	
	// energy reconstruction
	for(unsigned int t=0; t<=nBetaTubes; t++) {
		for(EventType tp=TYPE_0_EVENT; tp<=TYPE_III_EVENT; ++tp) {
			
			if(tp>TYPE_0_EVENT && t!=nBetaTubes) continue;
			
			hTubesR[t][tp]->SetLineColor(4);
			dat->hTubesR[t][tp]->SetLineColor(2);
#ifdef PUBLICATION_PLOTS
			hTubesR[t][tp]->GetSumw2()->Set(0);
#endif
			myA->defaultCanvas->SetLogy(true);
			drawHistoPair(dat->hTubesR[t][tp], hTubesR[t][tp]);
			printCanvas(mySource.name()+"/Spectrum_Comparison_"+itos(t)+(tp==TYPE_0_EVENT?"":"_type_"+itos(tp)));
			// same, linear scale
			myA->defaultCanvas->SetLogy(false);
			hTubesR[t][tp]->SetMinimum(0);
			dat->hTubesR[t][tp]->SetMinimum(0);
			drawHistoPair(dat->hTubesR[t][tp], hTubesR[t][tp]);
			printCanvas(mySource.name()+"/Spectrum_Comparison_Lin_"+itos(t)+(tp==TYPE_0_EVENT?"":"_type_"+itos(tp)));
		}
	}
	
	// positions
	for(AxisDirection d = X_DIRECTION; d <= Y_DIRECTION; ++d) {
#ifdef PUBLICATION_PLOTS
		hitPosR[d]->SetLineStyle(2);
		hitPosR[d]->GetSumw2()->Set(0);
#endif
		drawHistoPair(dat->hitPosR[d],hitPosR[d]);
		printCanvas(mySource.name()+"/Hit_Pos_"+(d==X_DIRECTION?"x":"y"));
	}
}


//-----------------------------------------------------------

SourcePositionsPlugin::SourcePositionsPlugin(RunAccumulator* RA): AnalyzerPlugin(RA,"sourceHitPos") {
	TH2F hTemplate("srcHitPos","Hit Positions",400,-65,65,400,-65,65);
	for(Side s = EAST; s <= WEST; ++s) {
		hitPos[s] = (TH2F*)myA->registerSavedHist(sideSubst("srcHitPos_%c",s),hTemplate);
		hitPos[s]->GetXaxis()->SetTitle("x position [mm]");
		hitPos[s]->GetYaxis()->SetTitle("y position [mm]");
		hitPos[s]->SetTitle(sideSubst("%s Hit Positions",s).c_str());
	}
}

void SourcePositionsPlugin::fillCoreHists(ProcessedDataScanner& PDS, double weight) {
	Side s = PDS.fSide;
	if(PDS.fType <= TYPE_III_EVENT && PDS.fPID == PID_BETA && (s==EAST || s==WEST))
		hitPos[s]->Fill(PDS.wires[s][X_DIRECTION].center,PDS.wires[s][Y_DIRECTION].center,weight);
}

void SourcePositionsPlugin::makePlots() {
	for(Side s = EAST; s <= WEST; ++s) {
		hitPos[s]->Draw("Col");
#ifndef PUBLICATION_PLOTS
		SourceHitsAnalyzer* SHA = (SourceHitsAnalyzer*)myA;
	 	for(std::vector<SourceHitsPlugin*>::const_iterator it = SHA->srcPlugins.begin(); it != SHA->srcPlugins.end(); it++)
			if((*it)->mySource.mySide==s)
				drawEllipseCut((*it)->mySource,4.0,(*it)->mySource.name());
#endif
		printCanvas(sideSubst("HitPos_%c",s));
	}
}

//-----------------------------------------------------------

SourceHitsAnalyzer::SourceHitsAnalyzer(OutputManager* pnt, const std::string& nm, const std::string& inflname):
RunAccumulator(pnt,nm,inflname), PCal(NULL) {
	ignoreMissingHistos = true;
	addPlugin(spos_plgn = new SourcePositionsPlugin(this));
	addPlugin(mwpcgain_plgn = new MWPCGainPlugin(this));
}

SegmentSaver* SourceHitsAnalyzer::makeAnalyzer(const std::string& nm,const std::string& inflname) {
	SourceHitsAnalyzer* SHA = new SourceHitsAnalyzer(this,nm,inflname);
	SHA->PCal = PCal;
	for(std::vector<SourceHitsPlugin*>::const_iterator it = srcPlugins.begin(); it != srcPlugins.end(); it++)
		SHA->addSource((*it)->mySource);
	return SHA;
}

void SourceHitsAnalyzer::addSource(const Source& s) {
	smassert(PCal);
	srcPlugins.push_back(new SourceHitsPlugin(this,s,PCal));
	addPlugin(srcPlugins.back());
	Stringmap m = s.toStringmap();
	for(unsigned int t=0; t<nBetaTubes; t++)
		m.insert("eta_"+itos(t),PCal->eta(s.mySide,t,s.x,s.y));
	qOut.insert("Source",m);
}

//-----------------------------------------------------------

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
	PMTCalibrator PCal(rn);
	RunInfo RI = CalDBSQL::getCDB()->getRunInfo(rn);
	
	OutputManager OMdat("NameUnused", getEnvSafe("UCNA_ANA_PLOTS")+"/SourceFitsData/"+replace(RI.groupName,' ','_')+"/");
	SourceHitsAnalyzer SHAdat(&OMdat,itos(rn)+"_"+RI.roleName);
	SHAdat.grouping = GROUP_RUN;
	SHAdat.PCal = &PCal;
	
	// set up analyzers for each expected source
	std::vector<Source> expectedSources = SourceDBSQL::getSourceDBSQL()->runSources(rn);
	printf("Expecting %i sources...\n",(int)expectedSources.size());
	for(std::vector<Source>::iterator it = expectedSources.begin(); it != expectedSources.end(); it++) {
		SHAdat.addSource(*it);
		// delete previous source fit data from DB
		SourceDBSQL::getSourceDBSQL()->clearPeaks(it->sID);
	}
	
	// collect source data points
	SHAdat.loadProcessedData(AFP_OTHER, GV_OPEN, *P);
	delete(P);
	SHAdat.makeOutput(true);
	
	
	// run simulations for each source
	std::string simOutName = getEnvSafe("UCNA_ANA_PLOTS")+"/SourceFitsSim/"+replace(RI.groupName,' ','_')+"/";
	if(rn < 16300) simOutName = getEnvSafe("UCNA_ANA_PLOTS")+"/SourceFitsSim2010/"+replace(RI.groupName,' ','_')+"/";
	OutputManager OMsim("NameUnused", simOutName);
	SourceHitsAnalyzer SHAsim(&OMsim,itos(rn)+"_"+RI.roleName);
	SHAsim.isSimulated = true;
	SHAsim.PCal = &PCal;
	SHAsim.copyTimes(SHAdat);
	SHAsim.grouping = GROUP_RUN;
	for(std::vector<SourceHitsPlugin*>::iterator it = SHAdat.srcPlugins.begin(); it != SHAdat.srcPlugins.end(); it++) {
		SHAsim.addSource((*it)->mySource);
	}
	for(std::vector<SourceHitsPlugin*>::iterator it = SHAsim.srcPlugins.begin(); it != SHAsim.srcPlugins.end(); it++) {
		
		// load appropriate simulation data
		Source& src = (*it)->mySource;
		Sim2PMT* g2p = NULL;
		printf("Loading source simulation data...\n");
		std::string g4dat = "/data2/mmendenhall/G4Out/2010/FixGeom_";
		if(src.t=="Bi207" || src.t=="Ce139" || src.t=="Sn113") g4dat = "/data2/mmendenhall/G4Out/2010/20120823_";
		if(rn > 16300) g4dat = getEnvSafe("UCNA_CALSRC_SIMS");
		
		if(src.t=="Ce139" || src.t=="Sn113" || src.t=="Bi207" ||
		   src.t=="Cd109" || src.t=="In114E" || src.t=="In114W" || src.t=="Cs137") {
			g2p = new G4toPMT();
			std::string simdat = g4dat + src.t + "/analyzed_*.root";
			printf("Loading data from '%s'\n",simdat.c_str());
			g2p->addFile(simdat);
			g2p->runCathodeSim();
		}
		if(!g2p || !g2p->getnFiles()) {
			printf("Unknown source '%s'!\n",src.t.c_str());
			continue;
		}
		
		printf("\n------------------------------\n");
		printf("Source %i:'%s' [%i events]\n",src.sID,src.name().c_str(),(int)src.nCounts);
		src.display("\t");
		printf("------------------------------\n");
		if(src.nCounts < 1000) {
			printf("\tToo few source events for consideration.\n");
			continue;
		}
		
		printf("Preparing to simulate source data...\n");
		g2p->setCalibrator(PCal);
		g2p->simSide = src.mySide;
		SourcedropPositioner SDP(src.x, src.y, src.t=="Ce139" ? 1.25 : src.t=="Cd109" ? 0.5 : 1.5 );
		g2p->SP = &SDP;
		
		float nRealCounts = 10000; //50000;
		g2p->basePhysWeight = src.nCounts/nRealCounts;
		
		SHAsim.setCurrentState(g2p->getAFP(),GV_OPEN);
		g2p->resetSimCounters();
		g2p->startScan((int)src.nCounts);
		while((*it)->counts() < src.nCounts) {
			g2p->nextPoint();
			SHAsim.loadSimPoint(*g2p);
		}
		printf("\n--Scan complete.--\n");
		
		//SHAsim.loadSimData(*g2p,nCounts);
		delete g2p;
	}
	
	SHAsim.makePlots();
	SHAsim.compareMCtoData(SHAdat);
	SHAsim.uploadAnaResults();
	SHAsim.write();
	SHAsim.setWriteRoot(true);
}





//-------------------------------------------------------------------------------------------------


void uploadRunSources(const std::string& rlogname) {
	
	std::vector<std::string> sources[2];
	std::string l;
	
	printf("Loading run log '%s'...\n",rlogname.c_str());
	std::ifstream fin(rlogname.c_str());
	
	while (fin.good()) {
		
		std::getline(fin,l);
		if(!l.size()) continue;
		std::vector<std::string> words = split(l);
		if(!words.size()) continue;
		
		if(l[0]=='*') {
			if(words.size() < 2)
				continue;
			RunNum rn;
			if(!sscanf(words[0].c_str(),"*%u",&rn) || words[1] != "SourcesCal")
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
					else if(*it == "Cs")
						sources[s].push_back("Cs137");
				}
			}
		}
	}
	fin.close();
}


