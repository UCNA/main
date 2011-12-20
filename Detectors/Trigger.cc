#include "Trigger.hh"
#include <map>
#include "GraphicsUtils.hh"
#include "PathUtils.hh"
#include "ManualInfo.hh"

void Trigger::specialize() {
	
	loadData("Sis00","","Sis00");		// 0	trigType
	loadData("S83028","","runclock");	// 1	runclock
	loadData("S8200","","beamclock");	// 2	protonClock
	loadData("Tdc016","","2of4_tdc");	// 3	E 2of4
	loadData("Tdc017","");				// 4	W 2of4
	loadData("S8300","");				// 5	E Beta Counter
	loadData("S8301","");				// 6	W Beta Counter
	loadData("S8302","");				// 7	UCN Mon 1 trigger counter
	loadData("S8303","");				// 8	UCN Mon 2 trigger counter
	loadData("S8304","");				// 9	UCN Mon 3 trigger counter
	
	if(thisRun->checkBranch("Delt0"))
		loadData("Delt0","","Delt0");		// 10	high resolution delta-t
	else {
		warn(SEVERE_WARNING,"Missing_Delt0");
		printf("...Proceeding without Delt0...\n");
	}
	if(!(thisRun->checkBranch("Clk0") && thisRun->checkBranch("Clk1"))) {
		warn(SEVERE_WARNING,"Missing_Blinding_Clocks");
		loadData("S83028","");				// 11   E clock UNBLINDED
		loadData("S83028","");				// 12   W clock UNBLINDED
	} else {
		loadData("Clk0","","clkE");			// 11   E clock
		loadData("Clk1","","clkW");			// 12   W clock
	}
	loadData("S8308","","gated_scaler");	// 13	E Gated trigger counter
	loadData("S8309","");					// 14	W Gated trigger counter
	
	loadData("Tdc02","TDCE1","tube_tdc");	// 15	E1 TDC
	loadData("Tdc03","TDCE2");				// 16	E2 TDC
	loadData("Tdc00","TDCE3");				// 17	E3 TDC
	loadData("Tdc01","TDCE4");				// 18	E4 TDC
	loadData("Tdc08","TDCW1");				// 19	W1 TDC
	loadData("Tdc09","TDCW2");				// 20	W2 TDC
	loadData("Tdc014","TDCW3");				// 21	W3 TDC
	loadData("Tdc011","TDCW4");				// 22	W4 TDC
}

Trigger::Trigger(RunManager* T): Subsystem(T,"Trigger",BOTH) {
	
		
	specialize();	
	fixTimes();
	// fill events
	float* f_rclock = getData("runclock");
	float* f_bclock = getData("runclock");
	float* f_delt0 = NULL;
	if(thisRun->checkBranch("Delt0"))
		f_delt0 = getData("Delt0");

	for(unsigned int e=0; e<nEvents; e++) {
		anEvent.runClock = f_rclock[e];
		anEvent.beamClock = f_bclock[e];
		anEvent.runNum = thisRun->RI.runNum;
		anEvent.e = e;
		anEvent.trigflags = getTrigFlags(e);
		if(f_delt0)
			anEvent.dt = f_delt0[e];
		else 
			anEvent.dt = 0;
		events.push_back(anEvent);
	}
	
	// set up event vetos
	for(unsigned int i=0; i<nEvents; i++)
		triggers[i] = true;
	if(thisRun->RI.type == ASYMMETRY_RUN) {
		examineBeamStructure();
	}
	
	calcRuntime();
	
	// determine timing cuts
	for(Side s=EAST; s<=WEST; s=nextSide(s)) {
		TH1F* ht1 = registeredTH1F(std::string("Timing_")+sideNames(s),"Scintillator Timing",300,10,4000);
		TH1F* ht2 = registeredTH1F(std::string("Oneside_Timing_")+sideNames(s),"Scintillator Timing",300,10,4000);
		ht2->SetLineColor(2);
		float* f_2tdc = getData("2of4_tdc",s);
		for(UInt_t e=0; e<nEvents; e++) {
			if(isCrud(e) || !isDataEvent(e))
				continue;
			ht1->Fill(f_2tdc[e]);
			if(!beta2of4(e,otherSide(s)))
				ht2->Fill(f_2tdc[e]);
		}
		scintTDCcuts[s] = ht2->GetMean()-4.0*ht2->GetRMS();
		printf("TDC timing cut %i: mu=%f, rms=%f, cut=%f (%i entries)\n",s,ht2->GetMean(),ht2->GetRMS(),scintTDCcuts[s],(int)ht1->GetEntries());
		
		if(thisRun->runMode & RUNMODE_MAKEPLOTS) {
			defaultCanvas->SetLogy(true);
			ht1->Draw();
			ht2->Draw("Same");
			drawVLine(scintTDCcuts[s],defaultCanvas);
			printCanvas(std::string("Scint_Timing_")+ctos(sideNames(s)));
			defaultCanvas->SetLogy(false);
		}
	}
	
	checkDeadtime();
	if(thisRun->runMode & RUNMODE_MAKEPLOTS)
		makeHistos();
	printf("%.1f minutes cut run time, %i events\n",runTime()/60.0,nEvents);
}

void Trigger::addOutBranches(TTree* T) { T->Branch(name.c_str(),&anEvent,"runClock/F:beamClock/F:runNum/I:event/I:flags/I:dt/F"); }

void Trigger::fillEvent(UInt_t e) { anEvent = events[e]; }

void Trigger::checkDeadtime() {
	
	std::map<unsigned int, std::map<int,unsigned int> > skips[2];
	std::map<unsigned int, std::map<int,unsigned int> >::iterator it0;
	std::map<int,unsigned int>::iterator it1;
	
	// count skipped events
	// and tabulate by preceding trigger type
	for(Side s = EAST; s <= WEST; s=nextSide(s)) {
		for(unsigned int e=1; e<nEvents; e++) {
			if(isCrud(e) || isCrud(e-1))
				continue;
			float* f_gs = getData("gated_scaler",s);
			int delta = int(f_gs[e]-f_gs[e-1]);
			if(delta < 0)
				skips[s][sis00(e-1)][10000-delta]++;
			if(beta2of4(e,s))
				skips[s][sis00(e-1)][delta-1]++;
			else
				skips[s][sis00(e-1)][int(-delta-1)]++;
		}
	}
	
	// record skipped events to output files
	makePath(thisRun->dataPath+"/Deadtime/");
	std::string skipfile = thisRun->dataPath+"Deadtime/Skips_"+itos(thisRun->RI.runNum)+".txt";
	FILE* skiphistos = fopen(skipfile.c_str(),"w");
	fprintf(skiphistos,"Runtime: %f\n",runTime());
	for(Side s = EAST; s <= WEST; s=nextSide(s)) {
		for(it0 = skips[s].begin(); it0 != skips[s].end(); it0++) {
			unsigned int ntot, nskipped;
			ntot = nskipped = 0;
			for(it1 = it0->second.begin(); it1 != it0->second.end(); it1++) {
				ntot += it1->second;
				if(it1->first > 0)
					nskipped += it1->second;
			}
			fprintf(skiphistos,"Side: %c\tFollows: %i\tEntries: %i\tTotal: %i\tSkipped: %i\n",sideNames(s),it0->first,(int)it0->second.size(),ntot,nskipped);
			for(it1 = it0->second.begin(); it1 != it0->second.end(); it1++)
				fprintf(skiphistos,"\t%i\t%i\n",it1->first,it1->second);
		}
	}
	fclose(skiphistos);
}


void Trigger::makeHistos() {
	
	defaultCanvas->cd();
	defaultCanvas->SetLogy(true);
	
	// force 8s bins for UCN monitors
	float cyclength =  forParent.getDefault("medianBeamPulseSpacing", 8.0);
	int ncycles = int(totalTime()/cyclength)+1;
	float newtime = cyclength*ncycles;
	TH1F hMon1 = TH1F("hMon1","Event Rates",ncycles,0,newtime);
	TH1F hMon2 = TH1F("hMon2","Event Rates",ncycles,0,newtime);
	TH1F hMon3 = TH1F("hMon3","Event Rates",ncycles,0,newtime);
	TH1F hMon4 = TH1F("hMon4","Event Rates",ncycles,0,newtime);
	TH1F hBetaE = TH1F("hBetaE","Event Rates",200,0,newtime);
	TH1F hBetaW = TH1F("hBetaW","Event Rates",200,0,newtime);
	TH1F hBeam = TH1F("hBeam","Event Rates",ncycles,0,newtime);
	TH1F hGMS = TH1F("hGMS","Event Rates",1000,0,newtime);
	
	hMon1.SetLineColor(1);	// black
	hMon2.SetLineColor(2);	// red
	hMon3.SetLineColor(3);	// green
	hMon4.SetLineColor(4);	// blue
	hBetaE.SetLineColor(5);	// yellow
	hBetaW.SetLineColor(6);	// magenta
	hBeam.SetLineColor(7);	// cyan
	hGMS.SetLineColor(28);	// brown
	
	for(unsigned int e=0; e<nEvents; e++) {
		if(isCrud(e)) {
			hBeam.Fill(eventTime(e));
			continue;
		}
		if(gmsTrigger(e) || pulserTrigger(e)) {
			hGMS.Fill(eventTime(e));
			continue;
		}
		if(isMon1(e)) hMon1.Fill(eventTime(e));
		if(isMon2(e)) hMon2.Fill(eventTime(e));
		if(isMon3(e)) hMon3.Fill(eventTime(e));
		if(isMon4(e)) hMon4.Fill(eventTime(e));
		if(beta2of4(e,EAST)) hBetaE.Fill(eventTime(e));
		if(beta2of4(e,WEST)) hBetaW.Fill(eventTime(e));
		
	}
	
	hMon1.Scale(1.0/hMon1.GetBinWidth(1));
	hMon2.Scale(1.0/hMon2.GetBinWidth(1));
	hMon3.Scale(1.0/hMon3.GetBinWidth(1));
	hMon4.Scale(1.0/hMon4.GetBinWidth(1));
	hBetaE.Scale(1.0/hBetaE.GetBinWidth(1));
	hBetaW.Scale(1.0/hBetaW.GetBinWidth(1));
	hBeam.Scale(1.0/hBeam.GetBinWidth(1));
	hGMS.Scale(1.0/hGMS.GetBinWidth(1));
	
	hGMS.SetMinimum(1e-1);
	hGMS.SetMaximum(1e3);
	
	hGMS.Draw();
	hMon1.Draw("Same");
	hBeam.Draw("Same");
	hMon2.Draw("Same");
	hMon3.Draw("Same");
	hMon4.Draw("Same");
	hBetaE.Draw("Same");
	hBetaW.Draw("Same");
	
	for(std::vector<Blip>::const_iterator it = blips.begin(); it != blips.end(); it++)
		drawExcludedRegion(getData("runclock")[it->start], getData("runclock")[it->end], defaultCanvas);
	
	printCanvas("EventRates");
	
	for(Side s = EAST; s <= WEST; s = nextSide(s)) {
		for(unsigned int t=0; t<nBetaTubes; t++) {
			TH1F* hPMTTDC = registeredTH1F(std::string("PMT_TDC_")+sideNames(s)+itos(t),"PMT TDC",400,0,4500);
			hPMTTDC->SetLineColor(t+2);
			float* ptdc = getData("tube_tdc",nBetaTubes*s+t);
			for(unsigned int e=0; e<nEvents; e++)
				if(isDataEvent(e) && !isCrud(e))
					hPMTTDC->Fill(ptdc[e]);
			if(!t)
				hPMTTDC->Draw();
			else
				hPMTTDC->Draw("Same");
		}
		printCanvas(std::string("PMT_TDC_")+sideNames(s));
	}
	
	defaultCanvas->SetLogy(false);
}

void Trigger::fixTimes() {
	
	// fix scaler reset jumps
	float deltaT = 0;
	int njumps = 0;
	float* f_rclock = getData("runclock");
	float* f_clkE = getData("clkE");
	float* f_clkW = getData("clkW");
	for(unsigned int e=1; e<nEvents; e++) {
		if(f_rclock[e] < f_rclock[e-1]-deltaT-1e9) {
			deltaT += 4294967296.0;
			njumps++;
		}
		f_rclock[e] += deltaT;
		f_clkE[e] += deltaT;
		f_clkW[e] += deltaT;
	}
	
	if(njumps)
		printf("Fixed %i timing scaler overflows.\n",njumps);
	
	// convert us to s
	float* f_bclock = getData("beamclock");
	for(unsigned int e=0; e<nEvents; e++) {
		f_rclock[e] *= 1.0e-6;
		f_bclock[e] *= 1.0e-6;
		f_clkE[e] *= 1.0e-6;
		f_clkW[e] *= 1.0e-6;
	}
}


std::vector<Blip> Trigger::makeBlips(const bool* trg) const {
	std::vector<Blip> B;
	unsigned int e0=0;
	const float* f_rclock = getData("runclock");
	for(unsigned int e=1; e<nEvents; e++) {
		if(trg[e] && !trg[e-1]) {
			Blip b;
			b.start = e0;
			b.end = e;
			b.time = f_rclock[e]-f_rclock[e0];
			B.push_back(b);
		}
		else if(!trg[e] && trg[e-1])
			e0 = e;
	}
	if(!trg[nEvents-1]) {
		Blip b;
		b.start = e0;
		b.end = nEvents-1;
		b.time = f_rclock[nEvents-1]-f_rclock[e0];
		B.push_back(b);		
	}
	
	return B;
}

void Trigger::calcRuntime() {
	
	const char* scols[] = { "clkE","clkW","runclock" };	// clock columns for each side
	
	blips = makeBlips(triggers);
	
	for(Side s = EAST; s <= BOTH; s = nextSide(s)) {
		float* f_rclock = getData(scols[s]);
		float ttotal = f_rclock[nEvents-1] - f_rclock[0];
		float tlost = 0;
		for(std::vector<Blip>::const_iterator b = blips.begin(); b != blips.end(); b++)
			tlost += f_rclock[b->end] - f_rclock[b->start];
		runtimes[s] = ttotal-tlost;
		
		if(s==BOTH) {
			printf("Uncut runtime %g minutes; lost %g%% to %i blips\n",ttotal/60.0,100.0*tlost/ttotal,(int)blips.size());
			forParent.insert("lostTime",tlost);
			forParent.insert("dataBlips",blips.size());
		}
	}
	runtimes[NONE] = runtimes[BOTH];	
	
}


// comparison to sort blips by number of events
bool sortBlipNEvents(Blip a, Blip b) {
	return a.counts() < b.counts();
}

void Trigger::examineBeamStructure() {
	
	const float beamnoisethresh = 0.80;	// proportion allowed below median beam noise
	
	// remove segments during beam pulses
	for(unsigned int e=0; e<nEvents; e++)
		triggers[e] = !(isBeamnoise(e) || isBeamout(e));
	
	
	// remove manually tagged segments
	std::vector< std::pair<double,double> > manualCuts = ManualInfo::MI->getRanges(itos(thisRun->RI.runNum)+"_timecut");
	if(manualCuts.size()) {
		printf("Manually cutting %i ranges...\n",(int)manualCuts.size());
		const float* f_rclock = getData("runclock");
		for(unsigned int e=0; e<nEvents; e++)
			for(std::vector< std::pair<double,double> >::const_iterator it = manualCuts.begin(); it != manualCuts.end(); it++)
				triggers[e] = triggers[e] && !(it->first <= f_rclock[e] && f_rclock[e] <= it->second);
	}
	
	std::vector<Blip> beamBlips = makeBlips(triggers);
	unsigned int nBeamBlips = beamBlips.size();
	printf("Located %i beam pulses...\n",nBeamBlips);
	if(!nBeamBlips) {
		warn(SEVERE_WARNING,"No_Beam_Pulses_Found");
		return;
	}
	if(nBeamBlips<5) {
		warn(MODERATE_WARNING,"Few_Beam_Pulses_Found");
		for(unsigned int i=0; i<nBeamBlips; i++)
			beamBlips[i].display();
		return;
	}
	
	examineClusterTiming();
	
	// determine beam noise rate, cut unusually quiet segments
	std::vector<Blip> blipSorted = beamBlips;
	std::sort(blipSorted.begin(),blipSorted.end(),&sortBlipNEvents);
	unsigned int medianBeamNoise = blipSorted[nBeamBlips/2].counts();
	forParent.insert("medianBeamNoiseCounts",medianBeamNoise);
	printf("Median beam noise: %i counts\n",medianBeamNoise);
	if(medianBeamNoise > 100) {
		for(unsigned int i=0; i<nBeamBlips-1; i++)
			if(beamBlips[i].counts() < beamnoisethresh*medianBeamNoise)
				for(unsigned int e=beamBlips[i].start; e<beamBlips[i+1].start; e++)
					triggers[e] = false;
	} else {
		printf("Beam too quiet for beam noise elimination cuts.\n");
		warn(MODERATE_WARNING,"Beam_Unusually_Quiet");
	}
	
	// determine beam spacing TODO use beamclock, tie into beam out cuts
	std::vector<float> pulseSpaces;
	for(unsigned int i=0; i<nBeamBlips-1; i++)
		pulseSpaces.push_back(getData("runclock")[beamBlips[i+1].start]-getData("runclock")[beamBlips[i].start]);
	std::sort(pulseSpaces.begin(),pulseSpaces.end());
	float medianPulseSpacing = pulseSpaces[nBeamBlips/2];
	printf("Median beam pulse spacing: %g s\n",medianPulseSpacing);
	forParent.insert("medianBeamPulseSpacing",medianPulseSpacing);	
}

void Trigger::examineClusterTiming() {
	
	if(!thisRun->checkBranch("Delt0")) {
		warn(MODERATE_WARNING,"Missing_Delt0_For_Clusters");
		return;
	}
	
	const float trapid = 50.0;	// treshold for "rapid" event rates, in us
	TH1F* hRapid = registeredTH1F("Cluster_Event_Spacing","Cluster Event Spacing [log us]",100,1.0,5.0);
	TH1F* hRapidB = registeredTH1F("Cluster_Beam_Event_Spacing","Cluster Event Spacing [log us]",100,1.0,5.0);
	TH1F* hRapidN = registeredTH1F("Cluster_Normal_Event_Spacing","Cluster Event Spacing [log us]",100,1.0,5.0);
	hRapidB->SetLineColor(2);
	hRapidN->SetLineColor(4);
	float* f_delt0 = getData("Delt0");
	for(unsigned int e=1; e<nEvents-1; e++) {
		float maxdel = std::max(f_delt0[e],f_delt0[e+1]);
		hRapid->Fill(log(maxdel)/log(10.0));
		if(!triggers[e])
			hRapidB->Fill(log(maxdel)/log(10.0));
		else
			hRapidN->Fill(log(maxdel)/log(10.0));
		if(maxdel<trapid)
			triggers[e-1] = triggers[e] = triggers[e+1] = false;
	}
	defaultCanvas->SetLogy(true);
	hRapid->Draw();
	hRapidB->Draw("Same");
	hRapidN->Draw("Same");
	drawVLine(log(trapid)/log(10.0), defaultCanvas);
	printCanvas("ClusterEvents");
	defaultCanvas->SetLogy(false);
}
