#include "Wirechamber.hh"
#include <math.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TF1.h>
#include <TGraph.h>
#include <TLatex.h>
#include "GraphicsUtils.hh"

void Wirechamber::specialize() {
	
	wireSpacing = 10.16*sqrt(0.6);
	setName(sideSubst("Wires_%c",mySide)+((myDirection==X_DIRECTION)?"X":"Y"));
	
	// cathode gain normalization
	float EX_Norm[] = { 252.6, 514.3, 606.0, 579.8, 555.7, 540.5, 536.0, 509.8, 491.3, 477.4, 466.4, 437.2, 436.2, 422.1, 325.6, 218.3 };
	float EY_Norm[] = { 280.8, 309.1, 345.8, 388.2, 426.7, 424.7, 422.1, 435.0, 458.5, 457.2, 483.7, 470.9, 481.7, 431.8, 384.4, 256.3 };
	float WX_Norm[] = { 258.3, 447.1, 572.4, 562.6, 550.0, 537.8, 524.5, 494.1, 485.1, 458.8, 438.0, 432.6, 410.3, 393.4, 321.2, 225.2 };
	float WY_Norm[] = { 280.6, 346.1, 381.5, 404.5, 406.5, 417.5, 411.8, 428.2, 403.2, 415.6, 388.6, 401.1, 418.5, 396.9, 373.2, 307.3 };
	nWires = 16;
	for(unsigned int i=0; i<nWires; i++) {
		if(mySide == EAST && myDirection == X_DIRECTION)
			cathNorm[i] = EX_Norm[i];
		if(mySide == EAST && myDirection == Y_DIRECTION)
			cathNorm[i] = EY_Norm[i];
		if(mySide == WEST && myDirection == X_DIRECTION)
			cathNorm[i] = WX_Norm[i];
		if(mySide == WEST && myDirection == Y_DIRECTION)
			cathNorm[i] = WY_Norm[i];			
	}
	
	// Wirechamber variable numbers / dead wires
	// ordered from greatest to least in position
	assert(thisRun->RI.runNum >= 7000);
	
	// branches of interest
	std::string sensname = sideSubst("MWPC%c",mySide);
	sensname += (myDirection==X_DIRECTION)?"x":"y";
	
	
	std::vector<bool> liveWires = getLiveWires(thisRun->RI.runNum,mySide,myDirection);
	std::vector<unsigned int> padcNums = getPadcNumbers(thisRun->RI.runNum,mySide,myDirection);
	nWires = 0;
	for(unsigned int i=0; i<liveWires.size(); i++) {
		if(!liveWires[i])
			continue;
		if(mySide == EAST)
			loadData(std::string("Pdc")+itos(padcNums[nWires]),sensname+itos(i+1),"cathode");
		else
			loadData(std::string("Padc")+itos(padcNums[nWires]),sensname+itos(i+1),"cathode");
		nWires++;
	}
	
	if(thisRun->RI.runNum < 7000) {
		if(mySide == EAST)
			loadData("Pdc20","MWPCEAnode","anode");
		else
			loadData("Padc0","MWPCWAnode","anode");
	} else {
		if(mySide == EAST)
			loadData("Pdc30","MWPCEAnode","anode");
		else
			loadData("Pdc34","MWPCWAnode","anode");
	}	
}

// function for determining cathode pedestals
std::vector< std::pair<float,float> > wire_pedestal_finder(Subsystem* S, void* cdat) {
	Wirechamber* W = (Wirechamber*)S;
	float* cathdat = (float*)cdat;
	std::vector< std::pair<float,float> > v;
	for(unsigned int e=0; e<W->nEvents; e++)
		if(!W->TG->beta2of4(e,W->mySide) && !W->TG->isBeamnoise(e))
			v.push_back(std::pair<float,float>(W->TG->eventTime(e),cathdat[e]));
	return v;
}


Wirechamber::Wirechamber(RunManager* T, Trigger* tg, Side s, AxisDirection d):
Subsystem(T,"",s), TG(tg), myDirection(d) {
	
	wireSpacing = 10.16*sqrt(0.6);
	specialize();
	wirePositions = calcWirePositions(thisRun->RI.runNum,mySide,myDirection);
	lowerEdge = -8*wireSpacing;
	upperEdge = 8*wireSpacing;
	
	for(unsigned int c=0; c<nWires; c++)
		cathdat.push_back(getData("cathode",c));
	hits = new wireHit[nEvents];
	
	// pedestals
	for (UInt_t i = 0; i<nWires; i++ ) 
		if(!verifyPedestal(i))
			monitorPeak(wire_pedestal_finder, 60.0, 3000, (void*)cathdat[i], sensorNames[i]);
	
	if(!(thisRun->runMode & RUNMODE_POSTRACK))
		return;

	printf("Tracking %i events on %i wires...\n",nEvents,nWires);

	//------
	// reconstruct each event
	//------
	
	
	float* wirePeds = new float[nWires];
	float* wireValues = new float[nWires];	
	for(unsigned int e=0; e<nEvents; e++) {
		for(unsigned int c=0; c<nWires; c++) {
			//wireValues[c] = cathdat[c][e]*500.0/cathNorm[c];
			//wirePeds[c] = PC.getPedestal(sensorNames[c],TG->eventTime(e))*500.0/cathNorm[c];
			wireValues[c] = cathdat[c][e];
			wirePeds[c] = PC.getPedestal(sensorNames[c],TG->eventTime(e));
		}
		hits[e] = mpmGaussianPositioner(wirePositions, wireValues, wirePeds);
	}
	delete(wirePeds);
	delete(wireValues);
	
	// cathode max cut
	TH1F* hCathMax = registeredTH1F("MaxValue","Wirechamber Maximum Cathode Value",250,-200,4200);
	TH1F* hCathMaxPed = registeredTH1F("MaxValuePed","Wirechamber Maximum Cathode Value Pedestal",1000,-200,4200);
	hCathMaxPed->SetLineColor(2);
	for(unsigned int e=0; e<nEvents; e++) {
		hCathMax->Fill(hits[e].maxValue);
		if(!TG->beta2of4(e,mySide))
			hCathMaxPed->Fill(hits[e].maxValue);
	}
	hCathMax->Scale(1.0/(TG->runTime(mySide)*hCathMax->GetBinWidth(1)));
	TF1* g1 = new TF1("g1","gaus",0, 200); 
	hCathMaxPed->Fit(g1,"QR");
	cathMaxCut = g1->GetParameter(1)+4.0*fabs(g1->GetParameter(2));
	float pheight = g1->GetParameter(0);
	g1->SetRange(g1->GetParameter(1)-2.0*fabs(g1->GetParameter(2)),g1->GetParameter(1)+2.0*fabs(g1->GetParameter(2)));
	hCathMax->Fit(g1,"QR");
	hCathMaxPed->Scale(g1->GetParameter(0)/pheight);
	
	forParent.insert("CathMaxCut",cathMaxCut);
	printf("Cathode Max Cut @ %.2f\n",cathMaxCut);
	
	defaultCanvas->SetLogy(true);
	hCathMax->Draw();
	hCathMaxPed->Draw("Same");
	drawVLine(cathMaxCut,defaultCanvas);
	printCanvas("MaxValue");
	defaultCanvas->SetLogy(false);
	
	// set triggers by max cathode cut
	int ntriggers = 0;
	for(UInt_t e = 0; e<nEvents; e++) {
		triggers[e] = hits[e].maxValue > cathMaxCut;
		if(triggers[e])
			ntriggers++;
	}
	printf(">>> %i wirechamber events (%.2f %%)\n\n",ntriggers,ntriggers*100.0/nEvents);
	
	// hit sigmas histogram
	TH1F* hSigma = registeredTH1F("Hits_Sigma","Wirechamber Hits Sigma",1000,0.0,30.0);
	for(unsigned int e=0; e<nEvents; e++) {
		if(triggers[e])
			hSigma->Fill(hits[e].width);
	}
	defaultCanvas->SetLogy(true);
	hSigma->Draw();
	printCanvas("HitSigma");
	defaultCanvas->SetLogy(false);
	
	// Backup position calculation for "bad" hits
	for(UInt_t e=0; e<nEvents; e++) {
		Float_t x = bradPosition(e);
		if(!triggers[e]) {
			hits[e].center = x;
		}
	}
	
	genHistograms();
}

void Wirechamber::genHistograms() {
	
	defaultCanvas->SetLogy(true);
	
	//----------------------------
	// cathode sum
	//----------------------------
	
	if(0) {
		TH1F* h1 = registeredTH1F("CathodeSum","Wirechamber cathode sum",500,-500,4000);
		for(UInt_t e=0; e<nEvents; e++)
			h1->Fill(hits[e].cathodeSum);
		h1->Draw();
		drawVLine(csumPedestal,defaultCanvas);
		drawVLine(csumPedestal+4*csumPedwidth,defaultCanvas);
		printCanvas("Cathode_Sum");
	}
	
	//----------------------------
	// individual wire spectra
	//----------------------------
	
	if(0) {
		for(unsigned int i=0; i<nWires; i++) {
			TH1F* h2 = registeredTH1F(std::string("Wire_Spectrum_")+itos(i),"Wire ADC",520,-20,3900);
			for(UInt_t e=0; e<nEvents; e++)
				h2->Fill(cathdat[i][e]);
			
			h2->Draw();
			drawVLine(PC.getPedestal(sensorNames[i],0.0),defaultCanvas);
			printCanvas(std::string("Wire_")+itos(i));
		}
	}
	
	
	//----------------------------
	// wire rates
	//----------------------------
	
	std::vector<TH1F*> wireRates;
	TH1F totalRate("TotalRate","Wirechamber Rates",100,0,TG->totalTime());
	for(unsigned int i=0; i<nWires; i++) {
		wireRates.push_back(registeredTH1F(std::string("Wire_Rate_")+itos(i),"Wire Rates",40,0,TG->totalTime()));
		wireRates[i]->SetLineColor(i+1);
	}
	
	for(unsigned int e=0; e<nEvents; e++) {
		if(triggers[e]) {
			wireRates[hits[e].maxWire]->Fill(TG->eventTime(e));
			totalRate.Fill(TG->eventTime(e));
		}
	}
	
	totalRate.Scale(1.0/totalRate.GetBinWidth(1));
	for(unsigned int i=0; i<nWires; i++)
		wireRates[i]->Scale(1.0/wireRates[i]->GetBinWidth(1));
	
	totalRate.SetMinimum(1e-2);
	totalRate.SetMaximum(1e2);	
	if(totalRate.GetEntries()) {
		totalRate.Draw();
		for(unsigned int i=0; i<nWires; i++)
			wireRates[i]->Draw("Same");
		printCanvas("Rates");
	}
	
	//----------------------------
	// multiplicity
	//----------------------------
	
	TH1F* h3 = registeredTH1F("Multiplicity","Wirechamber Multiplicity",nWires+3,-1.5,nWires+1.5);
	for(UInt_t e=0; e<nEvents; e++)
		if(TG->beta2of4(e, mySide))
			h3->Fill(hits[e].multiplicity);
	if(h3->GetEntries()) {
		h3->Draw();
		printf("%i of %i events have multiplicity 1.\n",(int)h3->GetBinContent(3),(int)h3->GetEntries());
		printCanvas("Multiplicity");
	}
	
	//----------------------------
	// source positions & peak fits
	//----------------------------
	
	// positions
	TH1F* h6 = registeredTH1F("Positions","Wirechamber Positions",400,getLowerEdge(),getUpperEdge());
	for(UInt_t e=0; e<nEvents; e++)
		if(triggers[e] && TG->beta2of4(e, mySide))
			h6->Fill(hits[e].center);
	h6->Scale(1.0/(TG->runTime(mySide)*h6->GetBinWidth(1)));
	
	if(0) {
		std::vector<TF1*> sourcefits;
		std::vector<Source>::const_iterator it;
		
		// fit peaks
		for(it = thisRun->SI.sourcesBegin(mySide); it != thisRun->SI.sourcesEnd(mySide); it++) {
			float pc;
			if(myDirection == X_DIRECTION)
				pc = it->x;
			else
				pc = it->y;
			sourcefits.push_back(new TF1((std::string("source_")+it->t).c_str(),"gaus",pc-15,pc+15));
			h6->Fit(sourcefits.back(),"QR+");
			//thisRun->addObject(sourcefits.back());
		}
		
		h6->Draw();
		
		// label peaks
		unsigned int i = 0;
		for(it = thisRun->SI.sourcesBegin(mySide); it != thisRun->SI.sourcesEnd(mySide); it++) {
			float pc;
			if(myDirection == X_DIRECTION)
				pc = it->x;
			else
				pc = it->y;			
			TLatex lx;
			char tmp[1024];
			sprintf(tmp,"#color[2]{#sigma = %.2f}",sourcefits[i]->GetParameter(2));
			lx.DrawLatex(pc-10,100,tmp);
			sprintf(tmp,"#color[2]{#mu = %.2f}",sourcefits[i]->GetParameter(1));
			lx.DrawLatex(pc-10,500,tmp);
			i++;
		}
	}
	if(h6->GetEntries()) {
		h6->Draw();
		// wire locations
		for(unsigned int i=0; i<nWires; i++)
			drawVLine(wirePositions[i],defaultCanvas); 
		printCanvas("Positions");
	}
	
	//----------------------------
	// max wire number
	//----------------------------
	
	TH1F* h4 = registeredTH1F("MaxWire","Wirechamber Maximum Wire",nWires+2,-1.5,nWires+0.5);
	for(UInt_t e=0; e<nEvents; e++) {
		if(triggers[e] && TG->beta2of4(e, mySide))
			h4->Fill(hits[e].maxWire);
	}
	h4->Scale(1.0/h4->Integral());
	if(h4->GetEntries()) {
		h4->Draw();
		printCanvas("MaxWire");
	}
	
	//----------------------------
	// clipping
	//----------------------------
	
	TH1F* hClip = registeredTH1F("Clipping","Wirechamber Clipping",nWires+3,-1.5,nWires+1.5);
	for(UInt_t e=0; e<nEvents; e++)
		hClip->Fill(hits[e].nClipped);
	if(hClip->GetEntries()) {
		hClip->Scale(1.0/hClip->Integral());
		hClip->Draw();
		printCanvas("Clipping");
	}
	defaultCanvas->SetLogy(false);
}






Float_t Wirechamber::bradPosition(UInt_t e) {
	Float_t cathode_signal[nWires];
	
	// pedestal-subtracted values with minimum cuts
	for(unsigned int i=0; i<nWires; i++) {
		cathode_signal[i] = cathdat[i][e];
		if(cathdat[i][e] > 4079)
			cathode_signal[i] = PC.getPedestal(sensorNames[i],TG->eventTime(e));
		cathode_signal[i] -= PC.getPedestal(sensorNames[i],TG->eventTime(e));
		if(cathode_signal[i] < 70)
			cathode_signal[i] = 0;
	}
	
	// isolated wire removal
	if(cathode_signal[0] > 0 && cathode_signal[1] == 0)
		cathode_signal[0] = 0;
	for(unsigned int i=1; i<nWires-1; i++)
		if(cathode_signal[i] > 0 && cathode_signal[i-1]==0 && cathode_signal[i+1]==0)
			cathode_signal[i] = 0;
	if(cathode_signal[nWires-1] > 0 && cathode_signal[nWires-2] == 0)
		cathode_signal[nWires-1] = 0;
	
	// wire multiplicity and sums
	hits[e].multiplicity = 0;
	Float_t cathode_sum = 0;
	for(unsigned int j=0; j<nWires; j++) {
		if (cathode_signal[j] > 0) {
			hits[e].multiplicity++;
			cathode_sum += cathode_signal[j];
		}
	}							
	
	// Compute MWPC positions and project MWPC coordinates back to spectrometer coordinates.
	Float_t sum1 = 0;
	Float_t sum2 = 0;
	
	for(unsigned int j=0; j<nWires; j++) {
		sum1 += wirePositions[j] * cathode_signal[j];
		sum2 += cathode_signal[j];
	}
	Float_t x;
	if (hits[e].multiplicity > 0)
		x = sum1 / sum2 * sqrt(0.6);
	else
		x = -1000;
	
	return x;
}
