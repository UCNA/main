
#include "WirechamberStudy.hh"
#include "CalDBSQL.hh"
#include "MultiGaus.hh"
#include "PostOfficialAnalyzer.hh"
#include "GraphicsUtils.hh"
#include "GraphUtils.hh"
#include <TDirectory.h>
#include <TROOT.h>
#include <TStyle.h>
#include <time.h>

Stringmap anodeseg2sm(const AnodeSeg& a) {
	Stringmap m;
	m.insert("side",sideSubst("%c",a.s));
	m.insert("segment",a.segment);
	m.insert("fiterr",a.fiterr);
	m.insert("energybin",a.energybin);
	m.insert("energy",a.energy);
	m.insert("npts",a.npts);
	m.insert("constant",a.constant.x);
	m.insert("mpv",a.mpv.x);
	m.insert("sigma",a.sigma.x);
	m.insert("d_constant",a.constant.err);
	m.insert("d_mpv",a.mpv.err);
	m.insert("d_sigma",a.sigma.err);
	return m;
}

AnodeSeg sm2anodeseg(const Stringmap& m) {
	AnodeSeg a;
	a.s=(m.getDefault("side","E")=="E")?EAST:WEST;
	a.segment=m.getDefault("segment",0);
	a.fiterr=m.getDefault("fiterr",0);
	a.energybin=m.getDefault("energybin",0);
	a.npts=m.getDefault("npts",0);
	a.constant=float_err(m.getDefault("constant",0),m.getDefault("d_constant",0));
	a.mpv=float_err(m.getDefault("mpv",0),m.getDefault("d_mpv",0));
	a.sigma=float_err(m.getDefault("sigma",0),m.getDefault("d_sigma",0));
	return a;
}

Stringmap cathseg2sm(const CathodeSeg& c) {
	Stringmap m;
	m.insert("side",sideSubst("%c",c.s));
	m.insert("direction",c.d);
	m.insert("i",c.i);
	m.insert("height",c.height.x);
	m.insert("d_height",c.height.err);
	m.insert("width",c.width.x);
	m.insert("d_width",c.width.err);
	m.insert("position",c.pos);
	return m;
}

CathodeSeg sm2cathseg(const Stringmap& m) {
	CathodeSeg c;
	c.s=(m.getDefault("side","E")=="E")?EAST:WEST;
	c.d=AxisDirection(int(m.getDefault("direction",0)));
	c.i=m.getDefault("i",0);
	c.pos=m.getDefault("position",0);
	c.height=float_err(m.getDefault("height",0),m.getDefault("d_height",0));
	c.width=float_err(m.getDefault("width",0),m.getDefault("d_width",0));
	return c;
}

WirechamberAnalyzer::WirechamberAnalyzer(OutputManager* pnt, const std::string& nm, float r, unsigned int nr, const std::string& infl):
RunAccumulator(pnt,nm,infl), eMax(1000), nEnergyBins(40), sects(nr,r) {
	
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
	
	// set up anode histograms, data
	masterEnergySpectrum = registerFGBGPair("hMasterEnergy","Event Energy",nEnergyBins,0,eMax,AFP_OTHER,BOTH);
	for(Side s = EAST; s <= WEST; ++s) {
		for(unsigned int m=0; m<sects.nSectors(); m++) {
			anodeDat[s].push_back(std::vector<AnodeSeg>());
			anodeSpectra[s].push_back(std::vector<fgbgPair>());
			for(unsigned int e=0; e<nEnergyBins; e++) {
				anodeSpectra[s][m].push_back(registerFGBGPair(std::string("hAn_Sect_")+itos(m)+"_Bin_"+itos(e),"Anode ADC Spectrum",
															  100,0,2000,AFP_OTHER,s));
				anodeDat[s][m].push_back(AnodeSeg());
				anodeDat[s][m].back().s=s;
				anodeDat[s][m].back().segment=m;
				anodeDat[s][m].back().energybin=e;
				anodeDat[s][m].back().energy=masterEnergySpectrum.h[GV_OPEN]->GetBinCenter(e+1);
			}
		}
	}
	// set up cathode histograms, data
	for(Side s = EAST; s <= WEST; ++s) {
		for(AxisDirection d = X_DIRECTION; d <= Y_DIRECTION; ++d) {
			assert(false);
			//cathPos[s][d] = calcWirePositions(16000,s,AxisDirection(d));
			for(unsigned int i=0; i<cathPos[s][d].size(); i++) {
				TH2F hCathodeTemplate((sideSubst("hCath_%c",s)+(d==X_DIRECTION?"x_":"y_")+itos(i)).c_str(),"Normalized Cathode",100,-20,20,100,-0.5,5);
				cathHists[s][d].push_back(registerFGBGPair(hCathodeTemplate,AFP_OTHER,s));
				for(unsigned int n=0; n<2; n++)
					cathSlices[s][d][n].push_back(NULL);
				cathDat[s][d].push_back(CathodeSeg());
				cathDat[s][d].back().s=s;
				cathDat[s][d].back().d=d;
				cathDat[s][d].back().i=i;
				cathDat[s][d].back().pos=cathPos[s][d][i];
			}
		}
	}
	
	// load old data
	if(fIn) {
		QFile qOld(inflname+".txt");
		std::vector<Stringmap> sds = qOld.retrieve("anodeSeg");
		for(std::vector<Stringmap>::iterator it = sds.begin(); it != sds.end(); it++) {
			AnodeSeg a = sm2anodeseg(*it);
			assert(a.s<=WEST && a.segment<sects.nSectors() && a.energybin<nEnergyBins);
			anodeDat[a.s][a.segment][a.energybin] = a;
		}
		std::vector<Stringmap> cs = qOld.retrieve("cathodeSeg");
		for(std::vector<Stringmap>::iterator it = sds.begin(); it != sds.end(); it++) {
			CathodeSeg c = sm2cathseg(*it);
			assert(c.s<=WEST && c.d<=Y_DIRECTION && c.i<cathDat[c.s][c.d].size());
			cathDat[c.s][c.d][c.i] = c;
		}		
	}
}

void WirechamberAnalyzer::fillCoreHists(ProcessedDataScanner& PDS, double weight) {
	const Side s = PDS.fSide;
	if(!(PDS.fType == TYPE_0_EVENT && PDS.fPID == PID_BETA && (s==EAST||s==WEST))) return;
	
	if(!PDS.isSimulated()) {
		for(AxisDirection d = X_DIRECTION; d <= Y_DIRECTION; ++d)
			for(unsigned int i=0; i<cathPos[s][d].size(); i++)
				if(!(PDS.wires[s][d].errflags & (WIRES_NONE | WIRES_SINGLET | WIRES_CLIPPED | WIRES_DOUBLET)))
					cathHists[s][d][i].h[currentGV]->Fill(PDS.wires[s][d].center-cathPos[s][d][i],PDS.cathodes[s][d][i]/PDS.mwpcs[s].anode);
	}
	
	unsigned int m = sects.sector(PDS.wires[s][X_DIRECTION].center,PDS.wires[s][Y_DIRECTION].center);
	if(m>=sects.nSectors()) return;
	masterEnergySpectrum.h[currentGV]->Fill(PDS.getEtrue());
	unsigned int ebin = masterEnergySpectrum.h[currentGV]->FindBin(PDS.getEtrue())-1;
	if(ebin>=nEnergyBins) return;
	anodeSpectra[s][m][ebin].h[currentGV]->Fill(PDS.mwpcs[s].anode);
}


std::vector<TH1D*> replaceFitSlicesY(TH2* h, TF1* f1=NULL) {
		
	Int_t nbins  = h->GetXaxis()->GetNbins();
	std::vector<TH1D*> hlist;
	
	//default is to fit with a gaussian
	if (!f1) {
		double ymin=h->GetYaxis()->GetXmin();
		double ymax=h->GetYaxis()->GetXmax();
		f1 = (TF1*)gROOT->GetFunction("gaus");
		if(!f1) f1 = new TF1("gaus","gaus",ymin,ymax);
		else f1->SetRange(ymin,ymax);
	}
	Int_t npar = f1->GetNpar();
	if (npar <= 0) return hlist;
	Double_t *parsave = new Double_t[npar];
	f1->GetParameters(parsave);
	
	//Create one histogram for each function parameter
	const TArrayD* bins = h->GetXaxis()->GetXbins();
	for (Int_t ipar=0; ipar<=npar; ipar++) {
		std::string name = std::string(h->GetName())+"_"+itos(ipar);
		std::string title = std::string("Fit parameter ")+itos(ipar);
		delete gDirectory->FindObject(name.c_str());
		if (bins->fN == 0) {
			hlist.push_back(new TH1D(name.c_str(), title.c_str(), nbins, h->GetXaxis()->GetXmin(),  h->GetXaxis()->GetXmax()));
		} else {
			hlist.push_back(new TH1D(name.c_str(), title.c_str(), nbins,bins->fArray));
		}
		hlist.back()->GetXaxis()->SetTitle( h->GetXaxis()->GetTitle());
	}
	
	//Loop on all bins in X, generate a projection along Y
	for (Int_t bin=1; bin<=nbins; bin++) {
		TH1D *hpy = h->ProjectionY("_temp",bin,bin,"e");
		if (hpy == 0) continue;
		Int_t nentries = Int_t(hpy->GetEntries());
		if (nentries == 0) {delete hpy; continue;}
		f1->SetParameters(parsave);
		hpy->Fit(f1,"QNR");
		Int_t npfits = f1->GetNumberFitPoints();
		if (npfits > npar) {
			for (Int_t ipar=0; ipar<npar; ipar++) {
				hlist[ipar]->Fill(h->GetXaxis()->GetBinCenter(bin),f1->GetParameter(ipar));
				hlist[ipar]->SetBinError(bin,f1->GetParError(ipar));
			}
			hlist[npar]->Fill(h->GetXaxis()->GetBinCenter(bin),f1->GetChisquare()/(npfits-npar));
		}
		delete hpy;
	}
	
	delete [] parsave;
	return hlist;
}

void WirechamberAnalyzer::calculateResults() {
	
	// fit anode spectrum in each segment, energy
	printf("\n------------------ Fitting anode data...\n");
	TF1 fLandau("landauFit","landau",0,1000);
	fLandau.SetLineColor(2);
	for(Side s = EAST; s <= WEST; ++s) {
		for(unsigned int m = 0; m < sects.nSectors(); m++) {
			for(unsigned int e = 0; e < nEnergyBins; e++) {
				AnodeSeg& a = anodeDat[s][m][e];
				a.npts = anodeSpectra[s][m][e].h[GV_OPEN]->Integral();
				if(a.npts >= 100) {
					a.fiterr = anodeSpectra[s][m][e].h[GV_OPEN]->Fit(&fLandau,"QR");
					a.constant = float_err(fLandau.GetParameter(0),fLandau.GetParError(0));
					a.mpv = float_err(fLandau.GetParameter(1),fLandau.GetParError(1));
					a.sigma = float_err(fLandau.GetParameter(2),fLandau.GetParError(2));
				}
				qOut.insert("anodeSeg",anodeseg2sm(a));
			}
		}
	}
	
	// fit hits distribution around each cathode
	printf("\n-------------------- Fitting cathode data...\n");
	TF1 fGaus("fGaus","gaus",-5,5);
	fGaus.SetLineColor(2);
	for(Side s = EAST; s <= WEST; ++s) {
		for(AxisDirection d = X_DIRECTION; d <= Y_DIRECTION; ++d) {
			for(unsigned int i=0; i<cathPos[s][d].size(); i++) {
				if(!cathHists[s][d][i].h[GV_OPEN]->Integral()) continue;
				std::vector<TH1D*> slicefits = replaceFitSlicesY((TH2F*)(cathHists[s][d][i].h[GV_OPEN]));
				for(std::vector<TH1D*>::iterator it = slicefits.begin(); it != slicefits.end(); it++) {
					(*it)->SetDirectory(0);
					addObject(*it);
				}
				cathSlices[s][d][0][i] = slicefits[1];
				cathSlices[s][d][1][i] = slicefits[2];
				cathSlices[s][d][0][i]->Fit(&fGaus,"QR");
				CathodeSeg& c = cathDat[s][d][i];
				c.height = float_err(fGaus.GetParameter(0),fGaus.GetParError(0));
				c.width = float_err(fGaus.GetParameter(2),fGaus.GetParError(2));
				qOut.insert("cathodeSeg",cathseg2sm(c));
			}
		}
	}
}

void WirechamberAnalyzer::makePlots() {
	defaultCanvas->cd();
	
	// binned energy spectrum
	masterEnergySpectrum.h[GV_OPEN]->Draw();
	printCanvas("MasterEnergy");
	
	// anode spectra
	std::vector<TH1*> hToPlot;
	for(Side s = EAST; s <= WEST; ++s) {
		for(unsigned int m = 0; m < sects.nSectors(); m++) {
			hToPlot.clear();
			for(unsigned int e = 0; e < nEnergyBins; e++)
				hToPlot.push_back(anodeSpectra[s][m][e].h[GV_OPEN]);
			drawSimulHistos(hToPlot);
			printCanvas(sideSubst("Anodes/%c_",s)+itos(m));
		}
	}
	
	// cathode plots
	for(Side s = EAST; s <= WEST; ++s) {
		for(AxisDirection d = X_DIRECTION; d <= Y_DIRECTION; ++d) {
			for(unsigned int i=0; i<cathPos[s][d].size(); i++) {
				cathHists[s][d][i].h[GV_OPEN]->Draw("Col");
				if(cathSlices[s][d][0][i])
					cathSlices[s][d][0][i]->Draw("Same");
				if(cathSlices[s][d][1][i])
					cathSlices[s][d][1][i]->Draw("Same");
				printCanvas(sideSubst("Cathodes/Cathode_%c",s)+(d==X_DIRECTION?"x_":"y_")+itos(i));
			}
		}
	}
}

void WirechamberAnalyzer::compareMCtoData(RunAccumulator& OAdata) {
}

void runWirechamberAnalyzer(RunNum r0, RunNum r1, unsigned int nrings) {
	
	double fidRadius = 50;
	
	// scan data from each run
	std::vector<std::string> snames;
	OutputManager OM1("NameUnused",getEnvSafe("UCNA_ANA_PLOTS")+"/WirechamberMaps/SingleRuns/");
	for(RunNum r = r0; r <= r1; r++) {
		std::string singleName = std::string("WC_")+itos(r)+"_"+itos(nrings)+"_"+dtos(fidRadius);
		std::string prevFile = OM1.basePath+"/"+singleName+"/"+singleName;
		snames.push_back(singleName);
		if(r0==r1 || !fileExists(prevFile+".root")) {
			WirechamberAnalyzer WA1(&OM1, singleName, fidRadius, nrings);
			PostOfficialAnalyzer POA(true);
			POA.addRun(r);
			WA1.loadProcessedData(AFP_OTHER, GV_OPEN, POA);
			WA1.write();
			WA1.setWriteRoot(true);
		}
	}
	
	if(r0==r1) return;
	
	// reload data
	OutputManager OM("NameUnused",getEnvSafe("UCNA_ANA_PLOTS")+"/WirechamberMaps/");
	WirechamberAnalyzer WA(&OM, std::string("WC_")+itos(r0)+"-"+itos(r1), fidRadius, nrings);	
	for(std::vector<std::string>::iterator it = snames.begin(); it != snames.end(); it++) {
		std::string prevFile = OM1.basePath+"/"+*it+"/"+*it;
		WirechamberAnalyzer WA1(&OM1, *it, 0, 0, prevFile);
		WA.addSegment(WA1);
	}
	
	// finish and output
	WA.calculateResults();
	WA.makePlots();
	WA.write();
	WA.setWriteRoot(true);	
}
