#include "PositionStudies.hh"
#include "Enums.hh"
#include "OutputManager.hh"
#include "GraphicsUtils.hh"
#include "SectorCutter.hh"
#include "WirechamberReconstruction.hh"
#include <math.h>
#include <vector>
#include <TH1F.h>
#include <TF1.h>
#include <TDirectory.h>



void AnodeCalibration(ProcessedDataScanner& PDS, OutputManager& OM, unsigned int nRings) {
	
	const float fidRadius = 50.0;
	SectorCutter sects(nRings,fidRadius);
	const int eMax = 1000;
	const int nEnergyBins = 40;
	/// master energy spectrum histogram; really used for binning energy
	TH1F* masterEnergySpectrum = OM.registeredTH1F("MasterEnergy","Type 0 Event Energy",nEnergyBins,0,eMax);
	std::vector<TH1F*> anodeSpectra[2][nEnergyBins];	//< anode spectra for each [side][energy][position]
	float bmax = PDS.isSimulated()?20:4000;
	for(Side s = EAST; s <= WEST; s = nextSide(s))
		for(int e = 0; e < nEnergyBins; e++)
			for(unsigned int p = 0; p < sects.nSectors(); p++)
				anodeSpectra[s][e].push_back(OM.registeredTH1F(sideSubst("Anode_%c_",s)+itos(p)+"_"+itos(e),"Anode ADC",200,0,bmax));
	
	// collect data
	printf("Scanning data...\n");
	PDS.startScan();
	while(PDS.nextPoint()) {
		Side s = PDS.fSide;
		if(PDS.fType != TYPE_0_EVENT || PDS.fPID != PID_BETA || !(s == EAST || s == WEST) || PDS.radius(s) > fidRadius) continue;
		unsigned int p = sects.sector(PDS.wires[s][X_DIRECTION].center,PDS.wires[s][Y_DIRECTION].center);
		int e = masterEnergySpectrum->FindBin(PDS.getEtrue())-1;
		if(p<0 || p>=sects.nSectors() || e<0 || e>=nEnergyBins) continue;
		masterEnergySpectrum->Fill(PDS.getEtrue());
		anodeSpectra[s][e][p]->Fill(PDS.mwpcs[s].anode);
	}
	
	masterEnergySpectrum->Draw();
	OM.printCanvas("MasterEnergy");

	// fit each cell, save results
	printf("Fitting data...\n");
	TF1 fLandau("landauFit","landau",0,bmax/2);
	fLandau.SetLineColor(2);
	for(Side s = EAST; s <= WEST; s = nextSide(s)) {
		for(unsigned int p = 0; p < sects.nSectors(); p++) {
			std::vector<TH1*> hToPlot;
			for(int e = 0; e < nEnergyBins; e++) {
				int npts = anodeSpectra[s][e][p]->Integral();
				if(npts < 100) continue;
				int err = anodeSpectra[s][e][p]->Fit(&fLandau,"QR");
				if(!err)
					hToPlot.push_back(anodeSpectra[s][e][p]);
				float x,y;
				sects.sectorCenter(p,x,y);
				Stringmap m;
				m.insert("side",sideSubst("%c",s));
				m.insert("position",p);
				m.insert("fiterr",err);
				m.insert("x",x);
				m.insert("y",y);
				m.insert("energyBin",e);
				m.insert("energy",masterEnergySpectrum->GetBinCenter(e+1));
				m.insert("npts",npts);
				m.insert("constant",fLandau.GetParameter(0));
				m.insert("mpv",fLandau.GetParameter(1));
				m.insert("sigma",fLandau.GetParameter(2));
				m.insert("d_constant",fLandau.GetParError(0));
				m.insert("d_mpv",fLandau.GetParError(1));
				m.insert("d_sigma",fLandau.GetParError(2));
				OM.qOut.insert("anodeCal",m);
			}
			if(!hToPlot.size()) continue;
			OM.defaultCanvas->SetLogy(true);
			drawSimulHistos(hToPlot);
			OM.printCanvas(sideSubst("Anodes_%c_",s)+itos(p));
		}
	}
	
	Stringmap m;
	m.insert("nRings",sects.n);
	m.insert("radius",sects.r);
	m.insert("nSectors",sects.nSectors());
	OM.qOut.insert("SectorCutter",m);
	
	OM.write();
	OM.setWriteRoot(true);
}

void CathodeCalibration(PostOfficialAnalyzer& POA, OutputManager& OM) {
	
	/// nominal wire positions, corresponding histograms
	std::vector<float> cathPos[2][2];
	std::vector<TH2F*> cathHists[2][2];
	for(Side s = EAST; s <= WEST; s = nextSide(s)) {
		for(int d = X_DIRECTION; d <= Y_DIRECTION; d++) {
			cathPos[s][d] = calcWirePositions(16000,s,AxisDirection(d));
			for(unsigned int i=0; i<cathPos[s][d].size(); i++)
				cathHists[s][d].push_back(OM.registeredTH2F(sideSubst("hCath_%c",s)+(d==X_DIRECTION?"x_":"y_")+itos(i),
															"Normalized Cathode",100,-20,20,100,-0.5,5));
		}
	}
	
	/// gather data
	POA.startScan();
	while(POA.nextPoint()) {
		Side s = POA.fSide;
		if(POA.fType != TYPE_0_EVENT || POA.fPID != PID_BETA || !(s == EAST || s == WEST)) continue;
		for(int d = X_DIRECTION; d <= Y_DIRECTION; d++)
			for(unsigned int i=0; i<cathPos[s][d].size(); i++)
				if(!(POA.wires[s][d].errflags & (WIRES_NONE | WIRES_SINGLET | WIRES_CLIPPED | WIRES_DOUBLET)))
					cathHists[s][d][i]->Fill(POA.wires[s][d].center-cathPos[s][d][i],POA.cathodes[s][d][i]/POA.mwpcs[s].anode);
	}

	/// process and plot
	for(Side s = EAST; s <= WEST; s = nextSide(s)) {
		for(int d = X_DIRECTION; d <= Y_DIRECTION; d++) {
			for(unsigned int i=0; i<cathPos[s][d].size(); i++) {
				cathHists[s][d][i]->FitSlicesY(); // default to gaus fit
				TH1F* centers = NULL;
				for(unsigned int j=0; j<=2; j++) {
					TH1F* hFitz = (TH1F*)gDirectory->Get((std::string(cathHists[s][d][i]->GetName())+"_"+itos(j)).c_str());
					assert(hFitz);
					hFitz->SetDirectory(0);
					OM.addObject(hFitz);
					if(j==1)
						centers = hFitz;
				}
				TF1 fGaus("fGaus","gaus",-5,5);
				fGaus.SetLineColor(2);
				centers->Fit(&fGaus,"QR");
				Stringmap m;
				m.insert("side",sideSubst("%c",s));
				m.insert("direction",d);
				m.insert("number",i);
				m.insert("height",fGaus.GetParameter(0));
				m.insert("d_height",fGaus.GetParError(0));
				m.insert("width",fGaus.GetParameter(2));
				m.insert("d_width",fGaus.GetParError(2));
				m.insert("position",cathPos[s][d][i]);
				OM.qOut.insert("cathCal",m);
			
				cathHists[s][d][i]->Draw("Col");
				centers->Draw("Same");
				OM.printCanvas(sideSubst("Cathode_%c",s)+(d==X_DIRECTION?"x_":"y_")+itos(i));
			}
		}
	}
	
	OM.write();
	OM.setWriteRoot(true);
}

