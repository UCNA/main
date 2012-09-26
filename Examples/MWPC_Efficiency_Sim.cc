#include "BetaDecayAnalyzer.hh"
#include "GraphicsUtils.hh"
#include <TStyle.h>
#include <vector>

const unsigned int nThreshBins = 100;
const double maxThresh = 1.0;

void runSimulation() {
	
	unsigned int nToSim = 1e7;
	
	printf("Simulating MWPC threshold variation...\n");
	
	OutputManager OM("ThisNameIsNotUsedAnywhere",getEnvSafe("UCNA_ANA_PLOTS")+"/MWPC_Effic_Sim/");
	std::vector<BetaDecayAnalyzer*> AAs;
	for(unsigned int i=0; i<nThreshBins; i++) {
		AAs.push_back(new BetaDecayAnalyzer(&OM,"Threshold_"+itos(i)));
		AAs.back()->currentGV = GV_OPEN;
	}
	
	printf("Setting input file...\n");
	G4toPMT simData;
	simData.addFile("/home/mmendenhall/geant4/output/LivPhys_495_neutronBetaUnpol_geomC/analyzed_*.root");
	PMTCalibrator PCal(16000);
	simData.setCalibrator(PCal);
	
	printf("Scanning simulation data...\n");
	for(unsigned int p=0; p<nToSim; p++) {
		simData.nextPoint();
		for(AFPState afp=AFP_OFF; afp<=AFP_ON; ++afp) {
			simData.setAFP(afp);
			simData.calcReweight();
			for(unsigned int i=0; i<nThreshBins; i++) {
				simData.mwpcThresh[EAST]=simData.mwpcThresh[WEST]=maxThresh*float(i)/(nThreshBins-1);
				simData.classifyEvent();
				AAs[i]->currentAFP = afp;
				AAs[i]->setFillPoints(afp,GV_OPEN);
				AAs[i]->loadSimPoint(simData);
			}
		}
		if(!(p%(nToSim/20))) {
			printf("*");
			fflush(stdout);
		}
	}
	printf(" Done.\n");
	
	for(unsigned int i=0; i<nThreshBins; i++) {
		
		AAs[i]->runTimes += simData.runTimes;
		for(AFPState afp=AFP_OFF; afp<=AFP_ON; ++afp)
			AAs[i]->totalTime[afp][GV_OPEN] += simData.totalTime;
		
		AAs[i]->calculateResults();
		AAs[i]->makePlots();
		AAs[i]->write();
		AAs[i]->setWriteRoot(true);
		delete(AAs[i]);
	}
}

void SnMWPCEffic() {
	unsigned int nToSim = 1e5;
	
	printf("Simulating MWPC threshold variation...\n");
	OutputManager OM("ThisNameIsNotUsedAnywhere",getEnvSafe("UCNA_ANA_PLOTS")+"/MWPC_Effic_Sim/");
	
	printf("Setting input file...\n");
	G4toPMT simData;
	simData.addFile("/home/mmendenhall/geant4/output/20120823_Sn113/analyzed_*.root");
	PMTCalibrator PCal(16000);
	simData.setCalibrator(PCal);
	simData.setAFP(AFP_OTHER);
	
	TH1F* hSnCounts[2];
	TH1F* hSnMWPC[2];
	for(Side s = EAST; s <= WEST; ++s) {
		hSnCounts[s] = new TH1F(sideSubst("hSnEffic_%c",s).c_str(),"Sn Source Type 0 Detection Efficiency",nThreshBins,0,maxThresh);
		hSnCounts[s]->GetXaxis()->SetTitle("East MWPC Threshold");
		hSnCounts[s]->GetYaxis()->SetTitle("Type 0 Efficiency");
		hSnMWPC[s] = new TH1F(sideSubst("hSnMWPC_%c",s).c_str(),"Sn Source MWPC Spectrum",5*nThreshBins,0,5*maxThresh);
		hSnMWPC[s]->GetYaxis()->SetTitle("MWPC Energy");
	}
	
	printf("Scanning simulation data...\n");
	std::vector<double> thresholds;
	for(unsigned int i=0; i<nThreshBins; i++)
		thresholds.push_back(maxThresh*float(i)/(nThreshBins-1));
	for(unsigned int p=0; p<nToSim; p++) {
		simData.mwpcThresh[EAST]=0;
		simData.nextPoint();
		for(Side s = EAST; s <= WEST; ++s)
			if(simData.fType==TYPE_0_EVENT && simData.fPID == PID_BETA && simData.fSide==s)
				hSnMWPC[s]->Fill(simData.mwpcEnergy[s]);
		for(unsigned int i=0; i<nThreshBins; i++) {
			simData.mwpcThresh[EAST]=hSnCounts[EAST]->GetBinCenter(i+1);
			simData.classifyEvent();
			for(Side s = EAST; s <= WEST; ++s)
				if(simData.fType==TYPE_0_EVENT && simData.fPID == PID_BETA && simData.fSide==s)
					hSnCounts[s]->Fill(simData.mwpcThresh[EAST]);
		}
		if(!(p%(nToSim/20))) {
			printf("+");
			fflush(stdout);
		}
	}
	
	printf(" Done.\n");
	for(Side s = EAST; s <= WEST; ++s) {
		hSnCounts[s]->Scale(1.0/hSnCounts[s]->GetBinContent(1));
		hSnCounts[s]->GetYaxis()->SetRangeUser(.995,1.001);
	}
	
	drawHistoPair(hSnCounts[EAST],hSnCounts[WEST]);
	OM.printCanvas("SnEffic");
	drawHistoPair(hSnMWPC[EAST],hSnMWPC[WEST]);
	OM.printCanvas("SnMWPC");
}

void processSimulation() {
	OutputManager OM("MWPC_Effic_Compare",getEnvSafe("UCNA_ANA_PLOTS")+"/MWPC_Effic_Sim/");
	std::vector<BetaDecayAnalyzer*> AAs;
	for(unsigned int i=0; i<nThreshBins; i++) {
		std::string tname = "Threshold_"+itos(i);
		AAs.push_back(new BetaDecayAnalyzer(&OM,tname,OM.basePath+'/'+tname+'/'+tname));
		AAs.back()->currentGV = GV_OPEN;
	}
	
	std::vector<TH1*> hDeltaAsym;
	std::vector<TH1*> hDeltaSpectrum;
	for(unsigned int i=0; i<nThreshBins; i++) {
		AAs[i]->calculateResults();
		hDeltaAsym.push_back((TH1*)AAs[i]->myAsym->hAsym->Clone(("hAsym_"+itos(i)).c_str()));
		hDeltaAsym.back()->Add(AAs[0]->myAsym->hAsym,-1.0);
		hDeltaAsym.back()->Divide(AAs[0]->myAsym->hAsym);
		hDeltaAsym.back()->SetLineColor( (i%6)+1 );
		hDeltaAsym.back()->Scale(100.0);
		hDeltaAsym.back()->SetMinimum(-5.);
		hDeltaAsym.back()->SetMaximum(0);
		hDeltaAsym.back()->SetTitle("Asymmetry % Change");
		hDeltaAsym.back()->GetXaxis()->SetTitle("Energy [keV]");
		
		hDeltaSpectrum.push_back((TH1*)AAs[i]->myAsym->hSuperSum->Clone(("hSS_"+itos(i)).c_str()));
		hDeltaSpectrum.back()->Add(AAs[0]->myAsym->hSuperSum,-1.0);
		hDeltaSpectrum.back()->Divide(AAs[0]->myAsym->hSuperSum);
		hDeltaSpectrum.back()->Scale(-100.0);
		hDeltaSpectrum.back()->SetTitle("Type 0 % Inefficiency");
		hDeltaSpectrum.back()->SetLineColor( (i%6)+1 );
		hDeltaSpectrum.back()->GetXaxis()->SetTitle("Energy [keV]");
		hDeltaSpectrum.back()->SetMinimum(0.);
		hDeltaSpectrum.back()->SetMaximum(15.);
	}
	OM.defaultCanvas->cd();
	drawSimulHistos(hDeltaAsym,"HIST");
	OM.printCanvas("DeltaAsym");
	
	OM.defaultCanvas->cd();
	drawSimulHistos(hDeltaSpectrum,"HIST");
	OM.printCanvas("DeltaSpectrum");
	
	for(unsigned int i=0; i<nThreshBins; i++)
		delete(AAs[i]);
}


int main(int argc, char *argv[]) {
	
	gStyle->SetOptStat("");
	gStyle->SetPalette(1);
	gStyle->SetNumberContours(255);
	
	SnMWPCEffic();
	
	//if(argc<=1)
	//	runSimulation();
	//else
	//	processSimulation();
	
	return 0;
}
