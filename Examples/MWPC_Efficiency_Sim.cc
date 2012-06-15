#include "AsymmetryAnalyzer.hh"
#include "GraphicsUtils.hh"
#include <TStyle.h>
#include <vector>

const unsigned int nThreshBins = 20;
const double maxThresh = 0.75;

void runSimulation() {
	
	unsigned int nToSim = 1e7;
	
	printf("Simulating MWPC threshold variation...\n");
	
	OutputManager OM("ThisNameIsNotUsedAnywhere",getEnvSafe("UCNA_ANA_PLOTS")+"/MWPC_Effic_Sim/");
	std::vector<AsymmetryAnalyzer*> AAs;
	for(unsigned int i=0; i<nThreshBins; i++) {
		AAs.push_back(new AsymmetryAnalyzer(&OM,std::string("Threshold_")+itos(i)));
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

void processSimulation() {
	OutputManager OM("MWPC_Effic_Compare",getEnvSafe("UCNA_ANA_PLOTS")+"/MWPC_Effic_Sim/");
	std::vector<AsymmetryAnalyzer*> AAs;
	for(unsigned int i=0; i<nThreshBins; i++) {
		std::string tname = std::string("Threshold_")+itos(i);
		AAs.push_back(new AsymmetryAnalyzer(&OM,tname,OM.basePath+'/'+tname+'/'+tname));
		AAs.back()->currentGV = GV_OPEN;
	}
	
	std::vector<TH1*> hDeltaAsym;
	std::vector<TH1*> hDeltaSpectrum;
	for(unsigned int i=0; i<nThreshBins; i++) {
		AAs[i]->calculateResults();
		hDeltaAsym.push_back((TH1*)AAs[i]->hAsym->Clone((std::string("hAsym_")+itos(i)).c_str()));
		hDeltaAsym.back()->Add(AAs[0]->hAsym,-1.0);
		hDeltaAsym.back()->Divide(AAs[0]->hAsym);
		hDeltaAsym.back()->SetLineColor( (i%6)+1 );
		hDeltaAsym.back()->Scale(100.0);
		hDeltaAsym.back()->SetMinimum(-5.);
		hDeltaAsym.back()->SetMaximum(0);
		hDeltaAsym.back()->SetTitle("Asymmetry % Change");
		hDeltaAsym.back()->GetXaxis()->SetTitle("Energy [keV]");
		
		hDeltaSpectrum.push_back((TH1*)AAs[i]->hSuperSum->Clone((std::string("hSS_")+itos(i)).c_str()));
		hDeltaSpectrum.back()->Add(AAs[0]->hSuperSum,-1.0);
		hDeltaSpectrum.back()->Divide(AAs[0]->hSuperSum);
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
	
	if(argc<=1)
		runSimulation();
	else
		processSimulation();
	
	return 0;
}
