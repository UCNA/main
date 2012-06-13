#include "G4toPMT.hh"
#include "PenelopeToPMT.hh"
#include "OutputManager.hh"
#include "PathUtils.hh"
#include "GraphicsUtils.hh"
#include <TROOT.h>
#include <TStyle.h>

void mc_compare_plots(OutputManager& OM, Sim2PMT& SP1, Sim2PMT& SP2) {
	std::vector<Sim2PMT*> sps;
	sps.push_back(&SP1);
	sps.push_back(&SP2);
	
	std::vector<TH1F*> hEvis[TYPE_II_EVENT+1];
	std::vector<TH1F*> hCTScint;
	std::vector<double> t0norm;
	
	for(unsigned int i=0; i<sps.size(); i++) {
		Sim2PMT* SP = sps[i];
		
		// define histograms
		for(EventType t=TYPE_0_EVENT; t <= TYPE_II_EVENT; ++t) {
			hEvis[t].push_back(OM.registeredTH1F("hEvis_"+itos(t)+"_"+itos(i),
												 "Type "+itos(t)+" Energy Spectrum",
												 100,0,1000));
		}
		hCTScint.push_back(OM.registeredTH1F("hCTScint_"+itos(i),"Cos Theta entering scintillator",100,0.,1.0));
		
		// book histograms
		SP->startScan();
		while(SP->nextPoint()) {
			if(SP->fPID != PID_BETA) continue;
			if(SP->fType <= TYPE_II_EVENT)
				hEvis[SP->fType].back()->Fill(SP->getEnergy(),SP->physicsWeight);
			if(SP->fSide <= WEST)
				hCTScint.back()->Fill(SP->cosThetaInScint[SP->fSide],SP->physicsWeight);
		}
		
		t0norm.push_back(hEvis[TYPE_0_EVENT].back()->Integral());
	}
	
	/////////////
	// process data and make plots
	/////////////
	
	std::vector<TH1*> hToPlot;
	OM.defaultCanvas->cd();
	
	for(EventType t=TYPE_0_EVENT; t <= TYPE_II_EVENT; ++t) {
		hToPlot.clear();
		for(unsigned int i=0; i<sps.size(); i++) {
			hEvis[t][i]->SetLineColor(2+2*i);
			hEvis[t][i]->Scale(1.0/t0norm[0]);
			hToPlot.push_back(hEvis[t][i]);
			Stringmap m;
			m.insert("type",t);
			m.insert("sim",i);
			m.insert("mean",hEvis[t][i]->GetMean());
			m.insert("rms",hEvis[t][i]->GetRMS());
			OM.qOut.insert("evis",m);
		}
		drawSimulHistos(hToPlot);
		OM.printCanvas("Evis_Type_"+itos(t));
	}
	
	hToPlot.clear();
	for(unsigned int i=0; i<sps.size(); i++) {
		hCTScint[i]->SetLineColor(2+2*i);
		hCTScint[i]->Scale(1.0/t0norm[0]);
		hToPlot.push_back(hCTScint[i]);
	}
	drawSimulHistos(hToPlot);
	OM.printCanvas("Scint_Costheta");
}

int main(int argc, char *argv[]) {
	gROOT->SetStyle("Plain");
	gStyle->SetPalette(1);
	gStyle->SetNumberContours(255);
	gStyle->SetOptStat("");
	
	for(int l = 50; l <1000; l *= 2) {
		OutputManager OM("MC_Compare",getEnvSafe("UCNA_ANA_PLOTS")+"/test/MC_Compare/"+itos(l)+"_keV");
		
		PMTCalibrator PCal(15925);
		
		G4toPMT g2p(true);
		g2p.addFile("/home/mmendenhall/geant4/output/IsotLine_eGunRandMomentum_"+itos(l)+".0keV/analyzed_*");
		g2p.setCalibrator(PCal);
		
		mc_compare_plots(OM,g2p,g2p);
		OM.write();
		OM.setWriteRoot(true);
	}
	return 0;
}
