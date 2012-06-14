#include "G4toPMT.hh"
#include "PenelopeToPMT.hh"
#include "OutputManager.hh"
#include "PathUtils.hh"
#include "GraphicsUtils.hh"
#include <TROOT.h>
#include <TStyle.h>

void mc_compare_plots(OutputManager& OM, Sim2PMT& SP1, Sim2PMT& SP2, double emax) {
	std::vector<Sim2PMT*> sps;
	sps.push_back(&SP1);
	sps.push_back(&SP2);
	
	std::vector<TH1F*> hEvis[TYPE_II_EVENT+1];
	std::vector<TH1F*> hEquench[TYPE_II_EVENT+1];
	std::vector<TH1F*> hCTScint;
	std::vector<TH1F*> hCTMWPC;
	std::vector<double> t0norm;
	
	for(unsigned int i=0; i<sps.size(); i++) {
		Sim2PMT* SP = sps[i];
		
		// define histograms
		for(EventType t=TYPE_0_EVENT; t <= TYPE_II_EVENT; ++t) {
			hEvis[t].push_back(OM.registeredTH1F("hEvis_"+itos(t)+"_"+itos(i),
												 "Type "+itos(t)+" Energy Spectrum",
												 100,0,emax));
			hEquench[t].push_back(OM.registeredTH1F("hEquench_"+itos(t)+"_"+itos(i),
												 "Type "+itos(t)+" Quenched Energy Spectrum",
												 100,0,emax));
		}
		hCTScint.push_back(OM.registeredTH1F("hCTScint_"+itos(i),"cos(theta) entering scintillator",100,0.,1.0));
		hCTMWPC.push_back(OM.registeredTH1F("hCTMWPC_"+itos(i),"cos(theta) entering MWPC",100,0.,1.0));
		
		// fill histograms
		SP->startScan();
		while(SP->nextPoint()) {
			if(SP->fPID != PID_BETA) continue;
			if(SP->fType <= TYPE_II_EVENT) {
				hEvis[SP->fType].back()->Fill(SP->eDep[EAST]+SP->eDep[WEST],SP->physicsWeight);
				hEquench[SP->fType].back()->Fill(SP->eQ[EAST]+SP->eQ[WEST],SP->physicsWeight);
			}
			if(SP->fSide <= WEST) {
				hCTScint.back()->Fill(SP->cosThetaInScint[SP->fSide],SP->physicsWeight);
				hCTMWPC.back()->Fill(SP->cosThetaInWinIn[SP->fSide],SP->physicsWeight);
			}
		}
		
		// record normalization counts
		t0norm.push_back(hEvis[TYPE_0_EVENT].back()->Integral());
	}
	
	/////////////
	// process data and make plots
	/////////////
	
	// normaliztion factors to Type 0 counts
	for(unsigned int i=0; i<t0norm.size(); i++)
		t0norm[i] = 1.0/t0norm[i];
	
	std::vector<TH1*> hToPlot;
	OM.defaultCanvas->cd();
	
	gStyle->SetOptStat("emr");
	for(EventType t=TYPE_0_EVENT; t <= TYPE_II_EVENT; ++t) {
		hToPlot.clear();
		for(unsigned int i=0; i<sps.size(); i++) {
			Stringmap m;
			m.insert("type",t);
			m.insert("sim",i);
			m.insert("mean",hEvis[t][i]->GetMean());
			m.insert("rms",hEvis[t][i]->GetRMS());
			m.insert("counts",hEvis[t][i]->Integral());
			m.insert("normcounts",hEvis[t][i]->Integral()*t0norm[i]);
			hEvis[t][i]->SetLineColor(2+2*i);
			hEvis[t][i]->Scale(t0norm[i]/hEvis[t][i]->GetBinWidth(1));
			hToPlot.push_back(hEvis[t][i]);
			OM.qOut.insert("evis",m);
		}
		drawSimulHistos(hToPlot);
		OM.printCanvas("Evis_Type_"+itos(t));
		
		hToPlot.clear();
		for(unsigned int i=0; i<sps.size(); i++) {
			Stringmap m;
			m.insert("type",t);
			m.insert("sim",i);
			m.insert("mean",hEquench[t][i]->GetMean());
			m.insert("rms",hEquench[t][i]->GetRMS());
			m.insert("counts",hEquench[t][i]->Integral());
			m.insert("normcounts",hEquench[t][i]->Integral()*t0norm[i]);
			hEquench[t][i]->SetLineColor(2+2*i);
			hEquench[t][i]->Scale(t0norm[i]/hEquench[t][i]->GetBinWidth(1));
			hToPlot.push_back(hEquench[t][i]);
			OM.qOut.insert("equench",m);
		}
		drawSimulHistos(hToPlot);
		OM.printCanvas("Equench_Type_"+itos(t));
	}
	gStyle->SetOptStat("");
	
	hToPlot.clear();
	for(unsigned int i=0; i<sps.size(); i++) {
		hCTScint[i]->SetLineColor(2+2*i);
		hCTScint[i]->Scale(t0norm[i]/hCTScint[i]->GetBinWidth(1));
		hToPlot.push_back(hCTScint[i]);
	}
	drawSimulHistos(hToPlot);
	OM.printCanvas("Scint_Costheta");
	
	hToPlot.clear();
	for(unsigned int i=0; i<sps.size(); i++) {
		hCTMWPC[i]->SetLineColor(2+2*i);
		hCTMWPC[i]->Scale(t0norm[i]/hCTMWPC[i]->GetBinWidth(1));
		hToPlot.push_back(hCTMWPC[i]);
	}
	drawSimulHistos(hToPlot);
	OM.printCanvas("MWPC_In_Costheta");
}

int main(int argc, char *argv[]) {
	gROOT->SetStyle("Plain");
	gStyle->SetPalette(1);
	gStyle->SetNumberContours(255);
	gStyle->SetOptStat("");
	
	// list of simulated energies
	int enlist[] = {100,150,200,300,400,600,800};
	
	for(int i = 0; i < 7; i++) {
		int l = enlist[i];
		
		OutputManager OM("MC_Compare",getEnvSafe("UCNA_ANA_PLOTS")+"/test/MC_Compare/"+itos(l)+"_keV");
		Stringmap mcdat;
		mcdat.insert("energy",l);
		mcdat.insert("n_MC",2);
		
		PMTCalibrator PCal(15925);
		
		G4toPMT g2p(true);
		std::string fname = "/home/mmendenhall/geant4/output/IsotLine_eGunRandMomentum_"+itos(l)+".0keV/analyzed_*";
		mcdat.insert("MC_1",fname);
		g2p.addFile(fname);
		if(!g2p.getnFiles()) continue;
		g2p.setCalibrator(PCal);
		
		PenelopeToPMT p2p;
		fname = "/home/ucna/penelope_output/iso_line_sources/event_"+itos(l/50-1)+" _*.root";
		mcdat.insert("MC_2",fname);
		p2p.addFile(fname);
		if(!p2p.getnFiles()) continue;
		p2p.setCalibrator(PCal);
		
		mc_compare_plots(OM,g2p,p2p,2*l);
		OM.qOut.insert("mcdat",mcdat);
		OM.write();
		OM.setWriteRoot(true);
	}
	return 0;
}
