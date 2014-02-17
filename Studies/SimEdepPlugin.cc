#include "SimEdepPlugin.hh"

SimEdepPlugin::SimEdepPlugin(RunAccumulator* RA): AnalyzerPlugin(RA,"simedep") {

	myA->isSimulated = true;
	int nEnergyBins = 10;
	double energyMax = 1000;
	myA->ignoreMissingHistos = true;
	
	TProfile pEdepTemplate("eDep","eDep",nEnergyBins,0,energyMax);
	
	volNames[0] = "Eprim";
	volNames[1] = "Scint";
	volNames[2] = "ScintQ";
	volNames[3] = "Foil";
	volNames[4] = "Kevlar";
	volNames[5] = "WinOut";
	volNames[6] = "WC";
	volNames[7] = "DeadWC";
	volNames[8] = "Wires";
	volNames[9] = "WinIn";
	volNames[10] = "DeadScint";
	
	for(unsigned int i=0; i<NUM_EDEP_VOLS; i++) {
		hEdepCombo[i] = (TProfile*)myA->registerSavedHist(volNames[i],pEdepTemplate);
		for(Side s = EAST; s <= WEST; ++s) {
			for(EventType tp = TYPE_0_EVENT; tp <= TYPE_III_EVENT; ++tp) {
				hEdep[s][tp][i] = (TProfile*)myA->registerSavedHist(volNames[i]+sideSubst("_%c_", s)+itos(tp), pEdepTemplate);
			}
		}
	}
	
	myA->ignoreMissingHistos = false;
}

void SimEdepPlugin::fillCoreHists(ProcessedDataScanner& PDS, double weight) {
	Sim2PMT& S2P = (Sim2PMT&)PDS;
	Side s = S2P.fSide;
	EventType tp = S2P.fType;
	if(!(s==EAST || s==WEST) || tp >= TYPE_IV_EVENT) return;
	if(!S2P.passesPositionCut(s)) return;
	double E0 = S2P.ePrim;
	
	hEdep[s][tp][0]->Fill(E0,E0,weight);
	hEdep[s][tp][1]->Fill(E0,S2P.eDep[s],weight);
	hEdep[s][tp][2]->Fill(E0,S2P.eDep[s]-S2P.eQ[s],weight);
	hEdep[s][tp][3]->Fill(E0,S2P.edepFoils[s],weight);
	hEdep[s][tp][4]->Fill(E0,S2P.edepKevlar[s],weight);
	hEdep[s][tp][5]->Fill(E0,S2P.edepWinIn[s],weight);
	hEdep[s][tp][6]->Fill(E0,S2P.eW[s],weight);
	hEdep[s][tp][7]->Fill(E0,S2P.edepDeadMWPC[s],weight);
	hEdep[s][tp][8]->Fill(E0,S2P.edepWires[s],weight);
	hEdep[s][tp][9]->Fill(E0,S2P.edepWinOut[s],weight);
	hEdep[s][tp][10]->Fill(E0,S2P.edepDeadScint[s],weight);
}

void SimEdepPlugin::calculateResults() {
	for(unsigned int i=0; i<NUM_EDEP_VOLS; i++) {
		hEdepCombo[i]->Reset();
		for(Side s = EAST; s <= WEST; ++s) {
			for(EventType tp = TYPE_0_EVENT; tp <= TYPE_0_EVENT; ++tp) {
				hEdepCombo[i]->Add(hEdep[s][tp][i]);
			}
		}
	}
}

void SimEdepPlugin::makePlots() {
	for(unsigned int i=0; i<NUM_EDEP_VOLS; i++) {
		hEdepCombo[i]->GetYaxis()->SetTitle("energy deposition [keV]");
		hEdepCombo[i]->GetXaxis()->SetTitle("primary event energy [keV]");
		hEdepCombo[i]->Draw();
		printCanvas("p_"+volNames[i]);
	}
}

void SimEdepPlugin::makeBigTable() {
	
	//int bmax = hEdepCombo[0]->GetNbinsX();
	int bmax = hEdepCombo[0]->FindBin(782);
	
	printf("\\begin{tabular}{|l||");
	for(int b = 1; b<=bmax; b++)
		printf("c|");
	printf("}\n\t\t\\hline\n");
	
	
	printf("$E_\\text{prim}$ range\t");
	std::vector<float> ltot;
	std::vector<float> ltot_ns;
	std::vector<float> ltot_ns2;
	for(int b = 1; b<=bmax; b++) {
		ltot.push_back(0);
		ltot_ns.push_back(0);
		ltot_ns2.push_back(0);
		printf("& %.0f--%.0f keV\t",hEdepCombo[0]->GetBinLowEdge(b),hEdepCombo[0]->GetBinLowEdge(b+1));
	}
	printf("\\\\ \\hline \\hline\n");
	
	
	
	for(unsigned int i=0; i<NUM_EDEP_VOLS; i++) {
		printf("%s\t",volNames[i].c_str());
		hEdepCombo[i]->SetErrorOption("s");
		for(int b = 1; b<=bmax; b++) {
			float eloss = hEdepCombo[i]->GetBinContent(b);
			float erms = hEdepCombo[i]->GetBinError(b);
			if(i != 0 && i != 2) {
				ltot[b-1] += eloss;
				if(i!=1) {
					ltot_ns[b-1] += eloss;
					ltot_ns2[b-1] += erms*erms;
				}
			}
			printf("& %.2f (%.1f)\t",eloss,erms);
		}
		printf("\\\\ \\hline\n");
	}
	
	printf("total\t");
	for(int b = 0; b<bmax; b++)
		printf("& %.2f (%.1f)\t",ltot_ns[b],sqrt(ltot_ns2[b]));
		//printf("& %.2f\t",ltot[b]);
	printf("\\\\ \\hline\n");
}


//---------------------------------------------------------------

SimEdepAnalyzer::SimEdepAnalyzer(OutputManager* pnt, const std::string& nm, const std::string& inflName): RunAccumulator(pnt,nm,inflName) {
	addPlugin(myEdep = new SimEdepPlugin(this));
}
