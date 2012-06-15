#include "FierzOctetAnalyzer.hh"
#include "GraphicsUtils.hh"
#include "PathUtils.hh"
#include <TStyle.h>
#include <TRandom3.h>

TRandom3 fierz_rnd_src;


FierzOctetAnalyzer::FierzOctetAnalyzer(OutputManager* pnt, const string& nm, const string& inflname): OctetAnalyzer(pnt,nm,inflname) {
	for(Side s = EAST; s <= WEST; ++s) {
		// energy histograms
		qFullEnergySpectrum[s] = registerCoreHist("hFullEnergy", "Full Energy", 100, 0, 1000, s, &hFullEnergySpectrum[s]);
		// trigger threshold counts
		pTriggerThreshold[s][0] = registerFGBGPair("hTriggerAll", "Trigger threshold, all events",150,0,300,AFP_OTHER,s);
		pTriggerThreshold[s][1] = registerFGBGPair("hTriggerTrig", "Trigger threshold, triggered events",150,0,300,AFP_OTHER,s);
	}
}

void FierzOctetAnalyzer::fillCoreHists(ProcessedDataScanner& PDS, double weight) {
	Side s = PDS.fSide;
	if(!(s==EAST || s==WEST)) return;
	if(PDS.fType == TYPE_0_EVENT && PDS.fPID == PID_BETA) {
		// fill events spectrum with Etrue on both sides
		hFullEnergySpectrum[s]->Fill(PDS.getEtrue(), weight);
		// trigger threshold extraction histograms
		if(currentGV==GV_OPEN) {
			if(PGen.getCalibrator() != PDS.ActiveCal) PGen.setCalibrator(PDS.ActiveCal);
			PGen.setSide(s);
			PGen.setPosition(PDS.wires[s][X_DIRECTION].center,PDS.wires[s][Y_DIRECTION].center);
			double evis = PGen.generate(fierz_rnd_src.Uniform(350)).energy.x;
			double etrue = PDS.ActiveCal->Etrue(s,TYPE_0_EVENT,s==EAST?evis:0,s==WEST?evis:0);
			pTriggerThreshold[s][0].h[GV_OPEN]->Fill(etrue);
			pTriggerThreshold[s][1].h[GV_OPEN]->Fill(etrue,PGen.triggered());
		}
	}
}


void FierzOctetAnalyzer::calculateResults() {
	// form (blinded) super-ratio and super-sum of energy spectra
	hFullEnergySR = (TH1F*)calculateSR("Full_Energy_Asymmetry", qFullEnergySpectrum[EAST], qFullEnergySpectrum[WEST]);
	hFullEnergySR->SetMinimum(-0.20);
	hFullEnergySR->SetMaximum(0.0);
	hFullEnergySS = (TH1F*)calculateSuperSum("Full_Energy_SuperSum", qFullEnergySpectrum[EAST], qFullEnergySpectrum[WEST]);
	// Calculate trigger efficiency fraction
	for(Side s = EAST; s <= WEST; ++s) {
		gTrigCurve[s] = (TGraphAsymmErrors*)addObject(new TGraphAsymmErrors(pTriggerThreshold[s][0].h[GV_OPEN]->GetNbinsX()));
		gTrigCurve[s]->SetName(sideSubst("gTrigCurve_%c",s).c_str());
		gTrigCurve[s]->BayesDivide(pTriggerThreshold[s][1].h[GV_OPEN],pTriggerThreshold[s][0].h[GV_OPEN],"w");
	}
}

void FierzOctetAnalyzer::makePlots() {
	// draw the asymmetry to "Full_Energy_Asymmetry.pdf" in this analyzer's base directory
	hFullEnergySR->Draw();
	printCanvas("Full_Energy_Asymmetry");

	// also draw the Super Sum to "Full_Energy_SuperSum.pdf" in the same dir
	hFullEnergySS->Draw();
	printCanvas("Full_Energy_SuperSum");

	// and draw the raw spectra (with both sides / AFP states in same plot), in their own subfolder "MWPC_Energy"
	drawQuadSides(qFullEnergySpectrum[EAST], qFullEnergySpectrum[WEST], true, "Full_Energy");
	
	// draw the trigger efficiency curves in their own subfolder
	for(Side s = EAST; s <= WEST; ++s) {
		pTriggerThreshold[s][0].h[GV_OPEN]->SetLineColor(4);
		pTriggerThreshold[s][0].h[GV_OPEN]->Draw();
		pTriggerThreshold[s][1].h[GV_OPEN]->SetLineColor(2);
		pTriggerThreshold[s][1].h[GV_OPEN]->Draw("Same");
		printCanvas(sideSubst("TrigEff/Input_%c",s));

		gTrigCurve[s]->SetMinimum(-0.10);
		gTrigCurve[s]->SetMaximum(1.10);
		gTrigCurve[s]->Draw("AP");
		gTrigCurve[s]->SetTitle(sideSubst("%s Trigger Efficiency",s).c_str());
		gTrigCurve[s]->GetXaxis()->SetTitle("Energy [keV]");
		gTrigCurve[s]->GetXaxis()->SetLimits(0,300);
		gTrigCurve[s]->GetYaxis()->SetTitle("Efficiency");
		gTrigCurve[s]->Draw("AP");
		printCanvas(sideSubst("TrigEff/TrigEff_%c",s));
	}
}

void FierzOctetAnalyzer::compareMCtoData(RunAccumulator& OAdata) {
	// re-cast to correct type
	FierzOctetAnalyzer& dat = (FierzOctetAnalyzer&)OAdata;

	// draw comparison with data (red=data, blue=MC)
	drawHistoPair(dat.hFullEnergySS,hFullEnergySS);
	printCanvas("DataComparison/Full_Energy_SuperSum");

	// same for Super-Ratio, except normalization is not needed
	drawHistoPair(dat.hFullEnergySR,hFullEnergySR);
	printCanvas("DataComparison/Full_Energy_Asymmetry");
}

string FierzOctetAnalyzer::processedLocation = "";	// set this later depending on situtation

int main(int argc, char *argv[]) {
	
	if(argc != 2) {
		printf("To scan data, use:\n %s scan\n",argv[0]);
		printf("To simulate scanned data, use:\n %s sim\n",argv[0]);
		exit(0);
	}
	
	gStyle->SetOptStat("e");	
		
	// environment variable "UCNA_ANA_PLOTS" needs to be set to an ouput directory where you keep analysis results
	OutputManager OM("ThisNameIsNotUsedAnywhere",getEnvSafe("UCNA_ANA_PLOTS"));

	// after running once, we can use the results for errorbar estimation in low-count background bins on later scans
	FierzOctetAnalyzer::processedLocation = OM.basePath+"/Fierz_Processed_Analysis/Fierz_Processed_Analysis";

	if(string(argv[1])=="scan") {
		
		// results will be placed in "Fierz_Processed_Analysis/" subfolder of base directory
		FierzOctetAnalyzer OAE(&OM,"Fierz_Processed_Analysis");
		
		// get list of octets to process, available in Aux/OctetList_2010.txt
		vector<Octet> octs = Octet::loadOctets(QFile("Aux/OctetList_2010.txt"));
		
		// for this example, limit to first three octets
		while(octs.size() > 3)
			octs.pop_back();
		
		// process the octets; skip re-processing if already produced under one hour ago
		processOctets(OAE, octs, 3600);
		
	} else if(string(argv[1])=="sim") {
		
		// make another version of the analyzer, with a different output path for simulation, pointing to the data to clone
		//FierzOctetAnalyzer OAE_Sim(&OM, "Full_Energy_Asymmetry_Example_Simulated", FierzOctetAnalyzer::processedLocation);
		FierzOctetAnalyzer OAE_Sim(&OM,"Fierz_Processed_Simulated");

		// and set up a simulation data source from neutron beta decay simulation
		G4toPMT simData;
		simData.addFile("/home/mmendenhall/geant4/output/LivPhys_495_neutronBetaUnpol_geomC/analyzed_*.root");

		// point the simulation cloner to the real data; clone equal amounts of counts where not already simulated within the past hour
		OAE_Sim.simuClone(OM.basePath+"/Fierz_Processed_Analysis", simData, 1.0, 3600);		
	}
	
	return 0;
}
