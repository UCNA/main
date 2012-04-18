#include "FierzOctetAnalyzer.hh"
#include "GraphicsUtils.hh"
#include <TStyle.h>

FierzOctetAnalyzer::FierzOctetAnalyzer(OutputManager* pnt, const string& nm, const string& inflname): OctetAnalyzer(pnt,nm,inflname) {
	// set up histograms of interest
	for(Side s = EAST; s <= WEST; ++s)
		qFullEnergySpectrum[s] = registerCoreHist("hFullEnergy", "Full Energy", 100, 0, 1000, s, &hFullEnergySpectrum[s]);
}

void FierzOctetAnalyzer::fillCoreHists(ProcessedDataScanner& PDS, double weight) {
	// fill events spectrum for Type 0 beta events on each side
	if(!(PDS.fSide==EAST || PDS.fSide==WEST)) return;
	if(PDS.fType == TYPE_0_EVENT && PDS.fPID == PID_BETA)
		//hFullEnergySpectrum[PDS.fSide]->Fill(PDS.scints[PDS.fSide].energy.x, weight);
		hFullEnergySpectrum[PDS.fSide]->Fill(PDS.getEtrue(), weight);
}

void FierzOctetAnalyzer::calculateResults() {
	// form (blinded) super-ratio and super-sum of anode spectra
	hFullEnergySR = (TH1F*)calculateSR("Full_Energy_Asymmetry", qFullEnergySpectrum[EAST], qFullEnergySpectrum[WEST]);
	hFullEnergySR->SetMinimum(-0.20);
	hFullEnergySR->SetMaximum(0.0);
	hFullEnergySS = (TH1F*)calculateSuperSum("Full_Energy_SuperSum", qFullEnergySpectrum[EAST], qFullEnergySpectrum[WEST]);
}

void FierzOctetAnalyzer::makePlots() {
	// draw the asymmetry to "Full_Energy_Asymmetry.pdf" in this analyzer's base directory
	hFullEnergySR->Draw();
	printCanvas("Full_Energy_Asymmetry");

	// also draw the Super Sum
	hFullEnergySS->Draw();
	printCanvas("Full_Energy_SuperSum");

	// and draw the raw spectra (with both sides / AFP states in same plot), in their own subfolder "MWPC_Energy"
	drawQuadSides(qFullEnergySpectrum[EAST],qFullEnergySpectrum[WEST],true,"MWPC_Energy");
	//drawQuadSides(qFullEnergySpectrum[EAST],qFullEnergySpectrum[WEST],true,"E_vis");
}

void FierzOctetAnalyzer::compareMCtoData(RunAccumulator& OAdata, float simfactor) {
	// re-cast to correct type
	FierzOctetAnalyzer& dat = (FierzOctetAnalyzer&)OAdata;

	// scale Super-Sum by simulation factor draw comparison with data (red=data, blue=MC)
	hFullEnergySS->Scale(simfactor);
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
		
	// environment variable "UCNA_ANALYSIS_OUTPUT_DIR" needs to be set to an ouput directory where you keep analysis results
	OutputManager OM("ThisNameIsNotUsedAnywhere",getEnvSafe("UCNA_ANALYSIS_OUTPUT_DIR"));

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
		processOctets(OAE,octs,3600);
		
	} else if(string(argv[1])=="sim") {
		
		// make another version of the analyzer, with a different output path for simulation, pointing to the data to clone
		FierzOctetAnalyzer OAE_Sim(&OM,"Full_Energy_Asymmetry_Example_Simulated",FierzOctetAnalyzer::processedLocation);

		// and set up a simulation data source from neutron beta decay simulation
		G4toPMT simData;
		simData.addFile("/home/mmendenhall/geant4/output/LivPhys_495_neutronBetaUnpol_geomC/analyzed_*.root");

		// point the simulation cloner to the real data; clone equal amounts of counts where not already simulated within the past hour
		simuClone(OM.basePath+"/Fierz_Processed_Analysis", OAE_Sim, simData, 1.0, 3600);		
	}
	
	return 0;
}
