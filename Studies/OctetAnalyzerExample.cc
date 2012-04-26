#include "OctetAnalyzerExample.hh"
#include "GraphicsUtils.hh"
#include <TStyle.h>

OctetAnalyzerExample::OctetAnalyzerExample(OutputManager* pnt, const std::string& nm, const std::string& inflname): OctetAnalyzer(pnt,nm,inflname) {
	// set up histograms of interest
	for(Side s = EAST; s <= WEST; ++s)
		qAnodeSpectrum[s] = registerCoreHist("hAnode", "Wirechamber Energy", 100, 0, 20, s, &hAnodeSpectrum[s]);
}

void OctetAnalyzerExample::fillCoreHists(ProcessedDataScanner& PDS, double weight) {
	// fill wirechamber spectrum for Type 0 beta events on each side
	if(!(PDS.fSide==EAST || PDS.fSide==WEST)) return;
	if(PDS.fType == TYPE_0_EVENT && PDS.fPID == PID_BETA)
		hAnodeSpectrum[PDS.fSide]->Fill(PDS.mwpcEnergy[PDS.fSide],weight);
}

void OctetAnalyzerExample::calculateResults() {
	// form (blinded) super-ratio and super-sum of anode spectra
	hAnodeSR = (TH1F*)calculateSR("Wirechamber_Energy_Asymmetry",qAnodeSpectrum[EAST],qAnodeSpectrum[WEST]);
	hAnodeSR->SetMinimum(-0.20);
	hAnodeSR->SetMaximum(0.0);
	hAnodeSS = (TH1F*)calculateSuperSum("Wirechamber_Energy_SuperSum",qAnodeSpectrum[EAST],qAnodeSpectrum[WEST]);
}

void OctetAnalyzerExample::makePlots() {
	// draw the asymmetry to "Wirechamber_Energy_Asymmetry.pdf" in this analyzer's base directory
	hAnodeSR->Draw();
	printCanvas("Wirechamber_Energy_Asymmetry");
	// also draw the Super Sum
	hAnodeSS->Draw();
	printCanvas("Wirechamber_Energy_SuperSum");
	// and draw the raw spectra (with both sides / AFP states in same plot), in their own subfolder "MWPC_Energy"
	drawQuadSides(qAnodeSpectrum[EAST],qAnodeSpectrum[WEST],true,"MWPC_Energy");
}

void OctetAnalyzerExample::compareMCtoData(RunAccumulator& OAdata) {
	// re-cast to correct type
	OctetAnalyzerExample& dat = (OctetAnalyzerExample&)OAdata;
	// draw comparison with data (red=data, blue=MC)
	drawHistoPair(dat.hAnodeSS,hAnodeSS);
	printCanvas("DataComparison/Wirechamber_Energy_SuperSum");
	// same for Super-Ratio, except normalization is not needed
	drawHistoPair(dat.hAnodeSR,hAnodeSR);
	printCanvas("DataComparison/Wirechamber_Energy_Asymmetry");
}

std::string OctetAnalyzerExample::processedLocation = "";	// set this later depending on situtation

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
	OctetAnalyzerExample::processedLocation = OM.basePath+"/Anode_Asymmetry_Example/Anode_Asymmetry_Example";

	if(std::string(argv[1])=="scan") {
		
		// results will be placed in "Anode_Asymmetry_Example/" subfolder of base directory
		OctetAnalyzerExample OAE(&OM,"Anode_Asymmetry_Example");
		
		// get list of octets to process, available in Aux/OctetList_2010.txt
		vector<Octet> octs = Octet::loadOctets(QFile("Aux/OctetList_2010.txt"));
		
		// for this example, limit to first three octets
		while(octs.size() > 3)
			octs.pop_back();
		
		// process the octets; skip re-processing if already produced under one hour ago
		processOctets(OAE,octs,3600);
		
	} else if(std::string(argv[1])=="sim") {
		
		// make another version of the analyzer, with a different output path for simulation
		OctetAnalyzerExample OAE_Sim(&OM,"Anode_Asymmetry_Example_Simulated");

		// and set up a simulation data source from neutron beta decay simulation
		G4toPMT simData;
		simData.addFile("/home/mmendenhall/geant4/output/LivPhys_495_neutronBetaUnpol_geomC/analyzed_*.root");
		// point the simulation cloner to the real data; clone equal amounts of counts where not already simulated within the past hour
		OAE_Sim.simuClone(OM.basePath+"/Anode_Asymmetry_Example", simData, 1.0, 3600);
		
	}
	
	return 0;
}
