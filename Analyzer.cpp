#include <exception>
#include <TStyle.h>
#include <iostream>

#include "ControlMenu.hh"
#include "PathUtils.hh"
#include "PositionResponse.hh"
#include "PositionStudies.hh"
#include "AsymmetryAnalyzer.hh"
#include "EndpointStudy.hh"
#include "PostOfficialAnalyzer.hh"
#include "PlotMakers.hh"
#include "PMTGenerator.hh"
#include "ReSource.hh"
#include "RData.hh"
#include "G4toPMT.hh"


std::vector<RunNum> selectRuns(RunNum r0, RunNum r1, std::string typeSelect) {
	char tmp[1024];
	if(typeSelect=="ref") {
		sprintf(tmp,"run_number >= %i AND run_number <= %i AND run_type = 'SourceCalib' ORDER BY run_number ASC",r0,r1);
		std::vector<RunNum> sruns = CalDBSQL::getCDB()->findRuns(tmp);
		std::vector<RunNum> rruns;
		for(std::vector<RunNum>::iterator it = sruns.begin(); it != sruns.end(); it++)
			if(CalDBSQL::getCDB()->getGMSRun(*it) == *it)
				rruns.push_back(*it);
		return rruns;
	} else if(typeSelect=="all")
		sprintf(tmp,"run_number >= %i AND run_number <= %i ORDER BY run_number ASC",r0,r1);
	else if(typeSelect=="asym")
		sprintf(tmp,"run_number >= %i AND run_number <= %i AND run_type = 'Asymmetry' ORDER BY run_number ASC",r0,r1);
	else if(typeSelect=="LED")
		sprintf(tmp,"run_number >= %i AND run_number <= %i AND run_type = 'LEDCalib' ORDER BY run_number ASC",r0,r1);
	else if(typeSelect=="source")
		sprintf(tmp,"run_number >= %i AND run_number <= %i AND run_type = 'SourceCalib' ORDER BY run_number ASC",r0,r1);
	else if(typeSelect=="beta")
		sprintf(tmp,"run_number >= %i AND run_number <= %i AND run_type = 'Asymmetry' AND gate_valve = 'Open' ORDER BY run_number ASC",r0,r1);
	else if(typeSelect=="bg")
		sprintf(tmp,"run_number >= %i AND run_number <= %i AND run_type = 'Asymmetry' AND gate_valve = 'Closed' ORDER BY run_number ASC",r0,r1);
	else
		sprintf(tmp,"0 = 1");
	return CalDBSQL::getCDB()->findRuns(tmp);
}

void mi_VerifyCalperiods(std::deque<std::string>&, std::stack<std::string>&) {
	std::vector<RunNum> C = CalDBSQL::getCDB()->findRuns("1 ORDER BY run_number ASC");
	printf("Runs DB contains %i runs...\n",(int)C.size());
	for(std::vector<RunNum>::iterator it=C.begin(); it!=C.end(); it++) {
		RunNum gmsrun = CalDBSQL::getCDB()->getGMSRun(*it);
		if(!gmsrun)
			printf("*** WARNING: No calibrations found for run %i!\n",gmsrun);
	}
}

void mi_EndpointStudy(std::deque<std::string>&, std::stack<std::string>& stack) {
	unsigned int nr = streamInteractor::popInt(stack);
	RunNum r1 = streamInteractor::popInt(stack);
	RunNum r0 = streamInteractor::popInt(stack);
	process_xenon(r0,r1,nr);
}

void mi_EndpointStudySim(std::deque<std::string>&, std::stack<std::string>& stack) {
	RunNum r1 = streamInteractor::popInt(stack);
	RunNum r0 = streamInteractor::popInt(stack);
	simulate_xenon(r0,r1);
}

void mi_EndpointEnresPlot(std::deque<std::string>&, std::stack<std::string>&) {
	assert(false); //TODO
	//RunNum r0 = streamInteractor::popInt(stack);
	//energyResolution(r0);
}

void mi_PosmapPlot(std::deque<std::string>&, std::stack<std::string>& stack) {
	unsigned int pmid = streamInteractor::popInt(stack);
	if(CalDBSQL::getCDB()->isValid(13883)) {
		OutputManager OM("Foo",getEnvSafe("UCNA_ANA_PLOTS")+"/Eta/Posmap_"+itos(pmid));
		etaPlot(OM,CalDBSQL::getCDB()->getPositioningCorrectorByID(pmid),pmid<1000,pmid<1000?2.5:250);
	} else {
		printf("Invalid CalDB!\n");
	}
}

void mi_nPEPlot(std::deque<std::string>&, std::stack<std::string>& stack) {
	RunNum rn = streamInteractor::popInt(stack);
	PMTCalibrator PCal(rn,CalDBSQL::getCDB());
	OutputManager OM("Foo",getEnvSafe("UCNA_ANA_PLOTS")+"/nPE/Posmap_"+itos(rn));
	npePlot(OM,&PCal);
}

void mi_PostprocessSources(std::deque<std::string>&, std::stack<std::string>& stack) {
	RunNum r1 = streamInteractor::popInt(stack);
	RunNum r0 = streamInteractor::popInt(stack);
	std::vector<RunNum> C = selectRuns(r0,r1,"source");
	if(!C.size()) {
		printf("No source runs found in Analysis DB; attempting manual scan...\n");
		for(RunNum r = r0; r<= r1; r++)
			reSource(r);
		return;
	}
	printf("Found %i source runs...\n",(int)C.size());
	for(std::vector<RunNum>::iterator it=C.begin(); it!=C.end(); it++)
		reSource(*it);
}

void mi_PlotGMS(std::deque<std::string>&, std::stack<std::string>& stack) {
	std::string typeSelect = streamInteractor::popString(stack);
	RunNum r1 = streamInteractor::popInt(stack);
	RunNum r0 = streamInteractor::popInt(stack);
	plotGMScorrections(selectRuns(r0,r1,typeSelect));
}

void mi_dumpPosmap(std::deque<std::string>&, std::stack<std::string>& stack) {
	int pnum = streamInteractor::popInt(stack);
	dumpPosmap(getEnvSafe("UCNA_ANA_PLOTS")+"/PosmapDump/",pnum);
}

void mi_processOctet(std::deque<std::string>&, std::stack<std::string>& stack) {
	int octn = streamInteractor::popInt(stack);
	const std::string outputDir="OctetAsym_Offic";
	//const std::string outputDir="OctetAsym_10keV_Bins";
	AsymmetryAnalyzer::processedLocation = getEnvSafe("UCNA_ANALYSIS_OUTPUT_DIR")+"/"+outputDir+"/"+outputDir;
	
	if(octn==1000) {
		OutputManager OM("ThisNameIsNotUsedAnywhere",getEnvSafe("UCNA_ANALYSIS_OUTPUT_DIR"));
		AsymmetryAnalyzer AA(&OM,outputDir);
		processOctets(AA,Octet::loadOctets(QFile(getEnvSafe("UCNA_OCTET_LIST"))),365*24*3600);
	} else if(octn==-1000) {
		OutputManager OM("ThisNameIsNotUsedAnywhere",getEnvSafe("UCNA_ANALYSIS_OUTPUT_DIR"));
		AsymmetryAnalyzer AA_Sim(&OM,outputDir+"_Simulated",AsymmetryAnalyzer::processedLocation);
		G4toPMT simData;
		simData.addFile("/home/mmendenhall/geant4/output/Livermore_neutronBetaUnpol_geomC/analyzed_*.root");
		simuClone(getEnvSafe("UCNA_ANALYSIS_OUTPUT_DIR")+"/"+outputDir, AA_Sim, simData, 1.0, 365*24*3600);
	} else if(octn < 0) {
		Octet oct = Octet::loadOctet(QFile(getEnvSafe("UCNA_OCTET_LIST")),-octn-1);
		if(!oct.getNRuns()) return;
		OutputManager OM("ThisNameIsNotUsedAnywhere",getEnvSafe("UCNA_ANALYSIS_OUTPUT_DIR")+"/"+outputDir+"_Simulated");
		AsymmetryAnalyzer AA_Sim(&OM,oct.octName(),getEnvSafe("UCNA_ANALYSIS_OUTPUT_DIR")+"/"+outputDir+"/"+oct.octName()+"/"+oct.octName());		
		G4toPMT simData;
		simData.addFile("/home/mmendenhall/geant4/output/Livermore_neutronBetaUnpol_geomC/analyzed_*.root");
		simuClone(getEnvSafe("UCNA_ANALYSIS_OUTPUT_DIR")+"/"+outputDir+"/"+oct.octName(), AA_Sim, simData, 1.0, 24*3600);
	} else {
		Octet oct = Octet::loadOctet(QFile(getEnvSafe("UCNA_OCTET_LIST")),octn);
		if(!oct.getNRuns()) return;
		OutputManager OM("ThisNameIsNotUsedAnywhere",getEnvSafe("UCNA_ANALYSIS_OUTPUT_DIR")+"/"+outputDir);
		AsymmetryAnalyzer AA(&OM,oct.octName());
		processOctets(AA,oct.getSubdivs(oct.divlevel+1,false));
	}
}

void simulations_evis() {
	OutputManager OM("Evis2ETrue",getEnvSafe("UCNA_ANA_PLOTS")+"/Evis2ETrue/Livermore/");
	G4toPMT g2p;
	g2p.addFile("/home/mmendenhall/geant4/output/Livermore_neutronBetaUnpol_geomC/analyzed_*.root");
	PMTCalibrator PCal(16000,CalDBSQL::getCDB());
	g2p.setCalibrator(PCal);
	SimSpectrumInfo(g2p,OM);
	OM.setWriteRoot(true);
	OM.write();
}

void mi_Special(std::deque<std::string>&, std::stack<std::string>&) {
	
	uploadRunSources();
	return;
	
	simulations_evis();
	return;
	
	makeCorrectionsFile(getEnvSafe("UCNA_ANA_PLOTS")+"/SpectrumCorrection.txt");
	return;	
	
	if(1) {
		RunNum rmin = 15991; //14264;
		RunNum rmax = 16077; //14347;
		OutputManager OM("AnodeCal",std::string("../PostPlots/AnodeCal_")+itos(rmin)+"_"+itos(rmax));
		PostOfficialAnalyzer POA(true);
		for(unsigned int i=rmin; i<=rmax; i++)
			POA.addRun(i);
		AnodeCalibration(POA,OM,12);
		return;
	}
	
	if(0) {
		OutputManager OM("SimAnodeCal","../PostPlots/SimAnodeCal/");
		G4toPMT g2p;
		g2p.addFile("/home/mmendenhall/geant4/output/Baseline_20110914_uniformRandMomentum_geomC/analyzed_*.root");
		
		PMTCalibrator PCal(16000,CalDBSQL::getCDB());
		g2p.setCalibrator(PCal);
		
		AnodeCalibration(g2p,OM,1);
		return;
	}
	
	if(0) {
		OutputManager OM("CathodeCal","../PostPlots/CathodeCal_14269_14270/");
		PostOfficialAnalyzer POA(true);
		for(unsigned int i=14264; i<=14269; i++)
			POA.addRun(i);
		CathodeCalibration(POA,OM);
		return;
	}
}

void Analyzer(std::deque<std::string> args=std::deque<std::string>()) {
	
	gStyle->SetPalette(1);
	gStyle->SetNumberContours(255);
	gStyle->SetOptStat("e");	
	TCanvas defaultCanvas;
	defaultCanvas.SetFillColor(0);
	defaultCanvas.SetCanvasSize(300,300);
		
	inputRequester exitMenu("Exit Menu",&menutils_Exit);
	inputRequester peek("Show stack",&menutils_PrintStack);
	
	// selection utilities
	NameSelector selectRuntype("Run Type");
	selectRuntype.addChoice("All Runs","all");
	selectRuntype.addChoice("LED Runs","LED");
	selectRuntype.addChoice("Source Runs","source");
	selectRuntype.addChoice("Beta & BG asymmetry runs","asym");
	selectRuntype.addChoice("Beta Runs","beta");
	selectRuntype.addChoice("Background Runs","bg");
	selectRuntype.addChoice("GMS Reference Runs","ref");
	selectRuntype.setDefault("all");
		
	// postprocessing/plots routines
	inputRequester pm_posmap("Generate Position Map",&mi_EndpointStudy);
	pm_posmap.addArg("Start Run");
	pm_posmap.addArg("End Run");
	pm_posmap.addArg("n Rings","12");
	inputRequester pm_posmap_sim("Simulate Position Map",&mi_EndpointStudySim);
	pm_posmap_sim.addArg("Start Run");
	pm_posmap_sim.addArg("End Run");
	inputRequester pm_mi2("Energy Resolution Plots",&mi_EndpointEnresPlot);
	pm_mi2.addArg("Run Number");
	inputRequester pm_mi5("Verify calibration assignments",&mi_VerifyCalperiods);
	
	inputRequester plotGMS("Plot GMS corrections",&mi_PlotGMS);
	plotGMS.addArg("Start Run");
	plotGMS.addArg("End Run");
	plotGMS.addArg("","",&selectRuntype);
	
	inputRequester posmapPlot("Plot Position Map",&mi_PosmapPlot);
	posmapPlot.addArg("Posmap ID");
	
	inputRequester nPEPlot("Plot nPE/MeV",&mi_nPEPlot);
	nPEPlot.addArg("Run Number");
	
	inputRequester posmapDumper("Dump Posmap",&mi_dumpPosmap);
	posmapDumper.addArg("Posmap ID");
	
	inputRequester octetProcessor("Process Octet",&mi_processOctet);
	octetProcessor.addArg("Octet number");
	
	inputRequester specialJunk("Special Junk",&mi_Special);
	
	// Posprocessing menu
	OptionsMenu PostRoutines("Postprocessing Routines");
	PostRoutines.addChoice(&pm_posmap,"pmap");
	PostRoutines.addChoice(&pm_posmap_sim,"pmsim");
	PostRoutines.addChoice(&pm_mi2);	
	PostRoutines.addChoice(&pm_mi5);
	PostRoutines.addChoice(&plotGMS);
	PostRoutines.addChoice(&posmapPlot);
	PostRoutines.addChoice(&nPEPlot);
	PostRoutines.addChoice(&posmapDumper,"dpm");
	PostRoutines.addChoice(&octetProcessor,"oct");
	PostRoutines.addChoice(&specialJunk,"spec");
	PostRoutines.addChoice(&exitMenu,"x");	
	
	// special run processing
	inputRequester postSources("Reprocess Source Runs",&mi_PostprocessSources);
	postSources.addArg("Start Run");
	postSources.addArg("End Run");
	
	// main menu
	OptionsMenu OM("Analyzer Main Menu");
	OM.addChoice(&PostRoutines,"pr");
	OM.addChoice(&postSources,"sr");
	OM.addChoice(&exitMenu,"x");
	OM.addSynonym("x","exit");
	OM.addSynonym("x","quit");
	OM.addSynonym("x","bye");
	OM.addChoice(&peek,"peek",SELECTOR_HIDDEN);
	
	std::stack<std::string> stack;
	OM.doIt(args,stack);
	
	printf("\n\n\n>>>>> Goodbye. <<<<<\n\n\n");
}

int main(int argc, char *argv[]) {
	std::deque<std::string> args;
	for(int i=1; i<argc; i++)
		args.push_back(argv[i]);
	Analyzer(args);
	return 0;
}
