#include <exception>
#include <TStyle.h>
#include <iostream>

#include "ControlMenu.hh"
#include "PathUtils.hh"
#include "GraphicsUtils.hh"
#include "PositionResponse.hh"
#include "WirechamberStudy.hh"
#include "AsymmetryAnalyzer.hh"
#include "EndpointStudy.hh"
#include "PostOfficialAnalyzer.hh"
#include "PlotMakers.hh"
#include "PMTGenerator.hh"
#include "ReSource.hh"
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

void mi_EndpointStudy(std::deque<std::string>&, std::stack<std::string>& stack) {
	unsigned int nr = streamInteractor::popInt(stack);
	RunNum r1 = streamInteractor::popInt(stack);
	RunNum r0 = streamInteractor::popInt(stack);
	process_xenon(r0,r1,nr);
}

void mi_EndpointStudySim(std::deque<std::string>&, std::stack<std::string>& stack) {
	RunNum rsingle = streamInteractor::popInt(stack);
	RunNum r1 = streamInteractor::popInt(stack);
	RunNum r0 = streamInteractor::popInt(stack);
	simulate_xenon(r0,r1,rsingle);
}

void mi_WirechamberStudy(std::deque<std::string>&, std::stack<std::string>& stack) {
	unsigned int nr = streamInteractor::popInt(stack);
	RunNum r1 = streamInteractor::popInt(stack);
	RunNum r0 = streamInteractor::popInt(stack);
	runWirechamberAnalyzer(r0,r1,nr);
}

void mi_PosmapPlot(std::deque<std::string>&, std::stack<std::string>& stack) {
	unsigned int pmid = streamInteractor::popInt(stack);
	if(CalDBSQL::getCDB()->isValid(13883)) {
		OutputManager OM("Foo",getEnvSafe("UCNA_ANA_PLOTS")+"/PositionMaps/Posmap_"+itos(pmid));
		PosPlotter PP(&OM);
		PP.etaPlot(CalDBSQL::getCDB()->getPositioningCorrectorByID(pmid));
	} else {
		printf("Invalid CalDB!\n");
	}
}

void mi_nPEPlot(std::deque<std::string>&, std::stack<std::string>& stack) {
	RunNum rn = streamInteractor::popInt(stack);
	PMTCalibrator PCal(rn,CalDBSQL::getCDB());
	OutputManager OM("NPE",getEnvSafe("UCNA_ANA_PLOTS")+"/nPE/Run_"+itos(rn));
	PosPlotter PP(&OM);
	PP.npePlot(&PCal);
	OM.write();
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

void mi_delPosmap(std::deque<std::string>&, std::stack<std::string>& stack) {
	int pnum = streamInteractor::popInt(stack);
	CalDBSQL::getCDB(false)->deletePosmap(pnum);
}

void mi_listPosmaps(std::deque<std::string>&, std::stack<std::string>&) { CalDBSQL::getCDB()->listPosmaps(); }

void mi_processOctet(std::deque<std::string>&, std::stack<std::string>& stack) {
	int octn = streamInteractor::popInt(stack);
	const std::string outputDir="OctetAsym_Offic";
	AsymmetryAnalyzer::processedLocation = getEnvSafe("UCNA_ANA_PLOTS")+"/"+outputDir+"/"+outputDir;
	
	if(octn==1000) {
		OutputManager OM("ThisNameIsNotUsedAnywhere",getEnvSafe("UCNA_ANA_PLOTS"));
		AsymmetryAnalyzer AA(&OM,outputDir);
		processOctets(AA,Octet::loadOctets(QFile(getEnvSafe("UCNA_OCTET_LIST"))),365*24*3600);
	} else if(octn < 0) {
		G4toPMT simData;
		simData.addFile("/home/mmendenhall/geant4/output/LivPhys_495_neutronBetaUnpol_geomC/analyzed_*.root");
		if(octn==-1000) {
			OutputManager OM("ThisNameIsNotUsedAnywhere",getEnvSafe("UCNA_ANA_PLOTS"));
			AsymmetryAnalyzer AA_Sim(&OM,outputDir+"_Simulated");
			simuClone(getEnvSafe("UCNA_ANA_PLOTS")+"/"+outputDir, AA_Sim, simData, 1.0, 365*24*3600);
		} else {
			Octet oct = Octet::loadOctet(QFile(getEnvSafe("UCNA_OCTET_LIST")),-octn-1);
			if(!oct.getNRuns()) return;
			OutputManager OM("ThisNameIsNotUsedAnywhere",getEnvSafe("UCNA_ANA_PLOTS")+"/"+outputDir+"_Simulated");
			AsymmetryAnalyzer AA_Sim(&OM,oct.octName());		
			simuClone(getEnvSafe("UCNA_ANA_PLOTS")+"/"+outputDir+"/"+oct.octName(), AA_Sim, simData, 1.0, 24*3600);
		}
	} else {
		Octet oct = Octet::loadOctet(QFile(getEnvSafe("UCNA_OCTET_LIST")),octn);
		if(!oct.getNRuns()) return;
		OutputManager OM("ThisNameIsNotUsedAnywhere",getEnvSafe("UCNA_ANA_PLOTS")+"/"+outputDir);
		AsymmetryAnalyzer AA(&OM,oct.octName());
		processOctets(AA,oct.getSubdivs(oct.divlevel+1,false));
	}
}

void mi_evis2etrue(std::deque<std::string>&, std::stack<std::string>&) {
	OutputManager OM("Evis2ETrue",getEnvSafe("UCNA_ANA_PLOTS")+"/Evis2ETrue/Livermore/");
	G4toPMT g2p;
	g2p.addFile("/home/mmendenhall/geant4/output/LivPhys_495_neutronBetaUnpol_geomC/analyzed_*.root");
	PMTCalibrator PCal(16000,CalDBSQL::getCDB());
	g2p.setCalibrator(PCal);
	SimSpectrumInfo(g2p,OM);
	OM.setWriteRoot(true);
	OM.write();
}

void mi_sourcelog(std::deque<std::string>&, std::stack<std::string>&) { uploadRunSources(); }

void mi_radcor(std::deque<std::string>&, std::stack<std::string>&) { makeCorrectionsFile(getEnvSafe("UCNA_ANA_PLOTS")+"/SpectrumCorrection.txt"); }

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
	
	// position map routines and menu
	inputRequester pm_posmap("Generate Position Map",&mi_EndpointStudy);
	pm_posmap.addArg("Start Run");
	pm_posmap.addArg("End Run");
	pm_posmap.addArg("n Rings","12");
	inputRequester pm_posmap_sim("Simulate Position Map",&mi_EndpointStudySim);
	pm_posmap_sim.addArg("Start Run");
	pm_posmap_sim.addArg("End Run");
	pm_posmap_sim.addArg("Single Run","0");
	inputRequester posmapLister("List Posmaps",&mi_listPosmaps);
	inputRequester posmapPlot("Plot Position Map",&mi_PosmapPlot);
	posmapPlot.addArg("Posmap ID");
	inputRequester posmapDumper("Dump Posmap",&mi_dumpPosmap);
	posmapDumper.addArg("Posmap ID");
	inputRequester posmapDel("Delete Posmap",&mi_delPosmap);
	posmapDel.addArg("Posmap ID");
	inputRequester nPEPlot("Plot nPE/MeV",&mi_nPEPlot);
	nPEPlot.addArg("Run Number");
	OptionsMenu PMapR("Position Map Routines");
	PMapR.addChoice(&pm_posmap,"gen");
	PMapR.addChoice(&pm_posmap_sim,"sim");
	PMapR.addChoice(&posmapLister,"ls");
	PMapR.addChoice(&posmapPlot,"plot");
	PMapR.addChoice(&posmapDumper,"dump");
	PMapR.addChoice(&posmapDel,"del");
	PMapR.addChoice(&nPEPlot,"npe");
	PMapR.addChoice(&exitMenu,"x");
	
	// postprocessing/plots routines	
	inputRequester plotGMS("Plot GMS corrections",&mi_PlotGMS);
	plotGMS.addArg("Start Run");
	plotGMS.addArg("End Run");
	plotGMS.addArg("","",&selectRuntype);
	
	inputRequester octetProcessor("Process Octet",&mi_processOctet);
	octetProcessor.addArg("Octet number");
	
	// Posprocessing menu
	OptionsMenu PostRoutines("Postprocessing Routines");
	PostRoutines.addChoice(&plotGMS);
	PostRoutines.addChoice(&octetProcessor,"oct");
	PostRoutines.addChoice(&exitMenu,"x");	
	
	// sources
	inputRequester postSources("Fit source data",&mi_PostprocessSources);
	postSources.addArg("Start Run");
	postSources.addArg("End Run");
	inputRequester uploadSources("Upload runlog sources",&mi_sourcelog);
	// evis2etrue
	inputRequester evis2etrue("Caluculate eVis->eTrue curves",&mi_evis2etrue);
	// radiative corrections
	inputRequester radcor("Make radiative corrections table",&mi_radcor);
	// wirechamber calibration
	inputRequester anawc("Gather wirechamber calibration data",&mi_WirechamberStudy);
	anawc.addArg("Start Run");
	anawc.addArg("End Run");
	anawc.addArg("n Rings","6");
	
	// main menu
	OptionsMenu OM("Analyzer Main Menu");
	OM.addChoice(&PMapR,"pmap");
	OM.addChoice(&PostRoutines,"pr");
	OM.addChoice(&postSources,"sr");
	OM.addChoice(&uploadSources,"us");
	OM.addChoice(&evis2etrue,"ev");
	OM.addChoice(&radcor,"rc");
	OM.addChoice(&anawc,"wc");
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
