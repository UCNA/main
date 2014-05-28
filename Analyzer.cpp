#include <exception>
#include <iostream>

#include "StyleSetup.hh"
#include "ControlMenu.hh"
#include "PathUtils.hh"
#include "GraphicsUtils.hh"
#include "PositionResponse.hh"
#include "BetaDecayAnalyzer.hh"
#include "XenonAnalyzer.hh"
#include "PostOfficialAnalyzer.hh"
#include "PlotMakers.hh"
#include "PMTGenerator.hh"
#include "ReSource.hh"
#include "G4toPMT.hh"
#include "PenelopeToPMT.hh"
#include "NuclEvtGen.hh"
#include "AsymmetryCorrections.hh"
#include "OctetSimuCloneManager.hh"

std::vector<RunNum> selectRuns(RunNum r0, RunNum r1, std::string typeSelect) {
	if(typeSelect=="ref") {
		std::vector<RunNum> sruns = CalDBSQL::getCDB()->findRuns("run_type = 'SourceCalib'",r0,r1);
		std::vector<RunNum> rruns;
		for(std::vector<RunNum>::iterator it = sruns.begin(); it != sruns.end(); it++)
			if(CalDBSQL::getCDB()->getGMSRun(*it) == *it)
				rruns.push_back(*it);
		return rruns;
	} else if(typeSelect=="all")
		return CalDBSQL::getCDB()->findRuns("",r0,r1);
	else if(typeSelect=="asym")
		return CalDBSQL::getCDB()->findRuns("run_type = 'Asymmetry'",r0,r1);
	else if(typeSelect=="LED")
		return CalDBSQL::getCDB()->findRuns("run_type = 'LEDCalib'",r0,r1);
	else if(typeSelect=="source")
		return CalDBSQL::getCDB()->findRuns("run_type = 'SourceCalib'",r0,r1);
	else if(typeSelect=="beta")
		return CalDBSQL::getCDB()->findRuns("run_type = 'Asymmetry' AND gate_valve = 'Open'",r0,r1);
	else if(typeSelect=="bg")
		return CalDBSQL::getCDB()->findRuns("run_type = 'Asymmetry' AND gate_valve = 'Closed'",r0,r1);
	return CalDBSQL::getCDB()->findRuns("0 = 1",r0,r1);
}

void mi_EndpointStudy(StreamInteractor* S) {
	int nr = S->popInt();
	RunNum r1 = S->popInt();
	RunNum r0 = S->popInt();
	if(nr <= 0) {
		printf("%i is not a good number of rings! Canceling!",nr);
		return;
	}
	process_xenon(r0,r1,nr);
}

void mi_EndpointStudySim(StreamInteractor* S) {
	unsigned int nRings = S->popInt();
	RunNum r1 = S->popInt();
	RunNum r0 = S->popInt();
	if(nRings <= 0) {
		printf("%i is not a good number of rings! Canceling!",nRings);
		return;
	}
	if(r0==r1)
		simulate_one_xenon(r0, nRings, true);
	else
		combine_xenon_sims(r0, r1, nRings);
}

void mi_EndpointStudyReSim(StreamInteractor* S) {
	unsigned int nRings = S->popInt();
	RunNum r1 = S->popInt();
	RunNum r0 = S->popInt();
	xenon_posmap(r0,r1,nRings);
}

void mi_PosmapPlot(StreamInteractor* S) {
	unsigned int pmid = S->popInt();
	if(CalDBSQL::getCDB()->isValid(13883)) {
		OutputManager OM("Foo",getEnvSafe("UCNA_ANA_PLOTS")+"/PositionMaps/Posmap_"+itos(pmid));
		PosPlotter PP(&OM);
		PP.etaPlot(CalDBSQL::getCDB()->getPositioningCorrectorByID(pmid),0.,2.);
	} else {
		printf("Invalid CalDB!\n");
	}
}

void mi_nPEPlot(StreamInteractor* S) {
	RunNum rn = S->popInt();
	PMTCalibrator PCal(rn);
	OutputManager OM("NPE",getEnvSafe("UCNA_ANA_PLOTS")+"/nPE/Run_"+itos(rn));
	PosPlotter PP(&OM);
	PP.npePlot(&PCal);
	OM.write();
}

void mi_PostprocessSources(StreamInteractor* S) {
	RunNum r1 = S->popInt();
	RunNum r0 = S->popInt();
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

void mi_DumpCalInfo(StreamInteractor* S) {
	std::string typeSelect = S->popString();
	RunNum r1 = S->popInt();
	RunNum r0 = S->popInt();
	QFile QOut(getEnvSafe("UCNA_ANA_PLOTS")+"/test/CalDump.txt",false);
	dumpCalInfo(selectRuns(r0,r1,typeSelect),QOut);
}

void mi_dumpPosmap(StreamInteractor* S) {
	int pnum = S->popInt();
	dumpPosmap(getEnvSafe("UCNA_ANA_PLOTS")+"/PosmapDump/",pnum);
}

void mi_delPosmap(StreamInteractor* S) {
	int pnum = S->popInt();
	CalDBSQL::getCDB(false)->deletePosmap(pnum);
}

void mi_listPosmaps(StreamInteractor*) { CalDBSQL::getCDB()->listPosmaps(); }

void mi_displayOctetList(StreamInteractor*) { displayOctetList(); }

void mi_processOctet(StreamInteractor* S) {
	int octn = S->popInt();
	
	std::string octlist = getEnvSafe("UCNA_OCTET_LIST");
	char yrdigit = octlist[octlist.size()-5];
	std::string replaydirname = "Asym_2011";
	if(yrdigit == '3') replaydirname = "Asym_2012";
		
	OctetSimuCloneManager OSCM(replaydirname);
	OutputManager OM("ThisNameIsNotUsedAnywhere",getEnvSafe("UCNA_ANA_PLOTS"));
	
	// simulations input setup
	//OSCM.simFile = "/data2/mmendenhall/G4Out/2010/20120823_neutronBetaUnpol/analyzed_";
	OSCM.simFile = getEnvSafe("G4OUTDIR")+"/2011-2012geo_1mil_neutrons_unpol/analyzed_";
	//OSCM.simFile = getEnvSafe("G4OUTDIR")+"/20120824_MagF_neutronBetaUnpol/analyzed_";
	//std::string simFile="/home/mmendenhall/geant4/output/20120824_MagF_neutronBetaUnpol/analyzed_";
	std::string simFile="/home/mmendenhall/geant4/output/thinFoil_neutronBetaUnpol/analyzed_";
	//OSCM.simFile= getEnvSafe("G4OUTDIR")+"/endcap_180_150_neutronBetaUnpol/analyzed_";
	OSCM.simFactor = 1.0;
	OSCM.doPlots = true;
	
	/////////// Geant4 MagF
	OSCM.nTot = 104;
	OSCM.stride = 14;
	
	/////////// Geant4 0823, thinfoil
	//OSCM.nTot = 312;
	//OSCM.stride = 73;
	
	/////////// endcap_180_150
	//OSCM.nTot = 492;
	//OSCM.stride = 73;

	const std::string simOutName = "_Sim";
	const std::string simOutputDir=OSCM.outputDir+simOutName;
	
	if(octn < 0) {
		SimBetaDecayAnalyzer BDA_Sim(&OM,simOutputDir);
		BDA_Sim.simPerfectAsym = true;
		if(octn==-1000) {
			BetaDecayAnalyzer BDA(&OM,OSCM.outputDir,RunAccumulator::processedLocation);
			OSCM.combineSims(BDA_Sim,&BDA);
		} else if(octn==-1001) {
			BetaDecayAnalyzer BDA(&OM,OSCM.outputDir,RunAccumulator::processedLocation);
			SimBetaDecayAnalyzer BDA_MC(&OM,simOutputDir,OSCM.baseDir+"/"+simOutputDir+"/"+simOutputDir);
			BDA_MC.compareMCtoData(BDA);
		} else { OSCM.simOct(BDA_Sim,-octn-1); }
	} else {
		BetaDecayAnalyzer BDA(&OM,OSCM.outputDir);
		if(octn==1000) OSCM.combineOcts(BDA);
		else if(octn==1001) OSCM.recalcAllOctets(BDA,false);
		else { OSCM.scanOct(BDA, octn); }
	}
}

void mi_evis2etrue(StreamInteractor*) {
	OutputManager OM("Evis2ETrue",getEnvSafe("UCNA_ANA_PLOTS")+"/Evis2ETrue/20120810/");
	G4toPMT g2p;
	g2p.addFile(getEnvSafe(getEnvSafe("G4OUTDIR")+"/2011-2012geo_1mil_neutrons_unpol/analyzed_*.root"));
	PMTCalibrator PCal(16000);
	g2p.setCalibrator(PCal);
	SimSpectrumInfo(g2p,OM);
	OM.setWriteRoot(true);
	OM.write();
}

void mi_sourcelog(StreamInteractor*) { uploadRunSources(getEnvSafe("UCNA_RUN_LOG")); }

void mi_radcor(StreamInteractor* S) {
	float Ep = S->popFloat();
	int Z = S->popInt();
	int A = S->popInt();
	makeCorrectionsFile(A,Z,Ep);
}

void mi_showGenerator(StreamInteractor* S) {
	std::string sName = S->popString();
	OutputManager OMTest("test",getEnvSafe("UCNA_ANA_PLOTS")+"/test/EventGenerators/"+sName+"/");
	NucDecayLibrary NDL(getEnvSafe("UCNA_AUX")+"/NuclearDecays",1e-6);
	PMTCalibrator PCal(16000);
	showSimSpectrum(sName,OMTest,NDL,PCal);
	return;
}

void mi_showCal(StreamInteractor* S) {
	RunNum rn = S->popInt();
	PMTCalibrator PCal(rn);
}

void mi_makeSimSpectrum(StreamInteractor* S) {
	float eMax = S->popFloat();
	std::string simName = S->popString();
	
	RunNum rn = 16194;
	
	G4toPMT G2P;
	std::string fPath = getEnvSafe("G4OUTDIR")+"/"+simName;
	G2P.addFile(fPath+"/analyzed_*.root");

	PMTCalibrator PCal(rn);
	G2P.setCalibrator(PCal);
	
	OutputManager OM("SimSpectrum",getEnvSafe("UCNA_ANA_PLOTS")+"/SimSpectrum/");
	TH1F* hSpec = OM.registeredTH1F("hSpec","EventSpectrum",200,0,eMax);
	
	G2P.startScan(false);
	while(G2P.nextPoint()) {
		if(G2P.fType >= TYPE_IV_EVENT) continue;
		hSpec->Fill(G2P.getErecon());
	}
	
	double nOrigEvts = 3e6;
	
	hSpec->SetTitle(NULL);
	hSpec->GetXaxis()->SetTitle("Energy [keV]");
	hSpec->GetYaxis()->SetTitle("Events / keV / 1000 decays");
	hSpec->GetYaxis()->SetTitleOffset(1.2);
	hSpec->Scale(1.0e3/nOrigEvts/hSpec->GetBinWidth(1));
	hSpec->Draw();
	OM.printCanvas(simName);
}

void Analyzer(std::deque<std::string> args=std::deque<std::string>()) {
	
	ROOTStyleSetup();
	
	InputRequester exitMenu("Exit Menu",&menutils_Exit);
	InputRequester peek("Show stack",&menutils_PrintStack);
	
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
	InputRequester pm_posmap("Generate Position Map",&mi_EndpointStudy);
	pm_posmap.addArg("Start Run");
	pm_posmap.addArg("End Run");
	pm_posmap.addArg("n Rings","12");
	InputRequester pm_posmap_sim("Simulate Xe Position Map",&mi_EndpointStudySim);
	pm_posmap_sim.addArg("Start Run");
	pm_posmap_sim.addArg("End Run");
	pm_posmap_sim.addArg("n Rings","12");
	InputRequester pm_posmap_resim("Compare data/sim to create map",&mi_EndpointStudyReSim);
	pm_posmap_resim.addArg("Start Run");
	pm_posmap_resim.addArg("End Run");
	pm_posmap_resim.addArg("n Rings","12");
	InputRequester posmapLister("List Posmaps",&mi_listPosmaps);
	InputRequester posmapPlot("Plot Position Map",&mi_PosmapPlot);
	posmapPlot.addArg("Posmap ID");
	InputRequester posmapDumper("Dump Posmap",&mi_dumpPosmap);
	posmapDumper.addArg("Posmap ID");
	InputRequester posmapDel("Delete Posmap",&mi_delPosmap);
	posmapDel.addArg("Posmap ID");
	InputRequester nPEPlot("Plot nPE/MeV",&mi_nPEPlot);
	nPEPlot.addArg("Run Number");
	OptionsMenu PMapR("Position Map Routines");
	PMapR.addChoice(&pm_posmap,"gen");
	PMapR.addChoice(&pm_posmap_sim,"sim");
	PMapR.addChoice(&pm_posmap_resim,"comp");
	PMapR.addChoice(&posmapLister,"ls");
	PMapR.addChoice(&posmapPlot,"plot");
	PMapR.addChoice(&posmapDumper,"dump");
	PMapR.addChoice(&posmapDel,"del");
	PMapR.addChoice(&nPEPlot,"npe");
	PMapR.addChoice(&exitMenu,"x");
	
	// postprocessing/plots routines
	InputRequester dumpCalInfo("Dump calibration info to file",&mi_DumpCalInfo);
	dumpCalInfo.addArg("Start Run");
	dumpCalInfo.addArg("End Run");
	dumpCalInfo.addArg(&selectRuntype);
	
	InputRequester showCal("Show run calibration",&mi_showCal);
	showCal.addArg("Run");
	
	InputRequester showOcts("Show octet list",&mi_displayOctetList);
	
	InputRequester octetProcessor("Process Octet",&mi_processOctet);
	octetProcessor.addArg("Octet number");

	InputRequester showGenerator("Event generator test",&mi_showGenerator);
	showGenerator.addArg("Generator name");
	
	InputRequester makeSimSpectrum("Sim Spectrum",&mi_makeSimSpectrum);
	makeSimSpectrum.addArg("Sim name");
	makeSimSpectrum.addArg("energy range");
	
	// Posprocessing menu
	OptionsMenu PostRoutines("Postprocessing Routines");
	PostRoutines.addChoice(&showCal,"cal");
	PostRoutines.addChoice(&dumpCalInfo,"dcl");
	PostRoutines.addChoice(&showOcts,"sho");
	PostRoutines.addChoice(&octetProcessor,"oct");
	PostRoutines.addChoice(&showGenerator,"evg");
	PostRoutines.addChoice(&makeSimSpectrum,"mks");
	PostRoutines.addChoice(&exitMenu,"x");
	
	// sources
	InputRequester postSources("Fit source data",&mi_PostprocessSources);
	postSources.addArg("Start Run");
	postSources.addArg("End Run");
	InputRequester uploadSources("Upload runlog sources",&mi_sourcelog);
	// evis2etrue
	InputRequester evis2etrue("Caluculate eVis->eTrue curves",&mi_evis2etrue);
	// radiative corrections
	InputRequester radcor("Make radiative corrections table",&mi_radcor);
	radcor.addArg("A","1");
	radcor.addArg("Z","1");
	radcor.addArg("Endpoint",dtos(neutronBetaEp));
	
	// main menu
	OptionsMenu OM("Analyzer Main Menu");
	OM.addChoice(&PMapR,"pmap");
	OM.addChoice(&PostRoutines,"pr");
	OM.addChoice(&postSources,"sr");
	OM.addChoice(&uploadSources,"us");
	OM.addChoice(&evis2etrue,"ev");
	OM.addChoice(&radcor,"rc");
	OM.addChoice(&exitMenu,"x");
	OM.addSynonym("x","exit");
	OM.addSynonym("x","quit");
	OM.addSynonym("x","bye");
	OM.addChoice(&peek,"peek",SELECTOR_HIDDEN);
	
	std::stack<std::string> stack;
	OM.mydeque = &args;
	OM.mystack = &stack;
	OM.doIt();
	
	printf("\n\n\n>>>>> Goodbye. <<<<<\n\n\n");
}

int main(int argc, char *argv[]) {
	std::deque<std::string> args;
	for(int i=1; i<argc; i++)
		args.push_back(argv[i]);
	Analyzer(args);
	return 0;
}
