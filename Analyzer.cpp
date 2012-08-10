#include <exception>
#include <TStyle.h>
#include <TROOT.h>
#include <iostream>

#include "ControlMenu.hh"
#include "PathUtils.hh"
#include "GraphicsUtils.hh"
#include "PositionResponse.hh"
#include "BetaDecayAnalyzer.hh"
#include "EndpointStudy.hh"
#include "PostOfficialAnalyzer.hh"
#include "PlotMakers.hh"
#include "PMTGenerator.hh"
#include "ReSource.hh"
#include "G4toPMT.hh"
#include "LEDScans.hh"
#include "NuclEvtGen.hh"
#include "EnumerationFitter.hh"

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

void mi_EndpointStudy(std::deque<std::string>&, std::stack<std::string>& stack) {
	unsigned int nr = streamInteractor::popInt(stack);
	RunNum r1 = streamInteractor::popInt(stack);
	RunNum r0 = streamInteractor::popInt(stack);
	process_xenon(r0,r1,nr);
}

void mi_EndpointStudySim(std::deque<std::string>&, std::stack<std::string>& stack) {
	unsigned int nRings = streamInteractor::popInt(stack);
	RunNum rsingle = streamInteractor::popInt(stack);
	RunNum r1 = streamInteractor::popInt(stack);
	RunNum r0 = streamInteractor::popInt(stack);
	simulate_xenon(r0,r1,rsingle,nRings);
}

void mi_EndpointStudyReSim(std::deque<std::string>&, std::stack<std::string>& stack) {
	unsigned int nRings = streamInteractor::popInt(stack);
	RunNum r1 = streamInteractor::popInt(stack);
	RunNum r0 = streamInteractor::popInt(stack);
	xenon_posmap(r0,r1,nRings);
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
	PMTCalibrator PCal(rn);
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
	
	const std::string simOutputDir=outputDir+"_Sim_MagF_2";
	
	std::string simFile="/home/mmendenhall/geant4/output/WideKev_neutronBetaUnpol/analyzed_";
	
	unsigned int nTot = 52;
	unsigned int stride = 7;
	
	BetaDecayAnalyzer::processedLocation = getEnvSafe("UCNA_ANA_PLOTS")+"/"+outputDir+"/"+outputDir;
	
	if(octn==1000) {
		OutputManager OM("ThisNameIsNotUsedAnywhere",getEnvSafe("UCNA_ANA_PLOTS"));
		BetaDecayAnalyzer AA(&OM,outputDir);
		processOctets(AA,Octet::loadOctets(QFile(getEnvSafe("UCNA_OCTET_LIST"))),365*24*3600);
	} else if(octn < 0) {
		G4toPMT simData;
		for(unsigned int i=0; i<stride; i++)
			simData.addFile(simFile+itos((stride*abs(octn)+i)%nTot)+".root");
		simData.PGen[EAST].xscatter = simData.PGen[WEST].xscatter = 0.01;
		if(octn==-1000) {
			OutputManager OM("ThisNameIsNotUsedAnywhere",getEnvSafe("UCNA_ANA_PLOTS"));
			SimBetaDecayAnalyzer AA_Sim(&OM,simOutputDir);
			AA_Sim.simPerfectAsym = true;
			AA_Sim.simuClone(getEnvSafe("UCNA_ANA_PLOTS")+"/"+outputDir, simData, 1.0, 365*24*3600);
		} else {
			Octet oct = Octet::loadOctet(QFile(getEnvSafe("UCNA_OCTET_LIST")),-octn-1);
			if(!oct.getNRuns()) return;
			OutputManager OM("ThisNameIsNotUsedAnywhere",getEnvSafe("UCNA_ANA_PLOTS")+"/"+simOutputDir);
			SimBetaDecayAnalyzer AA_Sim(&OM,oct.octName());
			AA_Sim.simPerfectAsym = true;
			AA_Sim.simuClone(getEnvSafe("UCNA_ANA_PLOTS")+"/"+outputDir+"/"+oct.octName(), simData, 1.0, 0.*3600);
		}
	} else {
		Octet oct = Octet::loadOctet(QFile(getEnvSafe("UCNA_OCTET_LIST")),octn);
		if(!oct.getNRuns()) return;
		OutputManager OM("ThisNameIsNotUsedAnywhere",getEnvSafe("UCNA_ANA_PLOTS")+"/"+outputDir);
		BetaDecayAnalyzer AA(&OM,oct.octName());
		processOctets(AA,oct.getSubdivs(oct.divlevel+1,false),10.*24*3600);
	}
}

void mi_evis2etrue(std::deque<std::string>&, std::stack<std::string>&) {
	OutputManager OM("Evis2ETrue",getEnvSafe("UCNA_ANA_PLOTS")+"/Evis2ETrue/WideKev/");
	G4toPMT g2p;
	g2p.addFile("/home/mmendenhall/geant4/output/WideKev_neutronBetaUnpol/analyzed_*.root");
	PMTCalibrator PCal(16000);
	g2p.setCalibrator(PCal);
	SimSpectrumInfo(g2p,OM);
	OM.setWriteRoot(true);
	OM.write();
}

void mi_sourcelog(std::deque<std::string>&, std::stack<std::string>&) { uploadRunSources(); }

void mi_radcor(std::deque<std::string>&, std::stack<std::string>& stack) {
	float Ep = streamInteractor::popFloat(stack);
	int Z = streamInteractor::popInt(stack);
	int A = streamInteractor::popInt(stack);
	makeCorrectionsFile(A,Z,Ep);
}

void mi_misc(std::deque<std::string>&, std::stack<std::string>&) {
	
	if(false) {
		OutputManager OMTest("test",getEnvSafe("UCNA_ANA_PLOTS")+"/test");
		EnumerationFitter EF;
		TGraphErrors* g = EF.loadFitFile("/home/mmendenhall/BGExcess_Combo_W.txt");
		TF1* f = EF.getFitter();
		f->SetRange(0,g->GetN());
		for(unsigned int i=0; i<EF.getNParams(); i++)
			f->SetParLimits(i,0.,100.);
		g->Fit(f,"R");
		printf("Chi2/ndf = %g/%i\n",f->GetChisquare(),f->GetNDF());
		g->Draw("A*");
		OMTest.printCanvas("BG_Components_ComboW");
		return;
	}
		
	if(false) {
		OutputManager OMTest("test",getEnvSafe("UCNA_ANA_PLOTS")+"/test");
		NucDecayLibrary NDL(getEnvSafe("UCNA_AUX")+"/NuclearDecays",1e-6);
		PMTCalibrator PCal(16000);
		showSimSpectrum("Cu66",OMTest,NDL,PCal);
		return;
	}
	
	if(true) {
		ErrTables ET;
		ET.gainfluctsTable(0.000125);
		ET.pedShiftsTable(0.015);
		ET.muonVetoEfficTable(0.002);
		ET.NGBGTable(0.99,0.08,0.48,0.09, 0.25);
		return;
	}
	
	if(false) {
		OutputManager OM("NGBG",getEnvSafe("UCNA_ANA_PLOTS")+"/NGBG/");
		SimBetaDecayAnalyzer AH(&OM,"Combined");
		
		SimBetaDecayAnalyzer AH1(&OM,"ScintFace_nCaptH",OM.basePath+"/ScintFace_nCaptH/ScintFace_nCaptH");
		AH1.scaleData(0.126);
		SimBetaDecayAnalyzer AH2(&OM,"DetAl_nCaptAl",OM.basePath+"/DetAl_nCaptAl/DetAl_nCaptAl");
		AH2.scaleData(0.073);
		SimBetaDecayAnalyzer AH3(&OM,"DetAl_nCaptAlGamma",OM.basePath+"/DetAl_nCaptAlGamma/DetAl_nCaptAlGamma");
		AH3.scaleData(0.594);
		
		AH.addSegment(AH1);
		AH.addSegment(AH2);
		AH.addSegment(AH3);
		
		AH.calculateResults();
		AH.makePlots();
		AH.write();
		AH.setWriteRoot(true);
	}
	
	//NGBGSpectra("EndcapEdge_nCaptH");
	//NGBGSpectra("EndcapEdge_nCaptCu");
	//NGBGSpectra("TrapWall_Cu66");
	//NGBGSpectra("TrapWall_nCaptCu");
	//NGBGSpectra("EntryPort_Al28");
	//NGBGSpectra("EntryPort_nCaptAl");
	NGBGSpectra("DetAl_nCaptAl");
	NGBGSpectra("DetAl_nCaptAlGamma");
	NGBGSpectra("ScintFace_nCaptH");
	return;
	
	processWirechamberCal(14264,16077,20);
	return;
	
	//decomposeXenon(15991,true);
	//decomposeXenon(14282,false);
	//decomposeXenon(14347,false);
	//decomposeXenon(17224,false);
	//return;
	
	compareXenonSpectra();
	return;
	
	
	
	//spectrumGenTest();
		
	/*
	 OutputManager OMLS("PMTCorr",getEnvSafe("UCNA_ANA_PLOTS")+"/PMTCorr");
	 LEDScanScanner LSS;
	 std::vector<RunNum> bruns = selectRuns(16193, 16216, "beta");
	 //std::vector<RunNum> bruns = selectRuns(16097, 16216, "beta");
	 LSS.addRuns(bruns);
	 PMT_LED_Correlations(OMLS,LSS);
	 exit(0);	
	 */
}

void Analyzer(std::deque<std::string> args=std::deque<std::string>()) {
	
	gROOT->SetStyle("Plain");
	gStyle->SetPalette(1);
	gStyle->SetNumberContours(255);
	gStyle->SetOptStat("e");	
	TCanvas defaultCanvas;
	defaultCanvas.SetFillColor(0);
	defaultCanvas.SetCanvasSize(300,300);
	
	inputRequester exitMenu("Exit Menu",&menutils_Exit);
	inputRequester peek("Show stack",&menutils_PrintStack);
	inputRequester doMisc("Misc",&mi_misc);
	
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
	pm_posmap_sim.addArg("n Rings","12");
	inputRequester pm_posmap_resim("Reupload Position Map",&mi_EndpointStudyReSim);
	pm_posmap_resim.addArg("Start Run");
	pm_posmap_resim.addArg("End Run");
	pm_posmap_resim.addArg("n Rings","12");
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
	PMapR.addChoice(&pm_posmap_resim,"rup");
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
	OM.addChoice(&doMisc,"msc");
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
