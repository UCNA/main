#include <cassert>
#include "StyleSetup.hh"
#include "strutils.hh"
#include "GraphicsUtils.hh"
#include "OutputManager.hh"
#include "PlotMakers.hh"
#include "AsymmetryCorrections.hh"
#include "EnumerationFitter.hh"
#include "CathodeTuningAnalyzer.hh"
#include "SimEdepPlugin.hh"

int main(int argc, char *argv[]) {
	
	if(argc < 2) {
		printf("Have a nice day.\n");
		return 0;
	}
	const std::string rname = argv[1];
	
	printf("Running routine '%s'\n",rname.c_str());
	
	ROOTStyleSetup();
	
	if(rname=="anchoices_penelope") {
		//OutputManager OMdat("test",getEnvSafe("UCNA_ANA_PLOTS")+"/test/Anchoices/");
		//calcAnalysisChoices(OMdat, getEnvSafe("UCNA_ANA_PLOTS")+"/OctetAsym_Offic/OctetAsym_Offic");
		//return;
		std::string sim = "SimPen";
		OutputManager OMsim("test",getEnvSafe("UCNA_ANA_PLOTS")+"/test/Anchoices_"+sim+"/");
		calcAnalysisChoices(OMsim, getEnvSafe("UCNA_ANA_PLOTS")+"/OctetAsym_Offic_"+sim+"/OctetAsym_Offic_"+sim);
	}
	
	if(rname=="magf_mc_comparison") {
		std::string sim = "SimMagF_4x";
		OutputManager OM("test",getEnvSafe("UCNA_ANA_PLOTS")+"/test/MCChanges_"+sim+"/");
		compareMCs(OM, getEnvSafe("UCNA_ANA_PLOTS")+"/OctetAsym_Offic_Sim0823_4x/OctetAsym_Offic_Sim0823_4x",
				 getEnvSafe("UCNA_ANA_PLOTS")+"/OctetAsym_Offic_"+sim+"/OctetAsym_Offic_"+sim,"MagF");
	}
	
	std::string sim = "Sim0823_4x";
	//std::string sim = "SimPen";
	//std::string sim = "thinfoil";
	
	if(rname=="make_mc_correction") {
		// MC correction for data set
		OutputManager OM("test",getEnvSafe("UCNA_ANA_PLOTS")+"/test/MCCors_"+sim+"/");
		calcMCCorrs(OM, getEnvSafe("UCNA_ANA_PLOTS")+"/OctetAsym_Offic/OctetAsym_Offic",
				  getEnvSafe("UCNA_ANA_PLOTS")+"/OctetAsym_Offic_"+sim+"/OctetAsym_Offic_"+sim,
				  getEnvSafe("UCNA_AUX")+"/Corrections/", false);
	}
	
	if(rname=="make_mc_self_correction") {
		// MC self-correction
		sim = "thinfoil";
		std::string simOutNm = getEnvSafe("UCNA_ANA_PLOTS")+"/OctetAsym_Offic_"+sim+"/OctetAsym_Offic_"+sim;
		OutputManager OM("test",getEnvSafe("UCNA_ANA_PLOTS")+"/test/MCCors/"+sim+"/");
		calcMCCorrs(OM, simOutNm, simOutNm,getEnvSafe("UCNA_ANA_PLOTS")+"/test/MCCors/"+sim+"/");
		
		AnalysisChoice a = ANCHOICE_C;
		OctetAnalyzer OAdat(&OM, "DataCorrector_"+ctos(choiceLetter(a)), simOutNm);
		AsymmetryPlugin* AAdat = new AsymmetryPlugin(&OAdat);
		OAdat.addPlugin(AAdat);
		AAdat->anChoice = a;
		doFullCorrections(*AAdat,OM,OM.basePath);
		
		OM.write();
	}
	
	if(rname=="apply_mc_correction") {
		OutputManager OM("CorrectedAsym",getEnvSafe("UCNA_ANA_PLOTS")+"/test/CorrectAsym_"+sim+"/");
		for(AnalysisChoice a = ANCHOICE_A; a <= ANCHOICE_D; ++a) {
			if(a!=ANCHOICE_C) continue;
			OctetAnalyzer OAdat(&OM, "DataCorrector_"+ctos(choiceLetter(a)), getEnvSafe("UCNA_ANA_PLOTS")+"/OctetAsym_Offic/OctetAsym_Offic");
			AsymmetryPlugin* AAdat = new AsymmetryPlugin(&OAdat);
			OAdat.addPlugin(AAdat);
			AAdat->anChoice = a;
			doFullCorrections(*AAdat,OM);
		}
		OM.write();
	}
	
	if(rname=="sep_23") {
		separate23(sim);
	}
	
	if(rname=="sim_edep_scan") {
		G4toPMT G2P(true);
		G2P.addFile(getEnvSafe("G4OUTDIR")+"/20120824_MagF_neutronBetaUnpol/analyzed_*.root");
		PMTCalibrator PCal(16000);
		G2P.setCalibrator(PCal);
		
		OutputManager OM("SimSpectrum",getEnvSafe("UCNA_ANA_PLOTS")+"/SimSpectrum/SimEdep/");
		SimEdepAnalyzer SEA(&OM);
		
		G2P.startScan(false);
		while(G2P.nextPoint()) {
			SEA.fillCoreHists(G2P, G2P.physicsWeight);
		}
		
		SEA.calculateResults();
		SEA.myEdep->makeBigTable();
		SEA.write();
		SEA.setWriteRoot(true);
	}
	
	if(rname=="sim_edep") {
		OutputManager OM("SimSpectrum",getEnvSafe("UCNA_ANA_PLOTS")+"/SimSpectrum/SimEdep/");
		SimEdepAnalyzer SEA(&OM,"SimEdepAnalyzer",OM.basePath+"/SimEdepAnalyzer/SimEdepAnalyzer");
		SEA.calculateResults();
		SEA.makePlots();
		SEA.myEdep->makeBigTable();
	}

	
	if(rname=="fit_bg_excess") {
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
	}
	
	if(rname=="make_errtables") {
		ErrTables ET;
		//ET.gainfluctsTable(0.000125);
		//ET.pedShiftsTable(0.015);
		ET.eLinearityTable(2010);
		//ET.muonVetoEfficTable(0.002);
		//ET.NGBGTable(0.99,0.08,0.48,0.09, 0.25);
	}
	
	if(rname=="combo_ngbg_spectrum") {
		OutputManager OM("NGBG",getEnvSafe("UCNA_ANA_PLOTS")+"/NGBG/");
		NGBGAnalyzer AH(&OM,"Combined");
		
		NGBGAnalyzer AH1(&OM,"ScintFace_nCaptH",OM.basePath+"/ScintFace_nCaptH/ScintFace_nCaptH");
		AH1.scaleData(0.126);
		NGBGAnalyzer AH2(&OM,"DetAl_nCaptAl",OM.basePath+"/DetAl_nCaptAl/DetAl_nCaptAl");
		AH2.scaleData(0.073);
		NGBGAnalyzer AH3(&OM,"DetAl_nCaptAlGamma",OM.basePath+"/DetAl_nCaptAlGamma/DetAl_nCaptAlGamma");
		AH3.scaleData(0.594);
		
		AH.addSegment(AH1);
		AH.addSegment(AH2);
		AH.addSegment(AH3);
		
		NGBGAnalyzer AHdat(&OM,"NameUnused",getEnvSafe("UCNA_ANA_PLOTS")+"/OctetAsym_Offic/OctetAsym_Offic");
		AH.loadTotalTime(AHdat);
		AHdat.calculateResults();
		
		AH.calculateResults();
		AH.makePlots();
		
		AH.compareMCtoData(AHdat);
		
		TH1* h = AH.flipperSummedRate(AH.myAsym->qTotalSpectrum[EAST],GV_OPEN);
		h->Scale(1000000);
		h->SetTitle("MC neutron generated spectrum contribution");
		h->GetYaxis()->SetTitleOffset(1.25);
		h->GetYaxis()->SetTitle("event rate [uHz/keV]");
		h->Draw();
		AH.printCanvas("SuperSum");
		
		AH.write();
		AH.setWriteRoot(true);
	}
	
	if(rname=="ngbg_spectra") {
		//NGBGSpectra("EndcapEdge_nCaptH");
		//NGBGSpectra("EndcapEdge_nCaptCu");
		//NGBGSpectra("TrapWall_Cu66");
		//NGBGSpectra("TrapWall_nCaptCu");
		//NGBGSpectra("EntryPort_Al28");
		//NGBGSpectra("EntryPort_nCaptAl");
		NGBGSpectra("DetAl_nCaptAl");
		NGBGSpectra("DetAl_nCaptAlGamma");
		NGBGSpectra("ScintFace_nCaptH");
	}
	
	if(rname=="decompose_xenon") {
		//decomposeXenon(15991,true);
		//decomposeXenon(14282,false);
		//decomposeXenon(14347,false);
		//decomposeXenon(17224,false);
	}
	
	if(rname=="compare_xenon") {
		compareXenonSpectra();
	}
	
	if(rname=="interpl_err") {
		OutputManager OM("interpl_err",getEnvSafe("UCNA_ANA_PLOTS")+"/PositionMaps/interpl_err/");
		PosPlotter PP(&OM);
		PositioningCorrector* PCor = CalDBSQL::getCDB()->getPositioningCorrector(16000);
		
		PositioningCorrector P1;
		P1.interpType = CubiTerpolator::newCubiTerpolator;
		//P1.interpType =  Interpolator::newInterpolator;
		P1.loadData(PCor->getData());
		PositioningCorrector P2;
		P2.interpType = LinTerpolator::newLinTerpolator;
		P2.loadData(PCor->getData());
		
		PP.diffPlot(P1,P2,1.);
		
		PMTCalibrator PC1(16000);
		CalDBSQL::getCDB()->forgetPositioningCorrector(16000);
		PositioningCorrector::defaultInterpType = LinTerpolator::newLinTerpolator;
		PMTCalibrator PC2(16000);
		
		PP.npeDiffPlot(PC1, PC2);
	}
	
	if(rname=="sectors") {
		OutputManager OM("ThisNameNotUsedAnywhere",getEnvSafe("UCNA_ANA_PLOTS")+"/test/");
		SectorCutter SC(11,50);
		TH2F foo("foo","",10,-60,60,10,-60,60);
		foo.GetXaxis()->SetTitle("x position [mm]");
		foo.GetYaxis()->SetTitle("y position [mm]");
		foo.Draw();
		drawSectors(SC,1);
		labelSectors(SC,1);
		OM.printCanvas("Sectors");
	}
	
	if(rname=="gradposmap") {
		OutputManager OM("PosmapGradient",getEnvSafe("UCNA_ANA_PLOTS")+"/test/PosmapGradient/");
		PosPlotter PP(&OM);
		PMTCalibrator PCal(16000);
		PP.npeGradPlot(&PCal);
	}
	
	if(rname=="alpha_loss") {
		PMTCalibrator PCal(16000);
		QFile qOut;
		for(unsigned int i=0; i<20; i++) {
			G4toPMT G2P;
			G2P.setCalibrator(PCal);
			G2P.addFile(getEnvSafe("G4OUTDIR")+"/AlphaFoil_"+itos(i)+"_eGun_5485.6keV/analyzed_0.root");
			G2P.startScan(false);
			TH1F alphaLoss("alphaEnergy","alphaEnergy",6000,0,6000);
			while(G2P.nextPoint()) {
				alphaLoss.Fill(G2P.eDep[EAST]);
			}
			Stringmap m;
			m.insert("thickness",i);
			m.insert("mean",alphaLoss.GetMean());
			m.insert("rms",alphaLoss.GetRMS());
			qOut.insert("alpha_eloss",m);
			qOut.commit(getEnvSafe("UCNA_ANA_PLOTS")+"/test/AlphaEnergyLossSim.txt");
		}
		
	}
	
	return 0;
}
