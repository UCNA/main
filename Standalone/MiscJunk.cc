#include <cassert>
#include "StyleSetup.hh"
#include "strutils.hh"
#include "OutputManager.hh"
#include "PlotMakers.hh"
#include "AsymmetryCorrections.hh"
#include "EnumerationFitter.hh"
#include "CathodeTuningAnalyzer.hh"

int main(int argc, char *argv[]) {
	
	if(argc < 2) {
		printf("Have a nice day.\n");
		return 0;
	}
	std::string rname = argv[1];
	
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
		ET.gainfluctsTable(0.000125);
		ET.pedShiftsTable(0.015);
		ET.muonVetoEfficTable(0.002);
		ET.NGBGTable(0.99,0.08,0.48,0.09, 0.25);
	}
	
	if(rname=="combo_ngbg_spectrum") {
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
	
	return 0;
}
