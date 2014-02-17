#include "AsymmetryCorrections.hh"
#include "AsymmetryPlugin.hh"
#include "SimAsymmetryPlugin.hh"
#include "GraphUtils.hh"
#include "PathUtils.hh"
#include "BetaSpectrum.hh"
#include <TRandom.h>
#include <stdio.h>
#include <cassert>
#include <TLatex.h>
#include "SMExcept.hh"

TRandom3 acRndSrc;

//-------------------------------------------------------------//

AsymCorrFile::AsymCorrFile(const std::string& nm, const std::string& basepath): AsymCorr(nm) {
	std::string fname = basepath+"/"+nm+".txt";
	FILE* f = fopen(fname.c_str(),"r");
	if(!f) {
		SMExcept e("fileMissing");
		e.insert("fname",fname);
		throw(e);
	}
	float e0,e1,c,u;
	
	while(!feof(f)) {
		int n = fscanf(f,"%f %f %f %f",&e0,&e1,&c,&u);
		if(n==4) {
			unsigned int p = gCor.GetN();
			gCor.Set(p+1);
			gCor.SetPoint(p,0.5*(e0+e1),c);
			gUnc.Set(p+1);
			gUnc.SetPoint(p,0.5*(e0+e1),u);
		} else {
			n = fscanf(f,"%*c");
		}
	}
	fclose(f);
	
	printf("Loaded %i points for correction '%s'\n",gCor.GetN(),name.c_str());
}

void doFullCorrections(AsymmetryPlugin& AA, OutputManager& OM, std::string mcBase) {
	
	AA.calculateResults();
	AA.hCxn = (TH1F*)AA.hAsym->Clone("hCxn");
	AA.hCxn->Scale(0);
	AA.myA->defaultCanvas->cd();
	AA.myA->defaultCanvas->SetLeftMargin(0.11);
	AA.myA->defaultCanvas->SetRightMargin(0.04);
	
	// analysis energy window; 220-670 keV for 2010 data.
	double emin = 220;
	double emax = 670;
	
	TF1 lineFit("lineFit","pol0",emin,emax);
	TLatex lx;
	char tmp[1024];
	lx.SetTextColor(2);
	lx.SetTextSize(0.04);
	
	RunNum R0 = AA.myA->runCounts.counts.begin()->first;
	RunNum R1 = AA.myA->runCounts.counts.rbegin()->first;
	
	std::vector<AsymCorr*> ctable;
	if(!mcBase.size()) mcBase = getEnvSafe("UCNA_AUX")+"/Corrections/";
	ctable.push_back(new AsymCorrFile("Delta_2_0_"+ctos(choiceLetter(AA.anChoice)),mcBase));
	ctable.push_back(new AsymCorrFile("Delta_2_1_"+ctos(choiceLetter(AA.anChoice)),mcBase));
	ctable.push_back(new AsymCorrFile("Delta_2_2_"+ctos(choiceLetter(AA.anChoice)),mcBase));
	ctable.push_back(new AsymCorrFile("Delta_2_3_"+ctos(choiceLetter(AA.anChoice)),mcBase));
	ctable.push_back(new AsymCorrFile("Delta_3_"+ctos(choiceLetter(AA.anChoice)),mcBase));
	ctable.push_back(new AsymCorrFile("MuonEffic"));
	ctable.push_back(new AsymCorrFile("PedShifts"));
	ctable.push_back(new AsymCorrFile("GainFlucts"));
	ctable.push_back(new AsymCorrFile("EnergyLinearityUncertainty_2010"));
	ctable.push_back(new AsymCorrFile("NGBG"));
	ctable.push_back(new AsymCorrFile("MWPCeffic"));
	ctable.push_back(new AsymCorrFile("MagF"));
	if(R0==14077 && R1==14782)
		ctable.push_back(new ConstAsymCorr("Polarization",0.0044,0.00549));
	else if(R0==14888 && R1==16216)
		ctable.push_back(new ConstAsymCorr("Polarization",0.0099,0.00636));
	else
		ctable.push_back(new ConstAsymCorr("Polarization",0.0044*.59+.0099*.41,0.0057));
	ctable.push_back(new AsymCorrFile("Radiative_h-g"));
	ctable.push_back(new AsymCorrFile("RecoilOrder"));
	
	AA.hAsym->GetYaxis()->SetRangeUser(-0.1,0);
	AA.hAsym->SetTitle("Raw Asymmetry");
	AA.hAsym->Draw();
	AA.myA->printCanvas("RawAsym");
	TF1* asft = AA.hAsym->GetFunction("asymFit");
	if(asft) asft->SetBit(TF1::kNotDraw);
	
	// initial beta/2 correction
	for(int b=1; b<=AA.hAsym->GetNbinsX(); b++) {
		double e = AA.hAsym->GetBinCenter(b);
		double c = 2./beta(e);
		AA.hAsym->SetBinContent(b, AA.hAsym->GetBinContent(b)*c);
		AA.hAsym->SetBinError(b, AA.hAsym->GetBinError(b)*c);
	}
	AA.hAsym->GetYaxis()->SetRangeUser(-0.2,0);
	AA.hAsym->Fit(&lineFit,"NR");
	Stringmap m1;
	m1.insert("anChoice",ctos(choiceLetter(AA.anChoice)));
	m1.insert("A0",lineFit.GetParameter(0));
	m1.insert("dA0",lineFit.GetParError(0));
	m1.insert("chi2",lineFit.GetChisquare());
	m1.insert("ndf",lineFit.GetNDF());
	m1.insert("eMin",emin);
	m1.insert("eMax",emax);
	OM.qOut.insert("rawFit",m1);
	
	// save a copy of the uncorrected version
	TH1* hA_unc = (TH1*)AA.hAsym->Clone("hA_unc");
	hA_unc->SetTitle("Uncorrected Asymmetry");

	
	AA.hAsym->GetYaxis()->SetRangeUser(-0.2,0);
	AA.hAsym->SetTitle("Uncorrected Asymmetry");
	AA.hAsym->Draw();
	sprintf(tmp,"#splitline{A' = %.5f #pm %.5f,}{#chi^{2}/ndf = %.1f/%i}",lineFit.GetParameter(0),lineFit.GetParError(0),lineFit.GetChisquare(),lineFit.GetNDF());
	//lx.DrawLatex(100,-0.06,tmp);
	AA.myA->printCanvas("hAsym_A1");
	printf("\n*** Uncorrected %s ***\n",tmp);
	OM.addObject(AA.hAsym->Clone(("hAsym_Uncorrected_"+ctos(choiceLetter(AA.anChoice))).c_str()));
	
	// apply each correction
	
	for(int b=1; b<=AA.hAsym->GetNbinsX(); b++) {
		double e = AA.hAsym->GetBinCenter(b);
		double cTot = 0;
		double cErr = 0;
		Stringmap m;
		m.insert("KE",e);
		m.insert("anChoice",ctos(choiceLetter(AA.anChoice)));
		m.insert("Araw",AA.hAsym->GetBinContent(b));
		for(unsigned int i=0; i<ctable.size(); i++) {
			double c = 1+ctable[i]->getCor(e);
			if(!(0.02 < c && c < 50.)) continue;
			AA.hAsym->SetBinContent(b, AA.hAsym->GetBinContent(b)*c);
			AA.hAsym->SetBinError(b, AA.hAsym->GetBinError(b)*fabs(c));
			m.insert(ctable[i]->name,ctable[i]->getCor(e));
			m.insert("d_"+ctable[i]->name,ctable[i]->getUnc(e));
			if(i>=ctable.size()-2) continue; // exclude theory corrections from reported total
			if(i>=ctable.size()-3) continue; // exclude polarization from reported total
			cTot += ctable[i]->getCor(e);
			if(i <= TYPE_III_EVENT) cErr += ctable[i]->getUnc(e);
			if(i == TYPE_III_EVENT) {
				m.insert("Delta_2",cTot);
				m.insert("d_Delta_2",cErr);
				cErr *= cErr;
			}
			if(i > TYPE_III_EVENT) cErr += pow(ctable[i]->getUnc(e),2);
		}
		m.insert("cTot",cTot);
		m.insert("cErr",cErr);
		m.insert("A0",AA.hAsym->GetBinContent(b));
		OM.qOut.insert("corrPoint",m);
		AA.hCxn->SetBinContent(b,cTot);
		AA.hCxn->SetBinError(b,sqrt(cErr));
	}
	
	// integrate errors
	int b0 = AA.hAsym->FindBin(emin+0.5);
	int b1 = AA.hAsym->FindBin(emax-0.5);
	
	// root average fit
	AA.hAsym->Fit(&lineFit,"R");
	Stringmap m2;
	m2.insert("anChoice",ctos(choiceLetter(AA.anChoice)));
	m2.insert("A0",lineFit.GetParameter(0));
	m2.insert("dA0",lineFit.GetParError(0));
	m2.insert("chi2",lineFit.GetChisquare());
	m2.insert("ndf",lineFit.GetNDF());
	OM.qOut.insert("correctedFit",m2);
	
	AA.hAsym->GetYaxis()->SetRangeUser(-0.2,0);
	AA.hAsym->SetTitle("Asymmetry A_{0}");
	AA.hAsym->Draw();
	sprintf(tmp,"#splitline{A_{0} = %.5f #pm %.5f,}{#chi^{2}/ndf = %.1f/%i}",lineFit.GetParameter(0),lineFit.GetParError(0),lineFit.GetChisquare(),lineFit.GetNDF());
	//lx.DrawLatex(100,-0.06,tmp);
	AA.myA->printCanvas("hAsym_A3");
	printf("*** Corrected %s ***\n",tmp);
	OM.addObject(AA.hAsym->Clone(("hAsym_Corrected_"+ctos(choiceLetter(AA.anChoice))).c_str()));
	
	
	// corrected+uncorrected drawing
	AA.hAsym->GetYaxis()->SetRangeUser(-0.16,-0.08);
	AA.hAsym->GetYaxis()->SetTitleOffset(1.7);
	AA.hAsym->Draw("E0");
	hA_unc->SetMarkerStyle(20);
	hA_unc->SetMarkerSize(0.35);
	hA_unc->Draw("HIST SAME P");
	AA.myA->printCanvas("hAsym_Comparison");
	delete hA_unc;
	
	// manual bins averaging
	double statw = 0;
	double mu = 0;
	std::vector<double> sumc(ctable.size());
	std::vector<double> sumcerr(ctable.size());
	for(int b=b0; b<=b1; b++) {
		double err = AA.hAsym->GetBinError(b);
		mu += AA.hAsym->GetBinContent(b)/(err*err);
		statw += 1.0/(err*err);
		double e = AA.hAsym->GetBinCenter(b);
		for(unsigned int i=0; i<ctable.size(); i++) {
			sumc[i] += ctable[i]->getCor(e)/(err*err);
			sumcerr[i] += ctable[i]->getUnc(e)/(err*err);
		}
	}
	mu /= statw;
	for(unsigned int i=0; i<ctable.size(); i++) {
		sumc[i]/=statw;
		sumcerr[i]/=statw;
	}
	printf("\n\n---------------- %i - %i keV window ----------------\n",int(emin),int(emax));
	printf("\n%i-%i STATISTICS: %.6f +- %.6f\n",R0,R1,mu,1./sqrt(statw));
	Stringmap m3;
	m3.insert("anChoice",ctos(choiceLetter(AA.anChoice)));
	for(unsigned int i=0; i<ctable.size(); i++) {
		m3.insert(ctable[i]->name,sumc[i]);
		m3.insert("d_"+ctable[i]->name,sumcerr[i]);
		printf("%s:\t%.3f +- %.3f %%\n",ctable[i]->name.c_str(),100*sumc[i],100*sumcerr[i]);
	}
	
	// total systematics
	double systot = 0;
	for(unsigned int i=TYPE_0_EVENT; i<=TYPE_III_EVENT; i++)
		systot += sumcerr[i];
	systot *= systot;
	for(unsigned int i=TYPE_III_EVENT+1; i<ctable.size(); i++)
		systot += sumcerr[i]*sumcerr[i];
	systot = sqrt(systot);
	double corrtot = 0;
	for(unsigned int i=0; i<ctable.size(); i++)
		corrtot += sumc[i];
	
	m3.insert("sysTot",systot);
	m3.insert("corrTot",corrtot);
	printf("\nCORRECTIONS TOTAL = %.3f +- %.3f %%\n\n",100*corrtot,100*systot);
	OM.qOut.insert("systematics",m3);
}



//-------------------------------------------------------------//

void calcMCCorrs(OutputManager& OM, const std::string& datin, const std::string& simin, const std::string& outDir, bool oldCorr) {
	for(AnalysisChoice a = ANCHOICE_A; a <= ANCHOICE_E; ++a) {
		OctetAnalyzer OAdat(&OM, "DataCorrector", datin);
		AsymmetryPlugin* AAdat = new AsymmetryPlugin(&OAdat);
		AAdat->anChoice = a;
		OAdat.addPlugin(AAdat);
		
		OctetAnalyzer OAsim(&OM, "Corr_Anchoice_"+ctos(choiceLetter(a)), simin);
		AsymmetryPlugin* AAsim = new AsymmetryPlugin(&OAsim);
		AAsim->anChoice = a;
		OAsim.addPlugin(AAsim);
		SimAsymmetryPlugin* SAAsim = new SimAsymmetryPlugin(&OAsim);
		OAsim.addPlugin(SAAsim);
		
		std::vector<TH1*> asymStages = oldCorr?SAAsim->calculateCorrectionsOld(*AAsim):SAAsim->calculateCorrections(*AAdat,*AAsim);
		
		// re-bin for less awful stats
		for(unsigned int i=0; i<asymStages.size(); i++)
			asymStages[i]->Rebin(5);
		// back out each incremental correction stage
		for(unsigned int i = asymStages.size()-1; i>0; i--)
			asymStages[i]->Divide(asymStages[i-1]);
		
		// convert to graphs and draw
		for(unsigned int i=1; i<asymStages.size(); i++) {
			std::string sname = (i<=TYPE_III_EVENT+1?"Delta_2_"+itos(i-1)+"_":"Delta_3_")+ctos(choiceLetter(a));
			
			TGraph* g = new TGraph(asymStages[i]->GetNbinsX());
			for(int b=1; b<=asymStages[i]->GetNbinsX(); b++) {
				double c = asymStages[i]->GetBinContent(b);
				c = c>5.0?5.0:c<0.2?0.2:c==c?c:1;
				g->SetPoint(b-1,asymStages[i]->GetBinCenter(b),c-1);
			}
			if(i<=TYPE_III_EVENT+1)
				g->SetTitle(("#Delta_{2,"+itos(i-1)+"} Correction").c_str());
			else
				g->SetTitle("#Delta_{3}");
			g->Draw("AP");
			g->GetYaxis()->SetRangeUser(-0.05,0.05);
			g->GetXaxis()->SetRangeUser(0,800);
			g->SetMarkerStyle(24);
			g->SetMarkerSize(0.5);
			g->SetMarkerColor(2);
			g->Draw("ACP");
			OAsim.printCanvas(sname);
			
			// write correction file
			if(outDir.size()) {
				FILE* f = fopen((outDir+"/"+sname+".txt").c_str(),"w");
				if(i<=TYPE_III_EVENT+1)
					fprintf(f,"# MC Backscattering Correction Delta_{2,%i} for Analysis Choice %c\n",i-1,choiceLetter(a));
				else
					fprintf(f,"# MC Acceptance Correction Delta_3 for Analysis Choice %c\n",choiceLetter(a));
				fprintf(f,"#\n#E_lo\tE_hi\tcorrection\tuncertainty\n");
				for(unsigned int b=0; b<800; b+=10) {
					double e = b+5.;
					double dA = g->Eval(e);
					fprintf(f,"%i\t%i\t%g\t%g\n",b,b+10,dA,fabs(dA*0.25));
				}
				fclose(f);
			}
		}
	}
}

//-------------------------------------------------------------//

void compareMCs(OutputManager& OM, const std::string& sim0, const std::string& sim1, const std::string& fOut) {
	OctetAnalyzer OA0(&OM, "DataCorrector", sim0);
	AsymmetryPlugin* AA0 = new AsymmetryPlugin(&OA0);
	OA0.addPlugin(AA0);
	OctetAnalyzer OA1(&OM, "DataCorrector", sim1);
	AsymmetryPlugin* AA1 = new AsymmetryPlugin(&OA1);
	OA1.addPlugin(AA1);
	
	AA1->anChoice = AA0->anChoice;
	
	OA0.calculateResults();
	OA1.calculateResults();
	
	AA0->hAsym->Rebin(10);
	AA1->hAsym->Rebin(10);
	
	AA0->hAsym->Divide(AA1->hAsym);
	
	TGraph* g = new TGraph(AA0->hAsym->GetNbinsX());
	for(int b=1; b<=AA0->hAsym->GetNbinsX(); b++) {
		double c = AA0->hAsym->GetBinContent(b);
		c = c>5.0?5.0:c<0.2?0.2:c==c?c:1;
		g->SetPoint(b-1,AA0->hAsym->GetBinCenter(b),c-1);
	}
	g->SetTitle("MC Correction");
	g->Draw("AL");
	g->GetYaxis()->SetRangeUser(-0.01,0.01);
	g->GetXaxis()->SetRangeUser(0,800);
	g->SetMarkerStyle(24);
	g->SetMarkerSize(0.5);
	g->SetMarkerColor(2);
	g->Draw("ALP");
	OM.printCanvas("AsymDiff_"+ctos(choiceLetter(AA0->anChoice)));
	
	// write correction file
	if(fOut.size()>0) {
		FILE* f = fopen((getEnvSafe("UCNA_AUX")+"/Corrections/"+fOut+".txt").c_str(),"w");
		fprintf(f,"# MC correction\n");
		fprintf(f,"#\n#E_lo\tE_hi\tcorrection\tuncertainty\n");
		for(unsigned int b=0; b<800; b+=10) {
			double e = b+5.;
			double dA = g->Eval(e);
			fprintf(f,"%i\t%i\t%g\t%g\n",b,b+10,dA,fabs(dA*0.5));
		}
		fclose(f);
	}
}

//-------------------------------------------------------------//


ErrTables::ErrTables(const std::string& datset):
OM("nameUnused",getEnvSafe("UCNA_ANA_PLOTS")),
Adat(&OM, "nameUnused", getEnvSafe("UCNA_ANA_PLOTS")+"/"+datset+"/"+datset) {
	Adat.calculateResults();
	for(Side s = EAST; s <= WEST; ++s)
		for(AFPState afp = AFP_OFF; afp <= AFP_ON; ++afp)
			S[s][afp] = TH1toTGraph(*Adat.myAsym->qTotalSpectrum[s]->fgbg[afp]->h[GV_OPEN]);
}

ErrTables::~ErrTables() {
	for(Side s = EAST; s <= WEST; ++s)
		for(AFPState afp = AFP_OFF; afp <= AFP_ON; ++afp)
			delete S[s][afp];
}

 double ErrTables::energyErrorEnvelope(double e, unsigned int year) const {
 	assert(year==2010);
	double err = e*0.0125;
	if(err<2.5) return 2.5;
	if(err>500*0.0125) return 500*0.0125;
	return err;
 }

double ErrTables::getRexp(double e) const {
	return (S[EAST][AFP_OFF]->Eval(e)*S[WEST][AFP_ON]->Eval(e) /
			(S[EAST][AFP_ON]->Eval(e)*S[WEST][AFP_OFF]->Eval(e)));
}

double ErrTables::AofR(double R) {
	double A =  (1-sqrt(R))/(1+sqrt(R));
	return (A==A)?A:0;
}

double ErrTables::getAexp(double e) const {
	double A = AofR(getRexp(e));
	return fabs(A)>.01?A:-.01;
}

void ErrTables::gainfluctsTable(double delta) {
	FILE* f = fopen((getEnvSafe("UCNA_AUX")+"/Corrections/GainFlucts.txt").c_str(),"w");
	fprintf(f,"# Run-to-run gain fluctuations uncertainty for anticorrelated delta = %g\n",delta);
	fprintf(f,"#\n#E_lo\tE_hi\tcorrection\tuncertainty\n");
	for(unsigned int b=0; b<800; b+=10) {
		double e = b+5.;
		double A = getAexp(e);
		double Rp = (S[EAST][AFP_OFF]->Eval(e)*S[WEST][AFP_ON]->Eval(e/(1-delta)) /
					 (S[EAST][AFP_ON]->Eval(e/(1+delta))*S[WEST][AFP_OFF]->Eval(e)));
		double Ap = AofR(Rp);
		fprintf(f,"%i\t%i\t%g\t%g\n",b,b+10,0.,(A-Ap)/A);
	}
	fclose(f);
}

void ErrTables::pedShiftsTable(double delta) {
	FILE* f = fopen((getEnvSafe("UCNA_AUX")+"/Corrections/PedShifts.txt").c_str(),"w");
	fprintf(f,"# Run-to-run pedestal shifts uncertainty for anticorrelated delta/sqrt(N) = %g keV\n",delta);
	fprintf(f,"#\n#E_lo\tE_hi\tcorrection\tuncertainty\n");
	for(unsigned int b=0; b<800; b+=10) {
		double e = b+5.;
		double A = getAexp(e);
		double Rp = (S[EAST][AFP_OFF]->Eval(e)*S[WEST][AFP_ON]->Eval(e+delta) /
					 (S[EAST][AFP_ON]->Eval(e-delta)*S[WEST][AFP_OFF]->Eval(e)));
		double Ap = AofR(Rp);
		fprintf(f,"%i\t%i\t%g\t%g\n",b,b+10,0.,(A-Ap)/A);
	}
	fclose(f);
}

void ErrTables::eLinearityTable(unsigned int yr) {
	// this doesn't work! too much statistical scatter in asymmetry data. Use smooth fit to asymmetry, instead.
	FILE* f = fopen((getEnvSafe("UCNA_AUX")+"/Corrections/EnergyLinearity_"+itos(yr)+".txt").c_str(),"w");
	fprintf(f,"# Energy reconstruction errors for %i error envelope, using observed spectra\n",yr);
	fprintf(f,"#\n#E_lo\tE_hi\tcorrection\tuncertainty\n");
	for(unsigned int b=0; b<800; b+=10) {
		double e = b+5.;
		double A = getAexp(e);
		double ee = e+energyErrorEnvelope(e,yr)*0.01;
		double Rp = (S[EAST][AFP_OFF]->Eval(ee)*S[WEST][AFP_ON]->Eval(ee) /
					 (S[EAST][AFP_ON]->Eval(ee)*S[WEST][AFP_OFF]->Eval(ee)));
		double Ap = AofR(Rp);
		fprintf(f,"%i\t%i\t%g\t%g\n",b,b+10,0.,100*(A-Ap)/A);
	}
	fclose(f);
}

void ErrTables::muonVetoEfficTable(double delta) {
	// load spectrum data
	TGraphErrors* M[2][2];
	for(Side s = EAST; s <= WEST; ++s)
		for(AFPState afp = AFP_OFF; afp <= AFP_ON; ++afp)
			M[s][afp] = TH1toTGraph(*Adat.myMuons->qMuonSpectra[s][false]->fgbg[afp]->h[GV_OPEN]);
	
	// write uncertainty file
	FILE* f = fopen((getEnvSafe("UCNA_AUX")+"/Corrections/MuonEffic.txt").c_str(),"w");
	fprintf(f,"# Run-to-run muon veto efficiency shifts by anticorrelated delta/sqrt(N) = %g%%\n",delta*100);
	fprintf(f,"#\n#E_lo\tE_hi\tcorrection\tuncertainty\n");
	for(unsigned int b=0; b<800; b+=10) {
		double e = b+5.;
		double A = getAexp(e);
		double Rp = ((S[EAST][AFP_OFF]->Eval(e)+0.5*delta*M[EAST][AFP_OFF]->Eval(e)) *
					 (S[WEST][AFP_ON]->Eval(e)+0.5*delta*M[WEST][AFP_ON]->Eval(e)) /
					 (S[EAST][AFP_ON]->Eval(e)-0.5*delta*M[EAST][AFP_ON]->Eval(e)) /
					 (S[WEST][AFP_OFF]->Eval(e)-0.5*delta*M[WEST][AFP_OFF]->Eval(e)));
		double Ap = AofR(Rp);
		fprintf(f,"%i\t%i\t%g\t%g\n",b,b+10,0.,(A-Ap)/A);
	}
	fclose(f);
	
	// cleanup
	for(Side s = EAST; s <= WEST; ++s)
		for(AFPState afp = AFP_OFF; afp <= AFP_ON; ++afp)
			delete M[s][afp];
}

void ErrTables::efficShiftTable(double delta) {
	FILE* f = fopen((getEnvSafe("UCNA_AUX")+"/Corrections/EfficShifts.txt").c_str(),"w");
	fprintf(f,"# Uncertainty for uniform anticorrelated efficiency shift of %g\n",delta);
	fprintf(f,"#\n#E_lo\tE_hi\tcorrection\tuncertainty\n");
	for(unsigned int b=0; b<800; b+=10) {
		double e = b+5.;
		double A = getAexp(e);
		double Rp = getRexp(e)*pow((1+0.5*delta)/(1-0.5*delta),2);
		double Ap = AofR(Rp);
		fprintf(f,"%i\t%i\t%g\t%g\n",b,b+10,0.,(A-Ap)/A);
	}
	fclose(f);
}

void ErrTables::NGBGTable(double EScale, double dEScale, double WScale, double dWScale, double dAFPfrac) {
	FILE* f = fopen((getEnvSafe("UCNA_AUX")+"/Corrections/NGBG.txt").c_str(),"w");
	fprintf(f,"# Uncertainty for neutron-generated backgrounds\n");
	fprintf(f,"#\n#E_lo\tE_hi\tcorrection\tuncertainty\n");
	
	// load NGBG estimate
	OutputManager OM2("NGBG",getEnvSafe("UCNA_ANA_PLOTS")+"/NGBG/");
	SimBetaDecayAnalyzer AH(&OM2,"Combined",OM2.basePath+"/Combined/Combined");
	AH.calculateResults();
	TGraph* hNGBG = TH1toTGraph(*AH.myAsym->qTotalSpectrum[EAST]->fgbg[AFP_OFF]->h[GV_OPEN]);
	
	// AFP rates
	double rAFP[2] = {25.8,16.7};
	
	double de = 10;
	for(unsigned int b=0; b<800; b+=de) {
		
		double e = b+0.5*de;
		double A = getAexp(e);
		double NGBG = hNGBG->Eval(e);
		
		// MC for errors
		double sx = 0;
		double sxx = 0;
		unsigned int n = 500;
		for(unsigned int i=0; i<n; i++) {
			// scaling for side
			double sideScale[2] = {EScale+acRndSrc.Gaus(0.,dEScale),WScale+acRndSrc.Gaus(0.,dWScale)};
			// calculate scaling between On/Off
			double afpScale[2][2];
			for(Side s = EAST; s <= WEST; ++s) {
				double sumScale = 0;
				for(AFPState afp = AFP_OFF; afp <= AFP_ON; ++afp) {
					afpScale[s][afp] = rAFP[afp]*acRndSrc.Gaus(1.0,dAFPfrac);
					sumScale += afpScale[s][afp];
				}
				for(AFPState afp = AFP_OFF; afp <= AFP_ON; ++afp)
					afpScale[s][afp] /= sumScale;
			}
			// asymmetry + background
			double Rp = ((S[EAST][AFP_OFF]->Eval(e)+NGBG*afpScale[EAST][AFP_OFF]*sideScale[EAST]) *
						 (S[WEST][AFP_ON]->Eval(e)+NGBG*afpScale[WEST][AFP_ON]*sideScale[WEST]) /
						 ((S[EAST][AFP_ON]->Eval(e)+NGBG*afpScale[EAST][AFP_ON]*sideScale[EAST]) *
						  (S[WEST][AFP_OFF]->Eval(e)+NGBG*afpScale[WEST][AFP_OFF]*sideScale[WEST])));
			double Ap = AofR(Rp);
			sx += (A-Ap)/A;
			sxx += (A-Ap)/A*(A-Ap)/A;
		}
		sx /= n;
		sxx /= n;
		
		fprintf(f,"%i\t%i\t%g\t%g\n",(int)b,(int)(b+de),sx,sqrt(sxx-sx*sx));
	}
	fclose(f);
}
