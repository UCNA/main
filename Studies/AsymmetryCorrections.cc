#include "AsymmetryCorrections.hh"
#include "AsymmetryAnalyzer.hh"
#include "SimAsymmetryAnalyzer.hh"
#include "GraphUtils.hh"
#include "PathUtils.hh"
#include "BetaSpectrum.hh"
#include <TRandom.h>
#include <stdio.h>
#include <cassert>
#include "SMExcept.hh"

TRandom3 acRndSrc;

//-------------------------------------------------------------//

AsymCorr::AsymCorr(const std::string& nm): name(nm) {
	std::string fname = getEnvSafe("UCNA_AUX")+"/Corrections/"+nm+".txt";
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
			int foo = fscanf(f,"%*c");
		}
	}
	fclose(f);
	
	printf("Loaded %i points for correction '%s'\n",gCor.GetN(),name.c_str());
}

void doFullCorrections(AsymmetryAnalyzer& AA, OutputManager& OM) {
	
	AA.calculateResults();
	AA.myA->defaultCanvas->cd();
	
	double emin = 275;
	double emax = 625;
	TF1 lineFit("lineFit","pol0",emin,emax);
	
	std::vector<AsymCorr*> ctable;
	ctable.push_back(new AsymCorr("Delta_2_0_"+ctos(choiceLetter(AA.anChoice))));
	ctable.push_back(new AsymCorr("Delta_2_1_"+ctos(choiceLetter(AA.anChoice))));
	ctable.push_back(new AsymCorr("Delta_2_2_"+ctos(choiceLetter(AA.anChoice))));
	ctable.push_back(new AsymCorr("Delta_2_3_"+ctos(choiceLetter(AA.anChoice))));
	ctable.push_back(new AsymCorr("Delta_3_"+ctos(choiceLetter(AA.anChoice))));
	ctable.push_back(new AsymCorr("MuonEffic"));
	ctable.push_back(new AsymCorr("PedShifts"));
	ctable.push_back(new AsymCorr("GainFlucts"));
	ctable.push_back(new AsymCorr("EnergyLinearityUncertainty_2010"));
	ctable.push_back(new AsymCorr("NGBG"));
	ctable.push_back(new AsymCorr("MWPCeffic"));
	ctable.push_back(new AsymCorr("MagF"));
	ctable.push_back(new AsymCorr("Radiative_h-g"));
	ctable.push_back(new AsymCorr("RecoilOrder"));
	
	// initial beta/2 correction
	for(int b=1; b<=AA.hAsym->GetNbinsX(); b++) {
		double e = AA.hAsym->GetBinCenter(b);
		double c = 2./beta(e);
		AA.hAsym->SetBinContent(b, AA.hAsym->GetBinContent(b)*c);
		AA.hAsym->SetBinError(b, AA.hAsym->GetBinError(b)*c);
	}
	AA.hAsym->Fit(&lineFit,"R");
	Stringmap m1;
	m1.insert("anChoice",ctos(choiceLetter(AA.anChoice)));
	m1.insert("A0",lineFit.GetParameter(0));
	m1.insert("dA0",lineFit.GetParError(0));
	m1.insert("chi2",lineFit.GetChisquare());
	m1.insert("ndf",lineFit.GetNDF());
	OM.qOut.insert("rawFit",m1);
	
	AA.hAsym->GetYaxis()->SetRangeUser(-0.2,0);
	AA.hAsym->Draw();
	AA.myA->printCanvas("hAsym_A1");
	
	// apply each correction
	for(unsigned int i=0; i<ctable.size(); i++) {
		for(int b=1; b<=AA.hAsym->GetNbinsX(); b++) {
			double e = AA.hAsym->GetBinCenter(b);
			double c = 1+ctable[i]->gCor.Eval(e);
			if(!(1e-1 < c && c <1e1)) continue;
			AA.hAsym->SetBinContent(b, AA.hAsym->GetBinContent(b)*c);
			AA.hAsym->SetBinError(b, AA.hAsym->GetBinError(b)*fabs(c));
		}
	}
	
	// integrate errors
	int b0 = AA.hAsym->FindBin(emin);
	int b1 = AA.hAsym->FindBin(emax);
	
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
	AA.hAsym->Draw();
	AA.myA->printCanvas("hAsym_A3");
	
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
			sumc[i] += ctable[i]->gCor.Eval(e)/(err*err);
			sumcerr[i] += ctable[i]->gUnc.Eval(e)/(err*err);
		}
	}
	mu /= statw;
	for(unsigned int i=0; i<ctable.size(); i++) {
		sumc[i]/=statw;
		sumcerr[i]/=statw;
	}
	printf("\n\nSTATISTICS: %.6f +- %.6f\n",mu,1./sqrt(statw));
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
	m3.insert("sysTot",systot);
	printf("SYSTEMATICS TOTAL = +- %.3f %%\n\n",100*systot);
	OM.qOut.insert("systematics",m3);
}



//-------------------------------------------------------------//

void calcMCCorrs(OutputManager& OM, const std::string& datin, const std::string& simin, bool writeAux, bool oldCorr) {
	for(AnalysisChoice a = ANCHOICE_A; a <= ANCHOICE_E; ++a) {
		OctetAnalyzer OAdat(&OM, "DataCorrector", datin);
		AsymmetryAnalyzer* AAdat = new AsymmetryAnalyzer(&OAdat);
		AAdat->anChoice = a;
		OAdat.addPlugin(AAdat);
		
		OctetAnalyzer OAsim(&OM, "Corr_Anchoice_"+ctos(choiceLetter(a)), simin);
		AsymmetryAnalyzer* AAsim = new AsymmetryAnalyzer(&OAsim);
		AAsim->anChoice = a;
		OAsim.addPlugin(AAsim);
		SimAsymmetryAnalyzer* SAAsim = new SimAsymmetryAnalyzer(&OAsim);
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
			if(writeAux) {
				FILE* f = fopen((getEnvSafe("UCNA_AUX")+"/Corrections/"+sname+".txt").c_str(),"w");
				if(i<=TYPE_III_EVENT+1)
					fprintf(f,"# MC Backscattering Correction Delta_{2,%i} for Analysis Choice %c\n",i-1,choiceLetter(a));
				else
					fprintf(f,"# MC Acceptance Correction Delta_3 for Analysis Choice %c\n",choiceLetter(a));
				fprintf(f,"#\n#E_lo\tE_hi\tcorrection\tuncertainty\n");
				for(unsigned int b=0; b<800; b+=10) {
					double e = b+5.;
					double dA = g->Eval(e);
					fprintf(f,"%i\t%i\t%g\t%g\n",b,b+10,dA,fabs(dA*0.30));
				}
				fclose(f);
			}
		}
	}
}

//-------------------------------------------------------------//

void compareMCs(OutputManager& OM, const std::string& sim0, const std::string& sim1, const std::string& fOut) {
	OctetAnalyzer OA0(&OM, "DataCorrector", sim0);
	AsymmetryAnalyzer* AA0 = new AsymmetryAnalyzer(&OA0);
	OA0.addPlugin(AA0);
	OctetAnalyzer OA1(&OM, "DataCorrector", sim1);
	AsymmetryAnalyzer* AA1 = new AsymmetryAnalyzer(&OA1);
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
				afpScale[s][2] = 0;
				for(AFPState afp = AFP_OFF; afp <= AFP_ON; ++afp)
					afpScale[s][2] += (afpScale[s][afp] = rAFP[afp]*acRndSrc.Gaus(1.0,dAFPfrac));
				for(AFPState afp = AFP_OFF; afp <= AFP_ON; ++afp)
					afpScale[s][afp] /= afpScale[s][2];
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
