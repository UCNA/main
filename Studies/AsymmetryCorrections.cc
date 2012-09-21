#include "AsymmetryCorrections.hh"
#include "AsymmetryAnalyzer.hh"
#include "SimAsymmetryAnalyzer.hh"
#include "GraphUtils.hh"
#include "PathUtils.hh"
#include "BetaSpectrum.hh"
#include <TRandom.h>
#include <stdio.h>
#include <cassert>

TRandom3 acRndSrc;

//-------------------------------------------------------------//

AsymCorr::AsymCorr(const std::string& nm): name(nm) {
	FILE* f = fopen((getEnvSafe("UCNA_AUX")+"/Corrections/"+nm+".txt").c_str(),"r");
	assert(f);
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
			fscanf(f,"%*c");
		}
	}
	fclose(f);
	
	printf("Loaded %i points for correction '%s'\n",gCor.GetN(),name.c_str());
}

void doFullCorrections(AsymmetryAnalyzer& AA) {
	
	AA.calculateResults();
	
	std::vector<AsymCorr*> ctable;
	ctable.push_back(new AsymCorr("Delta_2_0"));
	ctable.push_back(new AsymCorr("Delta_2_1"));
	ctable.push_back(new AsymCorr("Delta_2_2"));
	ctable.push_back(new AsymCorr("Delta_2_3"));
	ctable.push_back(new AsymCorr("MuonEffic"));
	ctable.push_back(new AsymCorr("PedShifts"));
	ctable.push_back(new AsymCorr("GainFlucts"));
	ctable.push_back(new AsymCorr("EnergyLinearityUncertainty_2010"));
	ctable.push_back(new AsymCorr("NGBG"));
	ctable.push_back(new AsymCorr("Radiative_h-g"));
	ctable.push_back(new AsymCorr("RecoilOrder"));
	
	// initial beta/2 correction
	for(int b=1; b<=AA.hAsym->GetNbinsX(); b++) {
		double e = AA.hAsym->GetBinCenter(b);
		double c = 2./beta(e);
		AA.hAsym->SetBinContent(b, AA.hAsym->GetBinContent(b)*c);
		AA.hAsym->SetBinError(b, AA.hAsym->GetBinError(b)*c);
	}
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
	double emin = 275;
	double emax = 625;
	int b0 = AA.hAsym->FindBin(emin);
	int b1 = AA.hAsym->FindBin(emax);

	// root average fit
	TF1 lineFit("lineFit","pol0",emin,emax);
	AA.hAsym->Fit(&lineFit,"R");
	AA.hAsym->GetYaxis()->SetRangeUser(-0.2,0);
	AA.hAsym->Draw();
	AA.myA->printCanvas("hAsym_A3");
		
	double statw = 0;
	double mu = 0;
	std::vector<double> sumc(ctable.size());
	std::vector<double> sumcerr(ctable.size());
	for(unsigned int b=b0; b<=b1; b++) {
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
	for(unsigned int i=0; i<ctable.size(); i++)
		printf("%s:\t%.5f +- %.5f %%\n",ctable[i]->name.c_str(),sumc[i],sumcerr[i]);
	
	// total systematics
	double systot = 0;
	for(unsigned int i=TYPE_0_EVENT; i<=TYPE_III_EVENT; i++)
		systot += sumcerr[i];
	systot *= systot;
	for(int i=TYPE_III_EVENT+1; i<ctable.size(); i++)
		systot += sumcerr[i]*sumcerr[i];
	systot = sqrt(systot);
	printf("SYSTEMATICS TOTAL = +- %.5f %%\n\n",systot);
}



//-------------------------------------------------------------//

void calcMCCorrs(OutputManager& OM, const std::string& datin, const std::string& simin) {
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
		
		std::vector<TH1*> asymStages = SAAsim->calculateCorrections(*AAdat,*AAsim);
		
		// re-bin for less awful stats
		for(unsigned int i=0; i<asymStages.size(); i++)
			asymStages[i]->Rebin(5);
		// back out each incremental correction stage
		for(unsigned int i = asymStages.size()-1; i>0; i--)
			asymStages[i]->Divide(asymStages[i-1]);
		
		// convert to graphs and draw
		for(unsigned int i=1; i<asymStages.size(); i++) {
			TGraph* g = new TGraph(asymStages[i]->GetNbinsX());
			for(int b=1; b<=asymStages[i]->GetNbinsX(); b++) {
				double c = asymStages[i]->GetBinContent(b);
				c = c>1.5?1.5:c<0.5?0.5:c==c?c:1;
				g->SetPoint(b-1,asymStages[i]->GetBinCenter(b),c-1);
			}
			g->SetTitle(("#Delta_{2,"+itos(i-1)+"} Correction").c_str());
			g->Draw("AP");
			g->GetYaxis()->SetRangeUser(-0.05,0.05);
			g->GetXaxis()->SetRangeUser(0,800);
			g->SetMarkerStyle(24);
			g->SetMarkerSize(0.5);
			g->SetMarkerColor(2);
			g->Draw("ACP");
			OAsim.printCanvas("Delta_2_"+itos(i-1));
			
			// write correction file
			if(a==ANCHOICE_C) {
				FILE* f = fopen((getEnvSafe("UCNA_AUX")+"/Corrections/Delta_2_"+itos(i-1)+".txt").c_str(),"w");
				fprintf(f,"# MC Backscattering Correction Delta_{2,%i}\n",i-1);
				fprintf(f,"#\n#E_lo\tE_hi\tcorrection\tuncertainty\n");
				for(unsigned int b=0; b<800; b+=10) {
					double e = b+5.;
					double dA = g->Eval(e);
					fprintf(f,"%i\t%i\t%g\t%g\n",b,b+10,dA,fabs(dA/3.));
				}
				fclose(f);
			}
			
		}
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
