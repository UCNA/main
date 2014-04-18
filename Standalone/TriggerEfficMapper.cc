#include "TChainScanner.hh"
#include "PointCloudHistogram.hh"
#include "PathUtils.hh"
#include "OutputManager.hh"
#include "Enums.hh"
#include "strutils.hh"
#include "GraphUtils.hh"
#include "GraphicsUtils.hh"

#include <TGraph.h>
#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Math/Functor.h>

#include <cassert>

/// Reader for backscatter trigger events tree
class TrigTreeScanner: public TChainScanner {
public:
	/// constructor
	TrigTreeScanner(): TChainScanner("trig") {}
	
	Int_t trigFlags;				///< channel trigger bit flags
	Float_t adc[BOTH][nBetaTubes];	///< pedestal-subtracted ADC for each PMT
	Int_t nf[BOTH];					///< number of firing PMTs on each side
	
protected:
	
	/// over-write this in subclass to automaticlly set readout points on first loaded file
	virtual void setReadpoints() {
		SetBranchAddress("trigFlags",&trigFlags);
		for(Side s = EAST; s <= WEST; ++s) {
			SetBranchAddress(sideSubst("ADC%c",s), &adc[s]);
			//SetBranchAddress("nFiring", nf);
		}
	}
};


/// Class for optimizing line cut separating two event classes
class SeparationLineCutOptimizer {
public:
	/// constructor
	SeparationLineCutOptimizer(const PointCloudHistogram* lo, const PointCloudHistogram* hi, unsigned int a):
	hLo(lo), hHi(hi), gSep(NULL), axs(a), minFit("minFit","pol3",0,1) {
		assert(hLo && hHi);
		assert(hLo->getNdim() == hHi->getNdim());
		vProj = new float[hLo->getNdim()];
	}
	/// destructor
	//virtual ~SeparationLineCutOptimizer() {
	//	delete vProj;
	//	if(gSep) delete gSep;
	//}
	
	/// determine minimum of cubic fit
	double fitMinPt() const {
		double p1 = 3*minFit.GetParameter(3);
		double p2 = 2*minFit.GetParameter(2);
		double p3 = 1*minFit.GetParameter(1);
		double del = sqrt(p2*p2-4*p1*p3);
		double x1 = (-p2+del)/(2*p1);
		double x2 = (-p2-del)/(2*p1);
		return minFit.Eval(x1)<minFit.Eval(x2)?x1:x2;
	}
	
	/// optimized evaluation
	double operator()(const double* dProj) {
	
		double s = 0;
		unsigned int n=0;
		for(unsigned int i=0; i<hHi->getNdim(); i++) {
			if(i==axs) continue;
			s += dProj[n]*dProj[n];
			vProj[i] = dProj[n++];
		}
		vProj[axs] = sqrt(1.-s);
		
		TGraph gLo, gHi;
		hLo->project(vProj, gLo);
		hHi->project(vProj, gHi);
		//if(gSep) delete gSep;	// TODO
		gSep = graphsep(gHi, gLo);
		
		int nmin = TMath::LocMin(gSep->GetN(), gSep->GetY());
		double c = gSep->GetX()[nmin];
		minFit.SetRange(0.5*c,1.5*c);
		gSep->Fit(&minFit,"QRN");
		xdiv = fitMinPt();
		
		return minFit.Eval(xdiv);
	}
	
	const PointCloudHistogram* hLo;	///< events preferably below cut
	const PointCloudHistogram* hHi;	///< events preferably above cut
	float* vProj;					///< space for float version of vector
	TGraph* gSep;					///< separation optimization plot
	double xdiv;					///< optimized separation point
	unsigned int axs;				///< main projection axis
	TF1 minFit;						///< fitter for minimum point
	
};

int main(int, char**) {
	
	// set up histograms
	printf("Setting up histograms...\n");
	KDTreeSet Tpts(nBetaTubes);
	const float xlo[] = {-50., -50., -50., -50.};
	const float xhi[] = {150., 150., 150., 150.};
	Tpts.fillPointRange(30*30*30*30, xlo, xhi);
	const float xlo2[] = {-30., -30., -30., -30.};
	const float xhi2[] = {30., 30., 30., 30.};
	Tpts.fillPointRange(30*30*30*30, xlo2, xhi2);
	Tpts.finalize();
	
	PointCloudHistogram* hADC[BOTH][nBetaTubes][2];
	for(Side s = EAST; s <= WEST; ++s)
		for(unsigned int t=0; t<nBetaTubes; t++)
			for(unsigned int i=0; i<2; i++)
				hADC[s][t][i] = new PointCloudHistogram(&Tpts);
	
	// load and scan data
	printf("Loading data...\n");
	TrigTreeScanner T;
	T.addFile(getEnvSafe("UCNAOUTPUTDIR")+"/hists/pmt_trig_151*.root");
	T.startScan();
	while(T.nextPoint()) {
		for(Side s = EAST; s <= WEST; ++s) {
			// require 2-of-4 on opposite side
			if(!(T.trigFlags & (1 << (nBetaTubes+(nBetaTubes+1)*(1-s))))) continue;
			
			// require in data range
			unsigned int t=0;
			for(t=0; t<nBetaTubes; t++) if(!(xlo[t] <= T.adc[s][t] && T.adc[s][t] <= xhi[t])) break;
			if(t<nBetaTubes) continue;
			
			for(t=0; t<nBetaTubes; t++) {
				if(T.trigFlags & (1<<(t+(nBetaTubes+1)*s))) hADC[s][t][true]->Fill(T.adc[s]);
				else hADC[s][t][false]->Fill(T.adc[s]);
			}
		}
	}
	
	OutputManager OM("TrigEffic",getEnvSafe("UCNA_ANA_PLOTS")+"/test/TrigEffic");
	TH1F* hTrig[2];
	for(unsigned int i=0; i<2; i++) {
		hTrig[i] = new TH1F(("hTrig_"+itos(i)).c_str(),"Trigger Events", 100, -50, 150);
		hTrig[i]->GetXaxis()->SetTitle("projected ADC channels");
		hTrig[i]->SetLineColor(2+2*i);
	}

	// calculate optimum cut orientation
	ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Simplex");
	min->SetTolerance(1.5);
	min->SetPrintLevel(3);
	const double dzero[] = {0,0,0,0,0};
	for(Side s = EAST; s <= WEST; ++s) {
		for(unsigned int t=0; t<nBetaTubes; t++) {
			
			min->Clear();
			
			SeparationLineCutOptimizer SLCO(hADC[s][t][false], hADC[s][t][true], t);
			printf("Main channel only: separation error = %g\n", SLCO(dzero));
			// un-optimized separation
			/*
			for(unsigned int i=0; i<2; i++) {
				hTrig[i]->Reset();
				hADC[s][t][i]->project(SLCO.vProj, *hTrig[i]);
			}
			hTrig[false]->Draw();
			hTrig[true]->Draw("Same");
			drawVLine(SLCO.xdiv, OM.defaultCanvas);
			OM.printCanvas(sideSubst("PMT0_%c",s)+itos(t));
			*/
			
			ROOT::Math::Functor f(SLCO, nBetaTubes-1);
			min->SetFunction(f);
			
			char vname[] = "a";
			for(unsigned int i=0; i<nBetaTubes-1; i++) {
				vname[0] = 'x'+i;
				min->SetLimitedVariable(i, vname, 0., 0.05, -0.5, 0.5);
			}
			
			min->Minimize();
			printf("Cut optimized!\n");
			
			const double* xopt = min->X();
			assert(xopt);
			printf("\nOptimum separation:");
			for(unsigned int i=0; i<nBetaTubes-1; i++) printf("\t%g",xopt[i]);
			printf("\n\n");
			
			// separation graph plot
			SLCO(xopt);
			SLCO.gSep->Draw("AP");
			SLCO.gSep->SetTitle("Trigger threshold optimization");
			SLCO.gSep->GetXaxis()->SetTitle("trigger cut [ADC channels]");
			SLCO.gSep->GetYaxis()->SetTitle("miscategorized events");
			SLCO.gSep->Draw("AP");
			OM.printCanvas(sideSubst("PMT_Sep_%c",s)+itos(t));
			
			// separated histograms plot
			for(unsigned int i=0; i<2; i++) {
				hTrig[i]->Reset();
				hADC[s][t][i]->project(SLCO.vProj, *hTrig[i]);
			}
			hTrig[false]->Draw();
			hTrig[true]->Draw("Same");
			drawVLine(SLCO.xdiv, OM.defaultCanvas);
			OM.printCanvas(sideSubst("PMT_%c",s)+itos(t));
			
			
		}
	}
	
	return 0;
}
