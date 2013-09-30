#include "PMTGenerator.hh"
#include "G4toPMT.hh"
#include "PathUtils.hh"
#include "StyleSetup.hh"
#include <algorithm>
#include <TGraphErrors.h>

int main(int argc, char *argv[]) {

	ROOTStyleSetup();
	
	RunNum rn = 16194;		// run number to re-simulate for
	double eMax = 1600;		// histogram maximum
	double nOrigEvts = 1e7;	// number of events in each simulation
	OutputManager OM("GammaChart",getEnvSafe("UCNA_ANA_PLOTS")+"/SimSpectrum/Gammas/");
	OM.defaultCanvas->cd();
	OM.defaultCanvas->SetLeftMargin(0.12);
	PMTCalibrator PCal(rn);
	
	std::vector<TH1F*> hSpec;
	std::vector<TH2F*> hPos;
	std::vector<double> eGamma;
	for(int l = 100; l<=1600; l *= 2) eGamma.push_back(l);
	for(int l = 300; l<=1200; l *= 2) eGamma.push_back(l);
	std::sort(eGamma.begin(),eGamma.end());
	
	for(unsigned int i=0; i<eGamma.size(); i++) {
		double l = eGamma[i];
		G4toPMT G2P;
		G2P.setCalibrator(PCal);
		G2P.addFile(getEnvSafe("G4WORKDIR")+"/output/SourceHolderGammas_eGunRandMomentum_"+itos(l)+".0keV/analyzed_*.root");
		hSpec.push_back(OM.registeredTH1F(std::string("hSpec_")+itos(l),"EventSpectrum",50,0,eMax));
		hSpec.back()->SetTitle("Source Foil Comptons");
		hSpec.back()->GetXaxis()->SetTitle("Energy [keV]");
		hSpec.back()->GetYaxis()->SetTitle("Events / keV / 10^{6} gammas");
		hSpec.back()->GetYaxis()->SetTitleOffset(1.2);
		
		hPos.push_back(OM.registeredTH2F(std::string("hPos_")+itos(l), "detected Compton electron positions", 100, -30, 30, 100, -30, 30));
		hPos.back()->GetXaxis()->SetTitle("x position [mm]");
		hPos.back()->GetYaxis()->SetTitle("y position [mm]");
		hPos.back()->GetYaxis()->SetTitleOffset(1.2);
		
		G2P.startScan(false);
		while(G2P.nextPoint()) {
			if(G2P.fType >= TYPE_IV_EVENT) continue;
			Side s = G2P.fSide;
			hPos.back()->Fill(G2P.wires[s][X_DIRECTION].center,G2P.wires[s][Y_DIRECTION].center);
			if(G2P.radius(s)<7.0) continue;
			hSpec.back()->Fill(G2P.getEtrue());
		}
	}
	
	TGraphErrors g((int)eGamma.size());
	double gScale = 1e6/nOrigEvts;
	for(unsigned int i=0; i<eGamma.size(); i++) {
		g.SetPoint(i,eGamma[i],hSpec[i]->Integral()*gScale);
		g.SetPointError(i, 0, sqrt(hSpec[i]->Integral())*gScale);
		hSpec[i]->Scale(1.0e6/nOrigEvts/hSpec[i]->GetBinWidth(1));
		hSpec[i]->GetYaxis()->SetRangeUser(0,4.);
		hSpec[i]->SetLineStyle(1+(i%3));
		hSpec[i]->SetLineWidth(2);
		hSpec[i]->Draw(i?"L SAME":"L");
	}
	OM.printCanvas("SourceComptons_Outer");
	
	g.Draw("AP");
	g.SetTitle("Compton production from sealed source gammas");
	g.GetXaxis()->SetTitle("gamma energy [keV]");
	g.GetYaxis()->SetTitle("detected events per 10^{6} gammas");
	g.GetYaxis()->SetTitleOffset(1.5);
	g.SetLineWidth(2);
	g.Draw("AP");
	OM.printCanvas("SourceComptonsNumber_Outer");
	
	for(unsigned int i=0; i<eGamma.size(); i++) {
		hPos[i]->Draw();
		OM.printCanvas("SimGammaPositions_"+itos(eGamma[i]));
	}
}
