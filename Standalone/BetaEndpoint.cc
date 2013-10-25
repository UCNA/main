#include "PMTGenerator.hh"
#include "EnergyCalibrator.hh"
#include "G4toPMT.hh"
#include "PathUtils.hh"
#include "StyleSetup.hh"
#include "GraphicsUtils.hh"
#include "KurieFitter.hh"
#include <algorithm>
#include <TGraphErrors.h>

class PMTCalibratorMeddler: public PMTCalibrator {
public:
	PMTCalibratorMeddler(RunNum rn): PMTCalibrator(rn) {}
	void setDeltaL(Side s, unsigned int t, double x) { deltaL[s][t] = x; } // = delta E / sqrt(E)
};

/// iterate Kurie plot until convergence is reached
//float_err kurieIterator(TH1* spectrum, float iguess, TGraphErrors** tgout = NULL, float targetEP = neutronBetaEp,
//						float fitStart = 250, float fitEnd = 700, unsigned int nmax = 50, float etol = 0.02 );


double avg_pe(const PMTCalibrator& PCal, Side s, unsigned int t) {
	double p=0;
	int n=0;
	for(double x = -60; x<=60; x+=2) {
		for(double y=-60; y<=60; y+=2) {
			if(x*x+y*y<55*55) {
				p += PCal.nPE(s,t,1000,x,y);
				n++;
			}
		}
	}
	return p/n;
}

void NeutronEndpointSmearShift() {

	OutputManager OM("EndpointSmear",getEnvSafe("UCNA_ANA_PLOTS")+"/SimSpectrum/EndpointSmear/");
	OM.defaultCanvas->cd();
	OM.defaultCanvas->SetLeftMargin(0.12);
	PMTCalibratorMeddler PCal(16194);
	G4toPMT G2P;
	G2P.setCalibrator(PCal);
	G2P.addFile(getEnvSafe("G4OUTDIR")+"/thinFoil_neutronBetaUnpol/analyzed_*.root");

	for(Side s = EAST; s<=WEST; ++s) {
		for(unsigned int t=0; t<nBetaTubes; t++) {
			Stringmap m;
			m.insert("side",sideWords(s));
			m.insert("tube",t);
			m.insert("PEperMEVc",PCal.nPE(s,t,1000,0,0));
			m.insert("PEperMEV",avg_pe(PCal,s,t));
			OM.qOut.insert("PMTspec",m);
		}
	}
	OM.write();

	double dL0 = 100/sqrt(1000); //< middle-range PMT energy resolution
	double r0 = 0.5;
	double r1 = 2.0;
	unsigned int npts = 3;
	
	double fitStart = 150;
	unsigned int cycNum = 0;
	double fitEnd = 635;

	for(unsigned int j=0; j<100; j++) {
		for(unsigned int i=0; i<npts; i++) {
			// set range of resolutions for each PMT, and matching histograms
			unsigned int nn = 0;
			TH1F* hSpec[2][nBetaTubes];
			std::vector<TH1*> hToPlot;
			for(Side s = EAST; s<=WEST; ++s) {
				for(unsigned int t=0; t<nBetaTubes; t++) {
					double l = r0+(r1-r0)*double(8*i+nn)/double(8*npts-1);
					PCal.setDeltaL(s,t, l*dL0);
					nn++;
					hSpec[s][t] = OM.registeredTH1F("hSpec", sideSubst("betaSpec_%c",s)+itos(t)+"_"+itos(cycNum), 200, 0, 1000);
					hToPlot.push_back(hSpec[s][t]);
				}
			}
			cycNum++;
			
			// simulate spectrum events
			G2P.resetSimCounters();
			G2P.nToSim = 100000;
			G2P.startScan(true);
			while(G2P.nSimmed < G2P.nToSim) {
				G2P.nextPoint();
				if(!(G2P.fType == TYPE_0_EVENT && G2P.fPID == PID_BETA)) continue;
				Side s = G2P.fSide;
				for(unsigned int t=0; t<nBetaTubes; t++)
					hSpec[s][t]->Fill(G2P.scints[s].tuben[t].x, G2P.physicsWeight);
			}
			
			for(fitStart = 80; fitStart <= 240; fitStart += 20) {
			//for(fitEnd = 600; fitEnd <= 720; fitEnd += 20) {
				// fit endpoints
				for(Side s = EAST; s<=WEST; ++s) {
					for(unsigned int t=0; t<nBetaTubes; t++) {
						float_err ep = kurieIterator(hSpec[s][t], 800., NULL, neutronBetaEp, fitStart, fitEnd);
						Stringmap m;
						m.insert("fitStart",fitStart);
						m.insert("fitEnd",fitEnd);
						m.insert("endpoint",ep.x);
						m.insert("dendpoint",ep.err);
						m.insert("counts",hSpec[s][t]->Integral());
						m.display("--- ");
						m.insert("PEperMEVc",PCal.nPE(s,t,1000,0,0));
						m.insert("PEperMEV",avg_pe(PCal,s,t));
						OM.qOut.insert("kurieFit",m);
					}
				}
			}
		}
		OM.write();
	}
}


void XenonEnergyShift() {
	OutputManager OM("EndpointEnergyShift",getEnvSafe("UCNA_ANA_PLOTS")+"/SimSpectrum/XenonEndpoint/");
	OM.defaultCanvas->cd();
	OM.defaultCanvas->SetLeftMargin(0.12);
	PMTCalibrator PCal(16194);
	
	double fitStart = 450;
	double fitEnd = 750;
	unsigned int cycNum = 0;
	const unsigned int nShifts = 9;
	double shifts[nShifts];
	for(unsigned int sh=0; sh<nShifts; sh++)
		shifts[sh] = (float(sh)/float(nShifts-1)-0.5)*0.10;
	G4toPMT G2P;
	G2P.setCalibrator(PCal);
	G2P.addFile("/data2/mmendenhall/G4Out/2010/20120917_Xe135_3-2+/analyzed_*.root");
	
	for(unsigned int j=0; j<100; j++) {

		TH1F* hSpec[2][nBetaTubes][nShifts];
		for(Side s = EAST; s<=WEST; ++s)
			for(unsigned int t=0; t<nBetaTubes; t++)
				for(unsigned int sh=0; sh<nShifts; sh++)
				hSpec[s][t][sh] = OM.registeredTH1F("hSpec", sideSubst("betaSpec_%c",s)+itos(t)+"_"+itos(sh)+"_"+itos(cycNum), 200, 0, 1000);
		cycNum++;
		
		// simulate spectrum events
		G2P.resetSimCounters();
		G2P.nToSim = 500000;
		G2P.startScan(true);
		while(G2P.nSimmed < G2P.nToSim) {
			G2P.nextPoint();
			if(!(G2P.fType == TYPE_0_EVENT && G2P.fPID == PID_BETA)) continue;
			Side s = G2P.fSide;
			for(unsigned int sh=0; sh<nShifts; sh++) {
				double shft = 1+shifts[sh]*(G2P.eQ[s]-500.)/500.;
				for(unsigned int t=0; t<nBetaTubes; t++) {
					hSpec[s][t][sh]->Fill(G2P.scints[s].tuben[t].x, G2P.physicsWeight*shft);
				}
			}
		}
		
		// fit endpoints
		for(Side s = EAST; s<=WEST; ++s) {
			for(unsigned int t=0; t<nBetaTubes; t++) {
				for(unsigned int sh=0; sh<nShifts; sh++) {
					float_err ep = kurieIterator(hSpec[s][t][sh], 915., NULL, 915., fitStart, fitEnd);
					Stringmap m;
					m.insert("side",sideSubst("%c",s));
					m.insert("tube",itos(t));
					m.insert("fitStart",fitStart);
					m.insert("fitEnd",fitEnd);
					m.insert("endpoint",ep.x);
					m.insert("dendpoint",ep.err);
					m.insert("shift",shifts[sh]);
					m.insert("counts",hSpec[s][t][sh]->Integral());
					m.insert("fitcounts",hSpec[s][t][sh]->Integral(hSpec[s][t][sh]->FindBin(fitStart),hSpec[s][t][sh]->FindBin(fitEnd)));
					OM.qOut.insert("kurieFit",m);
				}
			}
		}
		OM.write();
	}
}

int main(int argc, char *argv[]) {

	ROOTStyleSetup();
	
	if(argc >= 2 && std::string(argv[1])=="xenon") {
		XenonEnergyShift();
	} else {
		NeutronEndpointSmearShift();
	}
}
