#include "EndpointStudy.hh"
#include "CalDBSQL.hh"
#include "CalDBFake.hh"
#include "MultiGaus.hh"
#include <TStyle.h>

PositionBinner::PositionBinner(OutputManager* pnt, const std::string& nm, float r, unsigned int nr, const std::string& infl):
RunAccumulator(pnt,nm,infl), sects(nr,r) {
	
	// load sector data
	if(fIn) {
		QFile qOld(inflname+".txt");
		Stringmap sct = qOld.getFirst("SectorCutter");
		sects = SectorCutter(int(sct.getDefault("nRings",0)),sct.getDefault("radius",0));
		assert(sects.n && sects.r);
	}
	
	// save sector data to file
	Stringmap ms;
	ms.insert("nRings",sects.n);
	ms.insert("radius",sects.r);
	ms.insert("nSectors",sects.nSectors());
	qOut.insert("SectorCutter",ms);
	
	// set up histograms
	for(Side s = EAST; s <= WEST; s = nextSide(s)) {
		energySpectrum[s] = registerFGBGPair("hEnergy","Combined Energy",200,-100,1000,AFP_OTHER,s);
		for(unsigned int m=0; m<getNSectors(); m++) {
			sectAnode[s].push_back(registerFGBGPair(std::string("hAnode_Sector_")+itos(m),std::string("Sector ")+itos(m)+" Anode",200,0,4000,AFP_OTHER,s));
			for(unsigned int t=0; t<nBetaTubes; t++)
				sectEnergy[s][t].push_back(registerFGBGPair(std::string("hTuben_")+itos(t)+"_Sector_"+itos(m),
															std::string("Sector ")+itos(m)+" Tube "+itos(t)+" Energy",200,-100,1000,AFP_OTHER,s));
		}
	}
}

void PositionBinner::fillCoreHists(ProcessedDataScanner& PDS, double weight) {
	const Side s = PDS.fSide;
	if(!(PDS.fType == TYPE_0_EVENT && PDS.fPID == PID_BETA && (s==EAST||s==WEST))) return;
	unsigned int m = getSector(PDS.wires[s][X_DIRECTION].center,PDS.wires[s][Y_DIRECTION].center);
	if(m>=getNSectors()) return;
	for(unsigned int t=0; t<nBetaTubes; t++)
		sectEnergy[s][t][m].h[currentFG]->Fill(PDS.scints[s].tuben[t].x,weight);
	if(PDS.radius2(s) <= 25*25)
		energySpectrum[s].h[currentFG]->Fill(PDS.scints[s].energy.x,weight);
	if(PDS.getEtrue() > 225)
		sectAnode[s][m].h[currentFG]->Fill(PDS.mwpcs[s].anode,weight);
}


void PositionBinner::fitEndpoints() {
	for(unsigned int m=0; m<getNSectors(); m++) {
		for(unsigned int t=0; t<nBetaTubes; t++) {
		}
	}
}

/*		
		printf("Starting job %i [%i]...\n",j.job_id,j.thread_id);
		
		// center position of this sector
		float x,y;
		unsigned int n = j.job_id;
		sects->sectorCenter(n,x,y);
		printf("Kurie plot for sector %i@(%.1f,%.1f)\n",n,x,y);
		
		// generate simulation spectrum
		PMTGenerator PGen(*BaseGen);
		TH1F* simSpectrum[nBetaTubes+1];
		for(unsigned int t=0; t<nBetaTubes+1; t++) {
			simSpectrum[t] = new TH1F(*refSpectrum);
			zero(simSpectrum[t]);
			simSpectrum[t]->SetLineColor(t+1);
		}		
		PGen.setSide(s);
		PGen.setPosition(x, y);
		assert(false); //TODO
		//PGen.simSpectrum(refSpectrum,400000);
		
		// process each PMT
		for(unsigned int t=0; t<=nBetaTubes; t++) {
			
			pthread_mutex_lock(&output_lock);
			
			printf("\tTube %i\n",t);
			
			float eta = RefPos->eval(s, t, x, y, true);
			Stringmap m;
			m.insert("side",ctos(sideNames(s)));
			m.insert("tube",t);
			m.insert("sector",n);
			m.insert("x",x);
			m.insert("y",y);
			m.insert("eta_prev",eta);
			
			// finalize histograms
			SpectrumHistos* SH = secthists[s*(nBetaTubes+1)+t][n];
			SH->finishHistos(FG->totalTime.t[s], BG->totalTime.t[s]);
			
			// find feature
			float rscale = SH->diff->GetBinCenter(SH->diff->GetMaximumBin());
			TF1* gfit = new TF1("gfit","gaus",0.5*rscale,2*rscale);
			iterGaus(SH->diff,gfit,4,rscale,0.5*rscale,0.75,0);
			float_err feature(gfit->GetParameter(1),gfit->GetParError(1));
			float fheight = gfit->GetParameter(0);
			// find simulation peak
			simSpectrum[t]->Fit(gfit,"Q","",80,150);
			simSpectrum[t]->Scale(fheight/gfit->GetParameter(0));
			iterGaus(simSpectrum[t],gfit,4,110,50,0.75,0.0);
			float_err featurePos(gfit->GetParameter(1),gfit->GetParError(1));
			
			// 492keV peak in combined data
			if(t==nBetaTubes && find492pk) {
				TF1* fit492 = new TF1("fit492","gaus",400,500);
				iterGaus(SH->diff, fit492, 4, 450.0, 75.0, 0.50, 0);
				float_err pk492(fit492->GetParameter(1),fit492->GetParError(1));
				iterGaus(simSpectrum[t], fit492, 4, 450.0, 75.0, 0.50, 0);
				float_err sim492(fit492->GetParameter(1),fit492->GetParError(1));
				Stringmap m492(m);
				m492.insert("light",pk492.x*eta);
				m492.insert("energy",sim492.x);
				m492.insert("denergy",sim492.err);
				m492.insert("dlight",pk492.err*eta);
				m492.insert("feature","Xe492Peak");
				OM->qOut.insert("PosmapPoint",m492);
			}
			
			// 915keV endpoint
			float_err ep915, epSim;
			if(find915ep) {
				TGraphErrors* tgData = NULL;
				ep915 = kurieIterator(SH->diff,875.0,&tgData,915.,400,800);
				TGraphErrors* tgSim = NULL;
				epSim = kurieIterator(simSpectrum[t],875.0,&tgSim,915.,400,800);
				Stringmap m915(m);
				m915.insert("light",ep915.x*eta);
				m915.insert("feature_en",ep915.x);
				m915.insert("dlight",ep915.err*eta);
				m915.insert("energy",epSim.x);
				m915.insert("denergy",epSim.err);
				m915.insert("feature","Xe915EP");
				OM->qOut.insert("PosmapPoint",m915);
				if(t<nBetaTubes) {
					(*pinf)[nBetaTubes*s+t].adc.push_back(ep915.x*eta);
					(*pinf)[nBetaTubes*s+t].energy.push_back(epSim.x);
				}
				if(tgData && tgSim) {
					SH->defaultCanvas->cd();
					tgSim->SetMarkerColor(4);
					tgSim->Draw("AP");
					tgSim->GetXaxis()->SetLimits(0,1200);
					tgSim->Draw("AP");
					tgData->SetMarkerColor(2);
					tgData->Draw("P");
					SH->printCanvas(sideSubst("SectionKurie_%c",s)+itos(t)+"_"+itos(n));
				}
				if(tgData) delete(tgData);
				if(tgSim) delete(tgSim);
			}
			
					
			// scale out ref posmap; save results
			Stringmap m1(m);
			m1.insert("light",feature.x*eta);
			m1.insert("energy",featurePos.x);
			m1.insert("denergy",featurePos.err);
			m1.insert("dlight",feature.err*eta);
			m1.insert("feature","XeLowPeak");
			OM->qOut.insert("PosmapPoint",m1);
			
				
		} // tube t
	} // run()
};
*/
