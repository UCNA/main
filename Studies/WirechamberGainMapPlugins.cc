#include "WirechamberGainMapPlugins.hh"
#include "Sim2PMT.hh"
#include <time.h>
#include <TProfile.h>

WirechamberGainMapPluginBase::WirechamberGainMapPluginBase(RunAccumulator* RA, unsigned int nr, const std::string& nm): PositionBinnedPlugin(RA,nm,nr), myChgPrx(CHARGE_PROXY_NONE) {
	TH1F hTemplate(("hEw_"+nm).c_str(),"Wirechamber energy",50,0,5);
	for(Side s = EAST; s <= WEST; ++s)
		sectHists[s] = allocateSegmentHistograms(hTemplate,AFP_OTHER,s);
	TProfile pTemplate(("pWGain_"+nm).c_str(),"Average applied gain correction",sects.nSectors(),-0.5,sects.nSectors()-0.5);
	for(Side s = EAST; s <= WEST; ++s) {
		sectGains[s] = registerFGBGPair(pTemplate, AFP_OTHER, s);
		sectGains[s]->doSubtraction = false;
	}
}

bool WirechamberGainMapPluginBase::evtLoc(const ProcessedDataScanner& PDS, Side& s, unsigned int& m, float& x, float& y) const {
	s = PDS.fSide;
	if(!(PDS.fType == TYPE_0_EVENT && PDS.fPID == PID_BETA && 400 <= PDS.getEnergy() && PDS.getEnergy() <= 600)) return false;
	x = PDS.wires[s][X_DIRECTION].center;
	y = PDS.wires[s][Y_DIRECTION].center;
	m = sects.sector(x,y);
	return m<sects.nSectors();
}

void WirechamberGainMapPluginBase::fillCoreHists(ProcessedDataScanner& PDS, double weight) {
	Side s; unsigned int m; float x,y,x0,y0;
	if(evtLoc(PDS,s,m,x,y)) {
		const MWPC_Ecal_Spec& mes = PDS.ActiveCal->getAltEcal(s,myChgPrx);
		
		double Q = PDS.ActiveCal->chargeProxy(s,myChgPrx,PDS.wires[s][X_DIRECTION],PDS.wires[s][Y_DIRECTION],PDS.mwpcs[s]);
		sectHists[s][m]->h[currentGV]->Fill(Q*mes.gain_factor/mes.pcorr->eval(s,0,x,y,true),weight);
		
		sects.sectorCenter(m,x0,y0);
		((TProfile*)(sectGains[s]->h[currentGV]))->Fill(m,mes.gain_factor/mes.pcorr->eval(s,0,x0,y0,true),weight);
	}
}

void WirechamberGainMapPluginBase::genPosmap(const std::string& pmapNameBase) const {
	// name for position map
	std::string pmapname = pmapNameBase+"_"+itos(myA->runCounts.counts.begin()->first)+"-"+itos(myA->runCounts.counts.rbegin()->first)+"/"+itos(time(NULL));
	// Calibration DB to write output to
	CalDBSQL* CDBout = CalDBSQL::getCDB(false);
	assert(CDBout);
	// reserve new position map ID number
	unsigned int pmid = CDBout->newPosmap(pmapname,sects.n,sects.r);
	
	// set up Landau function fitter
	const TH1* hTemplate = sectHists[EAST][0]->h[GV_OPEN];
	const double xmn = hTemplate->GetXaxis()->GetXmin();
	const double xmx = hTemplate->GetXaxis()->GetXmax();
	TF1 fLandau("landauFit","landau",xmn,xmx);
	fLandau.SetParameter(2,0.5);
	fLandau.SetParLimits(2,0.2,3.0);
	
	// fit each position bin
	for(Side s=EAST; s<=WEST; ++s) {
		for(unsigned int m=0; m<sects.nSectors(); m++) {
			printf("\n------------- %s %i ------------------\n",sideWords(s),m);
			float x,y;
			sects.sectorCenter(m,x,y);
			TH1* hSector = sectHists[s][m]->h[GV_OPEN];
			double c = hSector->GetBinCenter(hSector->GetMaximumBin());
			double mx = hSector->GetMaximum();
			printf("Max = %.3g at %.1f\n",mx,c);
			fLandau.SetParameter(0,mx*c/fLandau.GetParameter(2));
			fLandau.SetParameter(1,c);
			hSector->Fit(&fLandau,"QEMBN");
			c = fLandau.GetParameter(1);
			double w = fLandau.GetParameter(2);
			fLandau.SetRange(c-2*w,c+2*w);
			hSector->Fit(&fLandau,"REM");
			
			CDBout->addPosmapPoint(pmid,s,0,m,fLandau.GetParameter(1),sectGains[s]->h[GV_OPEN]->GetBinContent(m+1),x,y);
		}
	}
	
	// save posmap name/ID to output file
	Stringmap pmsm;
	pmsm.insert("pmid",pmid);
	pmsm.insert("name",pmapname);
	myA->qOut.insert(pmapNameBase+"_posmap",pmsm);
}

//------------------------------------

void WirechamberEdepMapPlugin::fillCoreHists(ProcessedDataScanner& PDS, double weight) {
	Side s; unsigned int m; float x,y;
	if(evtLoc(PDS,s,m,x,y)) {
		Sim2PMT* S2P = (Sim2PMT*)&PDS;
		sectHists[s][m]->h[currentGV]->Fill(S2P->eW[s],weight);
		((TProfile*)(sectGains[s]->h[currentGV]))->Fill(m,1.0,weight);
	}
}

