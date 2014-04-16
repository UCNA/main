#include "ucnaAnalyzerBase.hh"

/// mini-analyzer for PMT trigger region events
class TriggerTreeScanner: public ucnaAnalyzerBase {
public:
	/// constructor
	TriggerTreeScanner(RunNum R, std::string bp, CalDB* CDB): ucnaAnalyzerBase(R, bp, "pmt_trig", CDB) {
		
		// set up input
		
		// PMTs and TDCs, in quadrant order (+x towards SCS, +y up, +z from East to West)
		const int pmt_adc_nums[] = {2,3,0,1,4,5,6,7};
		for(Side s = EAST; s <= WEST; ++s)
			for(unsigned int t=0; t<nBetaTubes; t++)
				SetBranchAddress("Qadc"+itos(pmt_adc_nums[t+4*s]),&r_PMTADC[s][t]);
		
		// set up output
		openOutfile();
		TTrig = (TTree*)addObject(new TTree("trig","trigger-region scintillator events"));
		TTrig->SetMaxVirtualSize(1000000);
		
		TTrig->Branch("Sis00",&SIS00,"Sis00/I");
		TTrig->Branch("trigFlags",&pmtTrigFlags,"trigFlags/I");
		for(Side s = EAST; s <= WEST; ++s) {
			//TTrig->Branch(sideSubst("TDC%c",s).c_str(), &r_PMTTDC[s], "t[5]/F");
			TTrig->Branch(sideSubst("ADC%c",s).c_str(), &r_PMTADC[s], "q[4]/F");
			TTrig->Branch("nFiring", nf, "nFiring[2]/I");
		}
	}
	
	/// scan the data
	void analyze() {
		startScan();
		while(nextPoint()) {
		
			calibrateTimes();
			if(!passesBeamCuts()) continue;					// skip beam-burst events also vetoed in normal data
			if(!(trig2of4(EAST)||trig2of4(WEST))) continue;	// events with at least 1 2-fold trigger (use on opposite side)
			
			SIS00 = int(r_Sis00);
			pmtTrigFlags = 0;
			for(Side s = EAST; s <= WEST; ++s) {
				for(unsigned int t=0; t<nBetaTubes; t++) if(pmtFired(s,t)) pmtTrigFlags |= 1 << (t+(nBetaTubes+1)*s);
				if(trig2of4(s)) pmtTrigFlags |= 1 << (nBetaTubes + (nBetaTubes+1)*s);
				nf[s] = nFiring(s);
			}
			for(Side s = EAST; s <= WEST; ++s)
				PCal.pedSubtract(s, r_PMTADC[s], fTimeScaler[BOTH]);
		
			TTrig->Fill();
		}
		
		printf("Saved %lli events.\n", TTrig->GetEntries());
		setWriteRoot(true);
		write();
	}
	
	
	Int_t SIS00;							///< integer version of SIS00 trigger flags
	Int_t pmtTrigFlags;						///< PMT trigger flags
	Float_t r_PMTADC[BOTH][nBetaTubes];		///< PMT ADCs
	Int_t nf[BOTH];							///< number of PMTs firing on each side
	TTree* TTrig;							///< trigger events output tree
	
};


int main(int argc, char** argv) {
	
	if(argc<2) {
		printf("Use: %s <run number>\n", argv[0]);
		exit(1);
	}
	
	RunNum rn = atoi(argv[1]);
	std::string outDir = getEnvSafe("UCNAOUTPUTDIR");
	std::string inDir = getEnvSafe("UCNADATADIR");
	
	TriggerTreeScanner A(rn, outDir, CalDBSQL::getCDB());
	A.addFile(inDir+"/full"+itos(rn)+".root");
	A.analyze();
	
	return 0;
}
