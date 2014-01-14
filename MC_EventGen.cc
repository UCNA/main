#include "NuclEvtGen.hh"
#include "ControlMenu.hh"
#include "PathUtils.hh"
#include <Math/QuasiRandom.h>
#include <Math/Random.h>
#include <TFile.h>
#include <TTree.h>


using namespace ROOT::Math;

void mi_evtgen(std::deque<std::string>&, std::stack<std::string>& stack) {

	// load arguments
	const unsigned int nTrees = streamInteractor::popInt(stack);
	const unsigned int nPerTree = streamInteractor::popInt(stack);
	const std::string rtSelect = streamInteractor::popString(stack);
	const std::string vpSelect = streamInteractor::popString(stack);
	std::string outPath = streamInteractor::popString(stack);
	const std::string genName = streamInteractor::popString(stack);
	
	// load generators
	static NucDecayLibrary NDL(getEnvSafe("UCNA_AUX")+"/NuclearDecays/",1e-6);
	NucDecaySystem& NDS = NDL.getGenerator(genName);
	NDS.display();
	PositionGenerator* PosGen = vpSelect=="f" ?	(PositionGenerator*)(new CylPosGen(3.,2.3*0.0254)) :
								vpSelect=="g" ?	(PositionGenerator*)(new CylPosGen(4.3,.075)) :
								vpSelect=="c" ?	(PositionGenerator*)(new CubePosGen()) :
												(PositionGenerator*)(new FixedPosGen());
	enum QRndType {
		INDEP_RANDOM,
		QR_SOBOL,
		QR_NIED
	} qrt = (rtSelect=="s") ? QR_SOBOL : (rtSelect=="n"?QR_NIED:INDEP_RANDOM);
		
	outPath = outPath+"/"+genName+"_"+vpSelect+"_"+rtSelect+"/";
	makePath(outPath);
	
	// set up variables
	unsigned int evtn = 0;
	double vpos[Z_DIRECTION+1];
	for(AxisDirection d = X_DIRECTION; d <= Z_DIRECTION; ++d) vpos[d] = 0;
	const unsigned int posDF = PosGen ? PosGen->getNDF(): 0;
	const unsigned int decayDF = NDS.getNDF();
	const unsigned int totDF = posDF+decayDF+1;
	std::vector<double> rnd(totDF);
	NucDecayEvent tEvt;
	RandomMT r0;
	QuasiRandomSobol rSobol(totDF);
	QuasiRandomNiederreiter rNied(totDF);
	if(qrt==QR_NIED && totDF>12) {
		printf("Warning: ROOT's Niederreiter quasi-random generator invalid for >12 DF; switching to Sobol\n");
		qrt=QR_SOBOL;
	}
	
	printf("Generating events in '%s' with %i+%i+1 random DF\n",outPath.c_str(),posDF,decayDF);
	
	for(unsigned int tn=0; tn<nTrees; tn++) {
	
		printf("Tree %i/%i: %i events\n",tn+1,nTrees,nPerTree);
	
		TFile f((outPath+"/Evts_"+itos(tn)+".root").c_str(),"RECREATE");
		f.cd();
		
		TTree T("Evts","MC initial events");
		T.Branch("num",&tEvt.eid,"num/I");
		T.Branch("PID",&tEvt.d,"PID/I");
		T.Branch("KE",&tEvt.E,"KE/D");
		T.Branch("vertex",tEvt.x,"vertex[3]/D");
		T.Branch("direction",tEvt.p,"direction[3]/D");
		T.Branch("time",&tEvt.t,"time/D");
		T.Branch("weight",&tEvt.w,"weight/D");
		
		for(unsigned int i=0; i<nPerTree; i++) {
			std::vector<NucDecayEvent> evts;
			if(qrt==INDEP_RANDOM) r0.RndmArray(totDF,&rnd[0]);
			else if(qrt==QR_SOBOL) rSobol.Next(&rnd[0]);
			else if(qrt==QR_NIED) rNied.Next(&rnd[0]);
			
			NDS.genDecayChain(evts, &rnd[0]);
			if(PosGen) PosGen->genPos(vpos,&rnd[decayDF+1]);
			for(unsigned int i=0; i<evts.size(); i++) {
				tEvt = evts[i];
				tEvt.eid = evtn;
				for(AxisDirection d = X_DIRECTION; d <= Z_DIRECTION; ++d) tEvt.x[d] = vpos[d];
				T.Fill();
			}
			evtn++;
		}
		
		T.Write();
		f.Close();
	}
}


int main(int argc, char *argv[]) {

	inputRequester exitMenu("Exit Menu",&menutils_Exit);
	
	// discrete selectors
	NameSelector selectRandomType("Random Source");
	selectRandomType.addChoice("Pseudo-Random","p");
	selectRandomType.addChoice("Sobol","s");
	selectRandomType.addChoice("Niederretier","n");
	selectRandomType.setDefault("n");
	
	NameSelector selectVertexPos("Vertex Positioner");
	selectVertexPos.addChoice("decay trap fiducial","f");
	selectVertexPos.addChoice("origin","o");
	selectVertexPos.addChoice("spectrometer gas fill","g");
	selectVertexPos.addChoice("unit cube","c");
	selectVertexPos.setDefault("f");
	
	// event generator routine
	inputRequester run_evt_gen("Run event generator",&mi_evtgen);
	run_evt_gen.addArg("Generator name");
	run_evt_gen.addArg("Output path",getEnvSafe("G4EVTDIR"));
	run_evt_gen.addArg("","",&selectVertexPos);
	run_evt_gen.addArg("","",&selectRandomType);
	run_evt_gen.addArg("Events per TTree","10000");
	run_evt_gen.addArg("N. TTrees","100");
		
	// main menu
	OptionsMenu OM("Event Generator Menu");
	OM.addChoice(&run_evt_gen,"run");
	OM.addChoice(&exitMenu,"x");
	
	// load command line arguments
	std::deque<std::string> args;
	for(int i=1; i<argc; i++)
		args.push_back(argv[i]);
	std::stack<std::string> stack;
	OM.doIt(args,stack);
	
	printf("\n\n\n>>>>> Goodbye. <<<<<\n\n\n");
	return 0;
}
