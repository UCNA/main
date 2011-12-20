#ifndef DATASOURCE_HH
#define DATASOURCE_HH 1

#include "ProcessedDataScanner.hh"
#include "PostAnalyzer.hh"
#include "PostOfficialAnalyzer.hh"
#include "G4toPMT.hh"

/// sources for input data
enum InputDataSource {
	INPUT_NULL,			//< null/unknown data source
	INPUT_UNOFFICIAL,	//< unofficial replay points
	INPUT_OFFICIAL,		//< official replay data
	INPUT_GEANT4		//< geant4 simulation data for a run
};

/// get a data source with the approriate specifications
ProcessedDataScanner* getDataSource(InputDataSource src, bool withCalibrators);

#endif
