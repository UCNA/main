#include "DataSource.hh"

ProcessedDataScanner* getDataSource(InputDataSource src, bool withCalibrators) {
	if(src==INPUT_UNOFFICIAL)
		return new PostAnalyzer(withCalibrators);
	if(src==INPUT_OFFICIAL)
		return new PostOfficialAnalyzer(withCalibrators);
	return NULL;
}
