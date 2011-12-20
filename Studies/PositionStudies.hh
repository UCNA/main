#ifndef POSITIONSTUDIES_HH
#define POSITIONSTUDIES_HH 1

#include "ProcessedDataScanner.hh"
#include "PostOfficialAnalyzer.hh"

/// gather Anode energy calibration data
void AnodeCalibration(ProcessedDataScanner& PDS, OutputManager& OM, unsigned int nRings = 6);

/// gather Cathode calibration data
void CathodeCalibration(PostOfficialAnalyzer& POA, OutputManager& OM);

#endif
