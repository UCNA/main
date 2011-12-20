#ifndef PEAKTRACKER_HH
#define PEAKTRACKER_HH 1

#include <vector>
#include <string>
#include <TGraphErrors.h>
#include "Enums.hh"

/// load RunMonitor data from several runs, optionally creating merged plot
std::vector<TGraphErrors*> mergeloader(const std::vector<RunNum>& druns, std::string sname, std::string infotype, bool domerged);

/// track LED peaks, group runs together by time and discontinuities for later Co60/muon tracking 
void trackPeaks(RunNum startRun, RunNum endRun, bool makePlots = false, float tsplit = 3600);

/// make plots tracking GMS correction
void trackGMS(RunNum startRun, RunNum endRun);

/// track Co60 peaks in ref PMT
void trackCo60(RunNum startRun0, RunNum endRun0, Side s);

/// track muon peak positions in Beta PMTs
//void trackMuonPeak(RunNum startRun0, RunNum endRun0, Side s);

/// create wirechamber pedestal history plots for run range
void trackPeds(RunNum startRun, RunNum endRun);

#endif
