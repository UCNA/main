#####################################################################
#
#  Name:         Makefile
#  Created by:   Michael Mendenhall
#
#  Contents:     Makefile for UCNA Analysis code
#
#####################################################################


# assure correct shell is used
SHELL = /bin/sh
# apply implicit rules only for listed file types
.SUFFIXES:
.SUFFIXES: .c .cc .cpp .o
	 
# compiler command to use
CC = cc
CXX = g++

CXXFLAGS = -g -O3 -m32 -Wall `root-config --cflags` \
	-I. -IIOUtils -IRootUtils -IBaseTypes -IDetectors -IMathUtils -ICalibration -IAnalysis -IStudies
LDFLAGS = `root-config --libs` -lSpectrum 

ifdef PROFILER_COMPILE
	CXXFLAGS += -pg
	LDFLAGS += -pg
endif

#
# things to build
#

VPATH = ./:IOUtils/:RootUtils/:BaseTypes/:Detectors/:MathUtils/:Calibration/:Analysis/:Studies/

Utils = ControlMenu.o strutils.o PathUtils.o TSpectrumUtils.o QFile.o CGraph.o GraphUtils.o MultiGaus.o TagCounter.o \
	Enums.o Types.o Octet.o SpectrumPeak.o Source.o SQL_Utils.o GraphicsUtils.o OutputManager.o RData.o AutoThreader.o

Detectors = RunManager.o Subsystem.o CoTracker.o Trigger.o BetaScint.o LEDTracker.o MuonVeto.o \
	WirechamberReconstruction.o Wirechamber.o MWPC.o Detector.o EventParser.o

Calibration = AnalysisDB.o RunsDB.o PositionResponse.o SimNonlinearity.o SimCalibrations.o \
	EnergyCalibrator.o CalDBSQL.o GainStabilizer.o EvisConverter.o ManualInfo.o
	
Analysis = TChainScanner.o ProcessedDataScanner.o PostAnalyzer.o PostOfficialAnalyzer.o G4toPMT.o DataSource.o \
	SpectrumHistos.o KurieFitter.o PeakTracker.o EndpointStudy.o Triad.o ReSource.o EfficCurve.o BetaSpectrum.o

Studies = PlotMakers.o SRAsym.o AsymHists.o

objects = $(Utils) $(Detectors) $(Calibration) $(Analysis) $(Studies)

all:
	make UCNAnalyzer

UCNAnalyzer: Analyzer.cpp $(objects)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) Analyzer.cpp $(objects) -o UCNAnalyzer
	
PostSimulator: PostSimulator.cpp $(objects)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) PostSimulator.cpp $(objects) -o PostSimulator

CalibratorExample: CalibratorExample.cpp $(objects)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) CalibratorExample.cpp $(objects) -o CalibratorExample

DataScannerExample: DataScannerExample.cc $(objects)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) DataScannerExample.cc $(objects) -o DataScannerExample
	
ExtractFierzTerm: ExtractFierzTerm.cc $(objects)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) ExtractFierzTerm.cc $(objects) -o ExtractFierzTerm
	
QCalc: FPNCalc.cc $(Utils)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) IOUtils/FPNCalc.cc $(objects) -o QCalc
	
#
# documentation via Doxygen
#

doc : latex/refman.pdf

latex/refman.pdf: latex/ 
	cd latex; make
latex/ : documentationConfig
	doxygen documentationConfig

#
# cleanup
#
.PHONY: clean
clean:
	-rm UCNAnalyzer $(objects)
	-rm -r *.dSYM
	-rm -r latex/
	-rm -r html/
	
