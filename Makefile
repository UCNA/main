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

CXXFLAGS = -O3 -m32 -Wall `root-config --cflags` \
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

Utils = ControlMenu.o strutils.o PathUtils.o TSpectrumUtils.o QFile.o GraphUtils.o MultiGaus.o TagCounter.o \
	Enums.o Types.o Octet.o SpectrumPeak.o Source.o SQL_Utils.o GraphicsUtils.o OutputManager.o RData.o

Detectors = WirechamberReconstruction.o

Calibration = PositionResponse.o SimNonlinearity.o PMTGenerator.o \
	EnergyCalibrator.o WirechamberCalibrator.o CalDBSQL.o SourceDBSQL.o GainStabilizer.o EvisConverter.o ManualInfo.o
	
Analysis = TChainScanner.o ProcessedDataScanner.o PostAnalyzer.o PostOfficialAnalyzer.o G4toPMT.o TH1toPMT.o DataSource.o \
	KurieFitter.o EndpointStudy.o ReSource.o EfficCurve.o BetaSpectrum.o

Studies = PlotMakers.o SRAsym.o PositionStudies.o SegmentSaver.o RunAccumulator.o OctetAnalyzer.o AsymmetryAnalyzer.o

objects = $(Utils) $(Detectors) $(Calibration) $(Analysis) $(Studies)

all:
	make UCNAnalyzer

UCNAnalyzer: Analyzer.cpp $(objects)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) Analyzer.cpp $(objects) -o UCNAnalyzer

CalibratorExample: CalibratorExample.cpp $(objects)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) CalibratorExample.cpp $(objects) -o CalibratorExample

DataScannerExample: DataScannerExample.cc $(objects)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) DataScannerExample.cc $(objects) -o DataScannerExample

OctetAnalyzerExample: OctetAnalyzerExample.cc $(objects)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) Studies/OctetAnalyzerExample.cc $(objects) -o OctetAnalyzerExample
	
ExtractFierzTerm: ExtractFierzTerm.cc $(objects)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) ExtractFierzTerm.cc $(objects) -o ExtractFierzTerm
	
QCalc: FPNCalc.cc $(Utils)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) IOUtils/FPNCalc.cc $(objects) -o QCalc
	
FierzOctetAnalyzer: FierzOctetAnalyzer.cc $(objects)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) Studies/FierzOctetAnalyzer.cc $(objects) -o FierzOctetAnalyzer
	
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
	-rm -f UCNAnalyzer OctetAnalyzerExample DataScannerExample CalibratorExample ExtractFierzTerm Analyzer
	-rm -f *.o
	-rm -rf *.dSYM
	-rm -rf latex/
	-rm -rf html/
	
