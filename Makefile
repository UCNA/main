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

CXXFLAGS = -O3 -Wall `root-config --cflags` -I. \
	-IIOUtils -IRootUtils -IBaseTypes -IMathUtils -ICalibration -IAnalysis -IStudies -IPhysics
LDFLAGS = `root-config --libs` -lSpectrum  -L. -lUCNA

ifdef PROFILER_COMPILE
	CXXFLAGS += -pg
	LDFLAGS += -pg
endif

#
# things to build
#

VPATH = ./:IOUtils/:RootUtils/:BaseTypes/:MathUtils/:Calibration/:Analysis/:Studies/:Physics/

Physics = BetaSpectrum.o ElectronBindingEnergy.o NuclEvtGen.o

Utils = ControlMenu.o strutils.o PathUtils.o TSpectrumUtils.o QFile.o GraphUtils.o MultiGaus.o TagCounter.o \
		SectorCutter.o LinHistCombo.o Enums.o Types.o FloatErr.o SMExcept.o Octet.o SpectrumPeak.o Source.o \
		SQL_Utils.o GraphicsUtils.o OutputManager.o RollingWindow.o RData.o

Calibration = PositionResponse.o PMTGenerator.o CathSegCalibrator.o WirechamberCalibrator.o \
		EnergyCalibrator.o CalDBSQL.o SourceDBSQL.o GainStabilizer.o EvisConverter.o
	
Analysis = TChainScanner.o RunSetScanner.o ProcessedDataScanner.o PostOfficialAnalyzer.o Sim2PMT.o G4toPMT.o \
		PenelopeToPMT.o TH1toPMT.o KurieFitter.o EndpointStudy.o ReSource.o EfficCurve.o AnalysisDB.o

Studies = PlotMakers.o SegmentSaver.o RunAccumulator.o OctetAnalyzer.o AsymmetryAnalyzer.o WirechamberStudy.o LEDScans.o

objects = $(Utils) $(Calibration) $(Analysis) $(Studies) $(Physics)

all:
	make UCNAnalyzer
	
libUCNA.a: $(objects)
	ar rs libUCNA.a $(objects)
	
UCNAnalyzer: Analyzer.cpp libUCNA.a
	$(CXX) $(CXXFLAGS) Analyzer.cpp -o UCNAnalyzer -L./ -lUCNA $(LDFLAGS)

CalibratorExample: CalibratorExample.cpp libUCNA.a
	$(CXX) $(CXXFLAGS) CalibratorExample.cpp  -o CalibratorExample $(LDFLAGS)

DataScannerExample: DataScannerExample.cc libUCNA.a
	$(CXX) $(CXXFLAGS) DataScannerExample.cc -o DataScannerExample $(LDFLAGS)

OctetAnalyzerExample: OctetAnalyzerExample.cc libUCNA.a
	$(CXX) $(CXXFLAGS) Studies/OctetAnalyzerExample.cc -o OctetAnalyzerExample $(LDFLAGS)
	
ExtractFierzTerm: ExtractFierzTerm.cc libUCNA.a
	$(CXX) $(CXXFLAGS) ExtractFierzTerm.cc -o ExtractFierzTerm $(LDFLAGS)

ExtractCorrectBetaSpectrum: ExtractCorrectBetaSpectrum.cc libUCNA.a
	$(CXX) $(CXXFLAGS) ExtractCorrectBetaSpectrum.cc -o ExtractCorrectBetaSpectrum $(LDFLAGS)

MWPC_Efficiency_Sim: MWPC_Efficiency_Sim.cc libUCNA.a
	$(CXX) $(CXXFLAGS) Studies/MWPC_Efficiency_Sim.cc -o MWPC_Efficiency_Sim $(LDFLAGS)
	
QCalc: FPNCalc.cc libUCNA.a
	$(CXX) $(CXXFLAGS) IOUtils/FPNCalc.cc -o QCalc $(LDFLAGS)
	
FierzOctetAnalyzer: FierzOctetAnalyzer.cc libUCNA.a
	$(CXX) $(CXXFLAGS) Studies/FierzOctetAnalyzer.cc -o FierzOctetAnalyzer $(LDFLAGS)
	
MC_Comparisons: Studies/MC_Comparisons.cc libUCNA.a
	$(CXX) $(CXXFLAGS) Studies/MC_Comparisons.cc -o MC_Comparisons $(LDFLAGS)

#
# documentation via Doxygen
#

doc : latex/refman.pdf

latex/refman.pdf: latex/ 
	cd latex; make
latex/ : Doxyfile
	doxygen

#
# cleanup
#
.PHONY: clean
clean:
	-rm -f libUCNA.a UCNAnalyzer OctetAnalyzerExample DataScannerExample CalibratorExample
	-rm -f ExtractFierzTerm Analyzer MWPC_Efficiency_Sim
	-rm -f *.o
	-rm -rf *.dSYM
	-rm -rf latex/
	-rm -rf html/
	
