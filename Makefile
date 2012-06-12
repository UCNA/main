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

CXXFLAGS = -O3 -Wall `root-config --cflags` \
	-I. -IIOUtils -IRootUtils -IBaseTypes -IMathUtils -ICalibration -IAnalysis -IStudies -IPhysics
LDFLAGS = `root-config --libs` -lSpectrum 

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
	PenelopeToPMT.o TH1toPMT.o KurieFitter.o EndpointStudy.o ReSource.o EfficCurve.o 

Studies = PlotMakers.o SegmentSaver.o RunAccumulator.o OctetAnalyzer.o AsymmetryAnalyzer.o WirechamberStudy.o LEDScans.o

objects = $(Utils) $(Calibration) $(Analysis) $(Studies) $(Physics)

all:
	make UCNAnalyzer

UCNAnalyzer: Analyzer.cpp $(objects)
	$(CXX) $(CXXFLAGS) Analyzer.cpp $(objects) -o UCNAnalyzer $(LDFLAGS)

CalibratorExample: CalibratorExample.cpp $(objects)
	$(CXX) $(CXXFLAGS) $(LDFLAGS)  $(objects) -o CalibratorExample

DataScannerExample: DataScannerExample.cc $(objects)
	$(CXX) $(CXXFLAGS) DataScannerExample.cc $(objects) -o DataScannerExample $(LDFLAGS)

OctetAnalyzerExample: OctetAnalyzerExample.cc $(objects)
	$(CXX) $(CXXFLAGS) Studies/OctetAnalyzerExample.cc $(objects) -o OctetAnalyzerExample $(LDFLAGS)
	
ExtractFierzTerm: ExtractFierzTerm.cc $(objects)
	$(CXX) $(CXXFLAGS) ExtractFierzTerm.cc $(objects) -o ExtractFierzTerm $(LDFLAGS)

ExtractCorrectBetaSpectrum: ExtractCorrectBetaSpectrum.cc $(objects)
	$(CXX) $(CXXFLAGS) ExtractCorrectBetaSpectrum.cc $(objects) -o ExtractCorrectBetaSpectrum $(LDFLAGS)

MWPC_Efficiency_Sim: MWPC_Efficiency_Sim.cc $(objects)
	$(CXX) $(CXXFLAGS) Studies/MWPC_Efficiency_Sim.cc $(objects) -o MWPC_Efficiency_Sim $(LDFLAGS)
	
QCalc: FPNCalc.cc $(Utils)
	$(CXX) $(CXXFLAGS) IOUtils/FPNCalc.cc $(objects) -o QCalc $(LDFLAGS)
	
FierzOctetAnalyzer: FierzOctetAnalyzer.cc $(objects)
	$(CXX) $(CXXFLAGS) Studies/FierzOctetAnalyzer.cc $(objects) -o FierzOctetAnalyzer $(LDFLAGS)
	
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
	-rm -f UCNAnalyzer OctetAnalyzerExample DataScannerExample CalibratorExample ExtractFierzTerm Analyzer MWPC_Efficiency_Sim
	-rm -f *.o
	-rm -rf *.dSYM
	-rm -rf latex/
	-rm -rf html/
	
