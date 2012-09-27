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

CXXFLAGS = -O3 -Wall -fPIC `root-config --cflags` -I. \
	-IIOUtils -IRootUtils -IBaseTypes -IMathUtils -ICalibration -IAnalysis -IStudies -IPhysics
LDFLAGS =  -L. -lUCNA -lSpectrum `root-config --libs`

ifdef PROFILER_COMPILE
	CXXFLAGS += -pg
	LDFLAGS += -pg
endif

ifdef UNBLINDED
	CXXFLAGS += -DUNBLINDED
endif

#
# things to build
#

VPATH = ./:IOUtils/:RootUtils/:BaseTypes/:MathUtils/:Calibration/:Analysis/:Studies/:Physics/:Examples/

Physics = BetaSpectrum.o ElectronBindingEnergy.o NuclEvtGen.o

Utils = ControlMenu.o strutils.o PathUtils.o TSpectrumUtils.o QFile.o GraphUtils.o MultiGaus.o TagCounter.o \
		SectorCutter.o LinHistCombo.o Enums.o Types.o FloatErr.o SMExcept.o Octet.o SpectrumPeak.o Source.o \
		SQL_Utils.o GraphicsUtils.o OutputManager.o RollingWindow.o RData.o EnumerationFitter.o

Calibration = PositionResponse.o PMTGenerator.o CathSegCalibrator.o WirechamberCalibrator.o \
		EnergyCalibrator.o CalDBSQL.o SourceDBSQL.o GainStabilizer.o EvisConverter.o
	
Analysis = TChainScanner.o RunSetScanner.o ProcessedDataScanner.o PostOfficialAnalyzer.o Sim2PMT.o G4toPMT.o \
		PenelopeToPMT.o TH1toPMT.o KurieFitter.o ReSource.o EfficCurve.o AnalysisDB.o

Studies = SegmentSaver.o RunAccumulator.o OctetAnalyzer.o \
	MuonAnalyzer.o PositionAnalyzer.o WirechamberGainAnalyzer.o BGDecayAnalyzer.o HighEnergyExcess.o \
	AsymmetryAnalyzer.o SimAsymmetryAnalyzer.o BetaDecayAnalyzer.o \
	CathodeTweakAnalyzer.o PositionBinnedAnalyzer.o AnodePositionAnalyzer.o XenonAnalyzer.o \
	PlotMakers.o LEDScans.o AsymmetryCorrections.o

objects = $(Utils) $(Calibration) $(Analysis) $(Studies) $(Physics)

ExampleObjs = CalibratorExample DataScannerExample ExtractFierzTerm \
	QCalc MC_Comparisons MWPC_Efficiency_Sim FierzOctetAnalyzer OctetAnalyzerExample

all: UCNAnalyzer
	
libUCNA.a: $(objects)
	ar rs libUCNA.a $(objects)
	
UCNAnalyzer: Analyzer.cpp libUCNA.a
	$(CXX) $(CXXFLAGS) Analyzer.cpp $(LDFLAGS) -o UCNAnalyzer

ExtractFierzTerm: ExtractFierzTerm.cc libUCNA.a
	$(CXX) $(CXXFLAGS) ExtractFierzTerm.cc $(LDFLAGS) -o ExtractFierzTerm
	
examples: $(ExampleObjs)

CalibratorExample: CalibratorExample.cpp libUCNA.a
	$(CXX) $(CXXFLAGS) Examples/CalibratorExample.cpp $(LDFLAGS) -o CalibratorExample

DataScannerExample: DataScannerExample.cc libUCNA.a
	$(CXX) $(CXXFLAGS) Examples/DataScannerExample.cc $(LDFLAGS) -o DataScannerExample

OctetAnalyzerExample: OctetAnalyzerExample.cc libUCNA.a
	$(CXX) $(CXXFLAGS) Examples/OctetAnalyzerExample.cc $(LDFLAGS) -o OctetAnalyzerExample

ExtractCorrectBetaSpectrum: ExtractCorrectBetaSpectrum.cc libUCNA.a
	$(CXX) $(CXXFLAGS) Examples/ExtractCorrectBetaSpectrum.cc $(LDFLAGS) -o ExtractCorrectBetaSpectrum

MWPC_Efficiency_Sim: MWPC_Efficiency_Sim.cc libUCNA.a
	$(CXX) $(CXXFLAGS) Examples/MWPC_Efficiency_Sim.cc $(LDFLAGS) -o MWPC_Efficiency_Sim
	
QCalc: FPNCalc.cc libUCNA.a
	$(CXX) $(CXXFLAGS) Examples/FPNCalc.cc $(LDFLAGS) -o QCalc
	
FierzOctetAnalyzer: FierzOctetAnalyzer.cc libUCNA.a
	$(CXX) $(CXXFLAGS) Examples/FierzOctetAnalyzer.cc $(LDFLAGS) -o FierzOctetAnalyzer
	
MC_Comparisons: MC_Comparisons.cc libUCNA.a
	$(CXX) $(CXXFLAGS) Examples/MC_Comparisons.cc $(LDFLAGS) -o MC_Comparisons

ucnG4:
	mkdir -p g4build/
	cd g4build; cmake -DGeant4_DIR=~/geant4.9.5/geant4.9.5-install/lib/Geant4-9.5.0/ ../ucnG4_dev/; make

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
	-rm -f libUCNA.a UCNAnalyzer $(ExampleObjs)
	-rm -f *.o
	-rm -rf *.dSYM
	-rm -rf latex/
	-rm -rf html/
	
