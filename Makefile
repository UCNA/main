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
LDFLAGS =  -L. -lUCNA -lSpectrum -lMLP `root-config --libs` -lMathMore

ifdef PROFILER_COMPILE
	CXXFLAGS += -pg
	LDFLAGS += -pg
endif

ifdef UNBLINDED
	CXXFLAGS += -DUNBLINDED
endif

ifdef PUBLICATION_PLOTS
	CXXFLAGS += -DPUBLICATION_PLOTS
endif

ifdef TSPECTRUM_USES_DOUBLE
	CXXFLAGS += -DTSPECTRUM_USES_DOUBLE
endif

#
# things to build
#

VPATH = ./:IOUtils/:RootUtils/:BaseTypes/:MathUtils/:Calibration/:Analysis/:Studies/:Physics/:Examples/:Fierz/:Standalone/

Physics = BetaSpectrum.o ElectronBindingEnergy.o NuclEvtGen.o

Utils = ControlMenu.o strutils.o ManualInfo.o PathUtils.o TSpectrumUtils.o QFile.o GraphUtils.o MultiGaus.o TagCounter.o \
		SectorCutter.o LinHistCombo.o Enums.o Types.o FloatErr.o SMExcept.o Octet.o SpectrumPeak.o Source.o \
		SQL_Utils.o GraphicsUtils.o OutputManager.o RollingWindow.o EnumerationFitter.o StyleSetup.o TChainScanner.o

Calibration = PositionResponse.o PMTGenerator.o CathSegCalibrator.o WirechamberCalibrator.o \
		EnergyCalibrator.o CalDBSQL.o SourceDBSQL.o GainStabilizer.o EvisConverter.o EventClassifier.o
	
Analysis = RunSetScanner.o ProcessedDataScanner.o PostOfficialAnalyzer.o Sim2PMT.o G4toPMT.o \
		PenelopeToPMT.o LED2PMT.o TH1toPMT.o KurieFitter.o ReSource.o EfficCurve.o AnalysisDB.o

Studies = SegmentSaver.o RunAccumulator.o OctetAnalyzer.o OctetSimuCloneManager.o \
	MuonPlugin.o PositionsPlugin.o WirechamberEnergyPlugins.o BGDecayPlugin.o HighEnergyExcessPlugin.o \
	AsymmetryPlugin.o SimAsymmetryPlugin.o BetaDecayAnalyzer.o \
	CathodeTuningAnalyzer.o PositionBinnedPlugin.o WirechamberGainMapPlugins.o XenonAnalyzer.o \
	PlotMakers.o AsymmetryCorrections.o FierzFitter.o GravitySpectrometerPlugin.o SimEdepPlugin.o

objects = $(Utils) $(Calibration) $(Analysis) $(Studies) $(Physics)




all: UCNAnalyzer

UCNAnalyzer: Analyzer.cpp libUCNA.a
	$(CXX) $(CXXFLAGS) Analyzer.cpp $(LDFLAGS) -o UCNAnalyzer

libUCNA.a: $(objects)
	ar rs libUCNA.a $(objects)

# generic rule for everything else .cc linked against libUCNA
% : %.cc libUCNA.a
	$(CXX) $(CXXFLAGS) $< $(LDFLAGS) -o $@
	

ExampleObjs = CalibratorExample DataScannerExample ExtractFierzTerm CombinedAbFit \
	FPNCalc  MWPC_Efficiency_Sim FierzOctetAnalyzer OctetAnalyzerExample

examples: $(ExampleObjs)

StandaloneObjs = GammaComptons BetaEndpoint BetaOctetPositions MC_Comparisons LEDScans MiscJunk \
					MC_EventGen QuasiRandomTest MWPC_Energy_Cal

standalone: $(StandaloneObjs)




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
	-rm -f libUCNA.a UCNAnalyzer
	-rm -f $(ExampleObjs)
	-rm -f $(StandaloneObjs)
	-rm -f *.o
	-rm -rf *.dSYM
	-rm -rf latex/
	-rm -rf html/
	
