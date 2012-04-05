ObjSuf        = o
SrcSuf        = cc
DllSuf        = so
OutPutOpt     = -o  

ROOTGLIBS     = $(shell root-config --glibs)
ROOTLIBS      = $(shell root-config --libs)
ROOTCFLAGS    = $(shell root-config --cflags)

VPATH = ./:../IOUtils/:../RootUtils/:../BaseTypes/:../Detectors/:../MathUtils/:../Calibration/:../Analysis/:Studies/

CXX           = g++ 
CXXFLAGS      = -O3 -m32 -Wall -fPIC
LD            = g++
SOFLAGS       = -shared
INCLUDEPATH   = -I./ -I../IOUtils -I../RootUtils -I../BaseTypes -I../Detectors -I../MathUtils -I../Calibration -I../Analysis -I../Studies
CXXFLAGS     += $(ROOTCFLAGS) $(INCLUDEPATH)

LIBS          = $(ROOTLIBS) $(ROOTGLIBS) -lSpectrum
LIBS         += -lstdc++ -lz

RECON_OBJS   = EnergyCalibrator.o CalDBSQL.o Source.o SourceDBSQL.o GainStabilizer.o SQL_Utils.o PositionResponse.o EfficCurve.o EvisConverter.o \
				WirechamberReconstruction.o strutils.o QFile.o RData.o OutputManager.o PathUtils.o Enums.o Types.o GraphUtils.o \
				TChainScanner.o ManualInfo.o GraphicsUtils.o RollingWindow.o WirechamberCalibrator.o MultiGaus.o SpectrumPeak.o

ANA_OBJS	= ucnaDataAnalyzer_hists.o ucnaDataAnalyzer_tree.o ucnaDataAnalyzer_peds.o

all: ucnaDataAnalyzer11b

ucnaDataAnalyzer11b: $(RECON_OBJS) $(ANA_OBJS) ucnaDataAnalyzer11b.cc ucnaDataAnalyzer11b.hh
	$(CXX) $(CXXFLAGS) ucnaDataAnalyzer11b.cc $(RECON_OBJS) $(ANA_OBJS) $(LIBS) -o ucnaDataAnalyzer11b

######################################################################
#  Generic compilation and linking step to make an executable from
#  a single *.cc file
#
%: %.$(SrcSuf)
	$(CXX) $(CXXFLAGS) -o $@ $< $(LIBS)
	@echo "$@ done"

clean:
	@rm -f *.o *~
	@rm -f *.so
	@rm -f ucnaDataAnalyzer11b
