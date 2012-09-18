.SUFFIXES:
.SUFFIXES: .c .cc .cpp .o

ROOTGLIBS     = $(shell root-config --glibs)
ROOTLIBS      = $(shell root-config --libs)
ROOTCFLAGS    = $(shell root-config --cflags)

CXX           = g++ 
INCLUDEPATH   = -I./ -I../IOUtils -I../RootUtils -I../BaseTypes -I../Detectors -I../MathUtils -I../Calibration -I../Analysis -I../Studies
CXXFLAGS      = -O3 -Wall -fPIC $(ROOTCFLAGS) $(INCLUDEPATH)
LDFLAGS       = -L.. -lUCNA $(ROOTLIBS) $(ROOTGLIBS) -lSpectrum -lstdc++ -lz

ANA_OBJS      = ManualInfo.o ucnaDataAnalyzer_hists.o ucnaDataAnalyzer_tree.o ucnaDataAnalyzer_peds.o

all: ucnaDataAnalyzer11b

../libUCNA.a:
	cd ..; make libUCNA.a

ucnaDataAnalyzer11b: ucnaDataAnalyzer11b.cc ucnaDataAnalyzer11b.hh $(ANA_OBJS) ../libUCNA.a
	$(CXX) $(CXXFLAGS) ucnaDataAnalyzer11b.cc $(ANA_OBJS) $(LDFLAGS) -o ucnaDataAnalyzer11b

clean:
	@rm -f *.o *~
	@rm -f ucnaDataAnalyzer11b
