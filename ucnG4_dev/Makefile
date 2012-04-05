# --------------------------------------------------------------
# Makefile for UCNA Geant4 MC
# M. Mendenhall, 3/2012
# --------------------------------------------------------------

# settings for Geant4 make scripts
name := ucnG4_prod
G4TARGET := $(name)
G4EXLIB := true
ifndef G4INSTALL
	G4INSTALL = ../../..
endif
EXTRALIBS += `root-config --libs --glibs` -L MCEvent

# compiler flags
CXX = g++
CPPFLAGS = -O3 -Wall `root-config --cflags`

.PHONY: all prereq

all: 
	$(MAKE) dict
	$(MAKE) libMCEvent.so
	$(MAKE) lib  bin
	$(MAKE) UCNA_MC_Analyzer
	$(MAKE) SiDet_Analyzer
	
dict: include/bmMCEventDict.h include/bmTrackInfoDict.h \
	include/bmPrimaryInfoDict.h

include/bmMCEventDict.h: include/bmMCEvent.hh
	cd include && \
        rootcint -f bmMCEventDict.cc -c bmMCEvent.hh && \
	mv bmMCEventDict.cc ../src

include/bmTrackInfoDict.h: include/bmTrackInfo.hh
	cd include && \
        rootcint -f bmTrackInfoDict.cc -c bmTrackInfo.hh && \
	mv bmTrackInfoDict.cc ../src

include/bmPrimaryInfoDict.h: include/bmPrimaryInfo.hh
	cd include && \
        rootcint -f bmPrimaryInfoDict.cc -c bmPrimaryInfo.hh && \
	mv bmPrimaryInfoDict.cc ../src

# make a compact event library for ROOT clients
libMCEvent.so: dict
	rm -rf tmplibbuild &&\
	mkdir tmplibbuild && \
	cp include/bmTrackInfo.hh tmplibbuild &&\
	cp include/bmTrackInfoDict.h tmplibbuild &&\
	cp src/bmTrackInfo.cc tmplibbuild &&\
	cp src/bmTrackInfoDict.cc tmplibbuild &&\
	cp include/bmPrimaryInfo.hh tmplibbuild &&\
	cp include/bmPrimaryInfoDict.h tmplibbuild &&\
	cp src/bmPrimaryInfo.cc tmplibbuild &&\
	cp src/bmPrimaryInfoDict.cc tmplibbuild &&\
	cp include/bmMCEvent.hh tmplibbuild && \
	cp include/bmMCEventDict.h tmplibbuild && \
	cp src/bmMCEvent.cc tmplibbuild && \
	cp src/bmMCEventDict.cc tmplibbuild && \
	cp makeEvtLib.gmk tmplibbuild &&\
	cd tmplibbuild; ${MAKE} -f makeEvtLib.gmk &&\
	cd ..; mv -f tmplibbuild/libMCEvent.so .; rm -rf tmplibbuild

# Analyzer for UCNA simulation data
UCNA_MC_Analyzer: UCNA_MC_Analyzer.cc AnalyzerBase.o libMCEvent.so
	$(CXX) $(CPPFLAGS) UCNA_MC_Analyzer.cc $(EXTRALIBS) AnalyzerBase.o ./libMCEvent.so -o UCNA_MC_Analyzer

# Analyzer for simple silicon detector data
SiDet_Analyzer: SiDet_Analyzer.cc AnalyzerBase.o libMCEvent.so
	$(CXX) $(CPPFLAGS) SiDet_Analyzer.cc $(EXTRALIBS) AnalyzerBase.o ./libMCEvent.so -o SiDet_Analyzer


# cleanup
clean::
	rm -rf libMCEvent.so include/bmMCEventDict.h src/bmMCEventDict.cc \
		include/bmTrackInfoDict.h src/bmTrackInfoDict.cc \
		include/bmPrimaryInfoDict.h src/bmPrimaryInfoDict.cc
	rm -rf tmplibbuild
	rm -rf $(G4WORKDIR)/tmp/
	rm -rf $(G4WORKDIR)/bin/
	rm -f  UCNA_MC_Analyzer SiDet_Analyzer

include $(G4INSTALL)/config/binmake.gmk 
