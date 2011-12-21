# $Id: GNUmakefile,v 1.8 2011-10-12 19:31:15 mmendenhall Exp $
# --------------------------------------------------------------
# GNUmakefile for UCNA Geant4 MC
# --------------------------------------------------------------

name := ucnG4_prod
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

CPPFLAGS  += -O3 `root-config --cflags` #-m32
EXTRALIBS += `root-config --libs --glibs`
EXTRALIBS += "-L MCEvent"

.PHONY: all prereq

all: 
	$(MAKE) dict
	$(MAKE) libMCEvent.so
	$(MAKE) lib bin 

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

clean::
	rm -f libMCEvent.so include/bmMCEventDict.h src/bmMCEventDict.cc \
		include/bmTrackInfoDict.h src/bmTrackInfoDict.cc \
		include/bmPrimaryInfoDict.h src/bmPrimaryInfoDict.cc
	rm -rf $(G4WORKDIR)/tmp/
	rm -f *~
	rm -f src/*~
	rm -f include/*~
	rm -rf ../tmp/
	rm -f analyzer

include $(G4INSTALL)/config/binmake.gmk 
