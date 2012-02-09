#####################################################################
#
#  Name:         Makefile
#  Created by:   Stefan Ritt
#
#  Contents:     Makefile for MIDAS example frontend and analyzer
#
#  $Id: Makefile 3203 2006-07-31 21:39:02Z ritt $
#
#####################################################################
#
#--------------------------------------------------------------------
# The MIDASSYS should be defined prior the use of this Makefile
ifndef MIDASSYS
missmidas::
	@echo "...";
	@echo "Missing definition of environment variable 'MIDASSYS' !";
	@echo "...";
endif

# get OS type from shell
OSTYPE = $(shell uname)

#--------------------------------------------------------------------
# The following lines contain specific switches for different UNIX
# systems. Find the one which matches your OS and outcomment the 
# lines below.

ifndef CERNLIB
misscern::
	@echo "..."
	@echo "Missing definition of the environment variable 'CERNLIB' "
	@echo "..."
else
CERNLIB_PACK  = $(CERNLIB)/libpacklib.a
endif

#-----------------------------------------
# This is for Linux
ifeq ($(OSTYPE),Linux)
OSTYPE = linux
endif

ifeq ($(OSTYPE),linux)

OS_DIR = linux
OSFLAGS = -DOS_LINUX -Dextname
CFLAGS = -g -O2 -Wall
# add to compile in 32-bit mode
# OSFLAGS += -m32
CFLAGS += -D_LARGEFILE64_SOURCE
FF = g77
LIBS = -lm -lz -lutil -lnsl -lpthread
endif

#-----------------------
# MacOSX/Darwin is just a funny Linux
#
ifeq ($(OSTYPE),Darwin)
OSTYPE = darwin
endif

ifeq ($(OSTYPE),darwin)
OS_DIR = darwin
FF = cc
OSFLAGS = -DOS_LINUX -DOS_DARWIN -DHAVE_STRLCPY -DAbsoftUNIXFortran -fPIC -Wno-unused-function
LIBS = -lpthread
SPECIFIC_OS_PRG = $(BIN_DIR)/mlxspeaker
NEED_STRLCPY=
NEED_RANLIB=1
NEED_SHLIB=
NEED_RPATH=

endif

#-------------------------------------------------------------------
# The following lines define directories. Adjust if necessary
#                 
DRV_DIR   = $(MIDASSYS)/drivers/vme/sis3100/linux
INC_DIR   = $(MIDASSYS)/include 
LIB_DIR   = $(MIDASSYS)/lib
SRC_DIR   = $(MIDASSYS)/src
CFOR_INC_DIR    = $(CERN_ROOT)/include/cfortran

#-------------------------------------------------------------------
# List of analyzer modules
#
MODULES   = adccalib.o adcsum.o scaler.o tdccalib.o

#-------------------------------------------------------------------
# Hardware driver can be (camacnul, kcs2926, kcs2927, hyt1331)
#
DRIVER = sis3100_vme_calls

#-------------------------------------------------------------------
# Frontend code name defaulted to frontend in this example.
# comment out the line and run your own frontend as follow:
# gmake UFE=my_frontend
#
UFE = frontend
##### uncomment line above for normal running######
##### LEACOCK added this line to test changes to frontend#######
#UFE = frontend_leacock_test


####################################################################
# Lines below here should not be edited
####################################################################

# MIDAS library
LIB = $(LIB_DIR)/libmidas.a

# compiler
CC = gcc
CXX = g++
CFLAGS += -g -I${CFOR_INC_DIR} -I$(INC_DIR) -I$(DRV_DIR)
LDFLAGS +=

all: $(UFE) analyzer

#$(UFE): $(LIB) $(LIB_DIR)/mfe.o $(DRIVER).o $(UFE).c $(SRC_DIR)/cnaf_callback.c
$(UFE): $(LIB) $(LIB_DIR)/mfe.o $(DRIVER).o $(UFE).c
	$(CC) $(CFLAGS) $(OSFLAGS) -o $(UFE) $(UFE).c \
	$(DRIVER).o  $(LIB_DIR)/mfe.o $(LIB) \
	$(LDFEFLAGS) $(LIBS)

$(DRIVER).o: $(DRV_DIR)/$(DRIVER).c
	$(CC) $(CFLAGS) $(OSFLAGS) -c -o $(DRIVER).o $(DRV_DIR)/$(DRIVER).c

analyzer: $(LIB) $(LIB_DIR)/hmana.o analyzer.o $(MODULES)
	$(FF) $(FFLAGS) $(CFLAGS) -o analyzer $(LIB_DIR)/hmana.o analyzer.o \
	$(MODULES) $(CERNLIB_PACK) \
	$(LIB) $(LDFLAGS) $(LIBS)

%.o: %.c experim.h
	$(CC) $(USERFLAGS) $(CFLAGS) $(OSFLAGS) -o $@ -c $<

clean::
	rm -f *.o *~ \#*

#end file
