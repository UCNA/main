#!/bin/bash



#
# CUORE base code svn and git repository
#
if [ -z "${CUOREBASE}" ]
then
  export CUOREBASE=$(pwd)
fi

#while true; do
#  Select "Install in ${CUOREBASE}? (" yn
#  case $yn in
#  	[Yy]* ) break;;
#	[Nn]*

read -p "What is your svn username? " username
export CUOREUSER=$username

export CUORESVN=cuore-svn.roma1.infn.it/cuoremc/trunk
export CUOREGIT=git@github.com:hickerson/cuore.git
export CUOREPATCH=${CUOREBASE}/patch
export CUORESRC=${CUOREBASE}/sources
export CUOREBIN=${CUOREBASE}/bin
export CUOREBUILD=${CUOREBASE}/build
export PATH=${PATH}:${CUOREBIN}
export CUOREEXE="qshields g2tas g2root MCuoreZ"



#
# Geant4 download, system, version and release information
#
if [ -z "${G4SRC}" ]
then
	# defaults to putting geant4 in the current directory
	export G4SRC=$(pwd)/geant4
fi

export G4URL=http://geant4.cern.ch/support/source
#export G4VERSION=9.6
export G4VERSION=10
export G4RELEASE=01
export G4SYSTEM=Linux
export G4CC=g++
export G4CCDIR=${G4SYSTEM}-${G4CC}
export G4NAME=geant4.${G4VERSION}.${G4RELEASE}
export G4NPROC=$(nproc)



#
# Basic Geant4 environment variables
#
export G4VDIR=${G4SRC}/v${G4VERSION}
export G4PATH=${G4VDIR}/${G4NAME}
export G4PACKAGES=${G4SRC}/packages
export G4INSTALL=${G4PATH}
#export G4INSTALL=${G4PATH}-install
export G4BUILD=${G4PATH}-build
#export Geant4_USE_FILE=${G4BUILD}/Geant4Config.cmake
#export Geant4_DIR=${G4BUILD}/Geant4Config.cmake
export G4WORKDIR=${G4SRC}/scratch
export G4TMP=${G4WORKDIR}/tmp
export G4BIN=${G4SRC}/bin/v${G4VERSION}
export G4LIB=${G4BUILD}/outputs/library
export G4INCLUDE=${G4INSTALL}/include
export GSPACEDIR=${G4LIB}/${G4CCDIR}
export G4DATA=${G4BUILD}/data



#
# CLHEP environment variables
#
#export CLHEP_VERSION=2.1.3.1
#export CLHEP_BASE_DIR=/usr/local
if [ -z "${CLHEP_VERSION}" ]
then
	echo "Using default CLHEP included with ${G4NAME}"
	export CLHEP_BASE_DIR=${G4PATH}/source/externals/clhep
	export CLHEP_INCLUDE_DIR=${CLHEP_BASE_DIR}/include
	export LD_LIBRARY_PATH=${CLHEP_BASE_DIR}/lib:${G4LIB}/${G4CCDIR}:${LD_LIBRARY_PATH}
else
	echo "Manually installing CLHEP"
	#probably should put more tests here about wether the dir is there
	export CLHEP_BASE_DIR=${G4SRC}/CLHEP-${CLHEP_VERSION}
	export CLHEP_INCLUDE_DIR=${CLHEP_BASE_DIR}/include
	#export LD_LIBRARY_PATH=${CLHEP_BASE_DIR}/lib:${G4LIB}/Linux-g++:${LD_LIBRARY_PATH}
	export LD_LIBRARY_PATH=${CLHEP_BASE_DIR}/lib:${LD_LIBRARY_PATH}
fi



# Export path for vrmlviewer (the bin directory)
#export PATH=${PATH}:${CUOREBASE}/bin/vrmlviewer



#export HEP_ODBMS_DIR=${G4INSTALL}/CLHEP-2.1.3.1
#export HEP_ODBMS_DIR=${CLHEP_INCLUDE_DIR}



#
# Visualization environment variables
#
export OGLHOME=""
export OGLFLAGS=""
export OGLLIBS=""
export G4USE_STL=1
export G4VIS_BUILD_OGLSX_DRIVER=1
export G4VIS_USE_OGLSX=1
export G4VIS_BUILD_OGLIX_DRIVER=''
export G4VIS_USE_OGLIX=''
export G4VIS_BUILD_OPENGLX_DRIVER=''
export G4VIS_USE_OPENGLX=''
export G4VIS_BUILD_OPENGLXM_DRIVER=''
export G4VIS_USE_OPENGLXM=''
export G4VIS_BUILD_VRMLFILE_DRIVER=1
export G4VIS_BUILD_VRML_DRIVER=1
export G4VIS_USE_VRMLFILE=1
export G4VIS_USE_VRML=1
export G4VIS_BUILD_DAWNFILE_DRIVER=1
export G4VIS_BUILD_DAWN_DRIVER=1
export G4VIS_USE_DAWNFILE=1
export G4VIS_USE_DAWN=1
export G4VRMLFILE_VIEWER='vrmlview'
export G4DAWNFILE_VIEWER=dawn
export DAWN_PS_PREVIEWER=display



#
# Data environment variables
#
if [ -z "${G4DATA}" ]
then
	echo "Using data downloaded with ${G4NAME}"
else
	echo "Using data files in ${G4DATA}"
	# We should probably test that G4DATA directory exists and contains the data dirs
	export G4LEDATA=${G4DATA}/G4EMLOW6.32
	export G4LEVELGAMMADATA=${G4DATA}/PhotonEvaporation2.3
	export G4NEUTRONHPDATA=${G4DATA}/G4NDL4.2
	export G4RADIOACTIVEDATA=${G4DATA}/RadioactiveDecay3.6
	export G4ABLADATA=${G4DATA}/G4ABLA3.0
	export G4REALSURFACEDATA=${G4DATA}/RealSurface1.0
	export G4NEUTRONXSDATA=${G4DATA}/G4NEUTRONXS1.2
	export G4PIIDATA=${G4DATA}/G4PII1.3 
	#export G4SAIDDATA=${G4DATA}/G4SAIDDATA1.1
	export G4SAIDXSDATA=${G4DATA}/G4SAIDDATA1.1
	#export GENDATA=${G4DATA}/GENDATA
	export G4MUDATA=${G4DATA}/MUDATA
	export NeutronHPCrossSections=${G4DATA}/G4NDL4.2
fi

#
# ROOT path to Cmake files
#
export G4ROOT=${G4INSTALL}/cmake/Modules


#
# Library exports for building Geant4 code
#
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${G4INSTALL}/lib:${G4LIB}/${G4CCDIR}
export PATH=${G4BIN}:${PATH}
