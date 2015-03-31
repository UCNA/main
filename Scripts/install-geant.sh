#!/bin/sh

##
#
# Adapted into a sh script by Kevin Peter Hickerson
# From CUORE MC group install instructions
#
# Aug 20, 2013
#
##


#
# Environment preparation
#
echo "Checking that we have everything we need..."
#for CMD in svn cmake make wget tar nproc
for CMD in svn cmake make tar root
do
	command -v $CMD >/dev/null 2>&1 || { echo "Require $CMD but it's not installed. Aborting." >&2; exit 1; }
done

echo "Checking if the environment is setup..."
if [ -z "${G4SRC}" ]; then
	echo "Need to set G4SRC."
	echo "Did set your genat4 install path?"
	echo "Set it, maybe in your ~/.bashrc file, then rerun your geant4 export script."
	echo "Aborting install."
	exit 1
fi

if [ -z "${G4PATH}" ]; then
	echo "Need to set G4PATH."
	echo "Export script seems to not be working."
	echo "Aborting install."
	exit 1
fi

if [ -z "${G4NAME}" ]; then
	echo "Need to set G4NAME"
	echo "Did you run your export script?"
	echo "Aborting install."
	exit 1
fi



echo "Looks like we have everything we need."

echo "CUORE svn repository: ${CUORESVN}"
echo "CUORE git repository: ${CUOREGIT}"
echo "Geant4 source url: ${G4URL}"
echo "Geant4 version: ${G4NAME}"
echo "System: ${G4SYSTEM}"
echo "Install: ${G4SRC}"
echo "Compiler: ${G4CC}"
echo "Cores: ${G4NPROC}"


echo "About to install CUORE Geant4.${G4VERSION}.${G4RELEASE} software in ${G4SRC}"
echo "Creating install directories..."
mkdir -vp ${G4SRC}
mkdir -vp ${G4PATH}
#mkdir -vp ${G4DATA}
mkdir -vp ${G4PACKAGES}
#mkdir ${G4SRC}/CLHEP
#mkdir ${G4SRC}/CLHEP-2.1.2.3
#mkdir ${G4SRC}/packages/G4.9.6

#cd ${G4SRC}/v${G4RELEASE}

# if the installation is done in private location/PC $CUORE_BASE_DIR must be substituted with /localdir/ in geant4.9.6-exports.sh and in all the following steps in this procedure
# Check RELEASE definition in geant4.9.6-exports.sh


#
# Software download files from 
# http://geant4.cern.ch/support/download.shtml 
# or 
# ${G4URL}
#
echo "Downloading packages from ${G4URL}..."
cd ${G4PACKAGES}
if [ -f ${G4NAME}.tar.gz ];
then
	echo "${G4NAME}.tar.gz already exists. Replacing file..."
	rm ${G4NAME}.tar.gz
fi

case ${G4SYSTEM} in
	"Linux") wget --tries=inf --retry-connrefused ${G4URL}/${G4NAME}.tar.gz;;
	"OSX") curl "${G4URL}/${G4NAME}.tar.gz" -o "${G4NAME}.tar.gz";;
	*) echo "Sorry, don't know what to do with system ${G4SYSTEM}."; exit 1;;
esac

# don't need these anymore: all data files .tgz
# don't need these anymore: clhep-2.1.3.1.tgz
# both of these are included as of 9.6


#
# Data files
#
# Downloading data files are included in the geant install process
#cp ${G4PACKAGES}/${G4NAME}.tar.gz ${G4SRC}/data
#cd  ${G4SRC}/data
#for i in *; 
#	do tar -xzf  $i;
#done
#rm -f *.gz
#check if G4.9.6exports data-files definitions is compatible with the downloaded ones and modify it in the case some files are changed


#
# CLHEP installation
#
# CLHEP is now included with geant4
#tar zxvf clhep-2.1.2.3.tgz
#mv 2.1.2.3/CLHEP ${G4SRC}/v9.6/CLHEP/source
#cd  ${G4SRC}/v9.6/CLHEP
#mkdir build
#cd build
#cmake -DCMAKE_INSTALL_PREFIX=${G4SRC}/v9.6/CLHEP-2.1.2.3 ${G4SRC}/v9.6/CLHEP/source
#crio: cmake -DCMAKE_INSTALL_PREFIX=/usr/remote/geant4/v4.9.6/CLHEP-2.1.2.3 /usr/remote/geant4/v4.9.6/CLHEP/source
#make
#sudo make install
#cd ../..  #(go in ${G4SRC}/v9.6)


#
# Base GEANT4 installation
#
echo "Unpacking files downloaded from ${G4URL}..."
tar zxf ${G4PACKAGES}/${G4NAME}.tar.gz
if [ -d ${G4VDIR}/${G4NAME} ];
then
	rm -rf ${G4VDIR}/${G4NAME}
fi
mv ${G4PACKAGES}/${G4NAME} ${G4VDIR}/


#
# Patch all of geant4 with ${G4NAME}.diff
#
cp ${CUOREPATCH}/${G4NAME}.diff ${G4VDIR}
cd ${G4VDIR}
patch -p0 <${G4NAME}.diff
rm ${G4VDIR}/${G4NAME}.diff


#find source -name include -exec cp -r {}  .  \;
#edit cmake/Templates/geant4-config.in: add -lG4esources
# patch...
#edit cmake/Templates/Geant4Config.cmake.in: add G4esources${_geant4_lib_use_suffix} to 'kernel libraries'
# patch...
# patch CMakeLists.txt ${CUOREPATCH}/sources/CMakeLists.txt.diff
#cd ../../source


#
# Retrieve and install esources
#
if [ $G4VERSION -ge 10 ]; then
    echo "Skipping esources patch for this version of Geant4.$G4VERSION"
else
    echo "Retrieving esources from svn repository..."
    exit -1
    cd ${G4PATH}/source
    #echo "svn co svn+ssh://$CUOREUSER@$CUORESVN/esources"
    svn co svn+ssh://$CUOREUSER@$CUORESVN/esources
    cp $CUOREPATCH/sources.cmake $G4PATH/source/esources/	
    # patch to update esources to include more files in newer geant4
fi



#vi CMakeLists.txt 
#add line: add_subdirectory(esources)
# patch CMakeLists.txt $PATCH_DIR/G4.$G4VERSION.$G4RELEASE-CMakeLists-sources-patch.diff
#cd ../../
#mv geant4.9.6.b01  $G4SRC/v9.6
#cd $G4SRC/v$VERSION
#mkdir geant4.$VERSION.$RELEASE-build
#cd geant4.9.6.b01-build

if [ -d $G4BUILD ]; then
	echo "Removed old build directory."
	rm -rf $G4BUILD
fi
mkdir -vp $G4BUILD 
mkdir -vp $G4INSTALL
cd $G4BUILD
cmake -DGEANT4_INSTALL_DATA=ON -DCMAKE_INSTALL_PREFIX=$G4INSTALL $G4PATH
#(package expat-devel is needed. If not present and you get an error like "Could NOT find EXPAT (missing: EXPAT_LIBRARY EXPAT_INCLUDE_DIR)" make yum list expat-devel and install the outgoing package)
cd $G4BUILD

make -j"$G4NPROC"
#exit 0
make install






#
#
#  SCRIPT EXIT
#
#
exit 0










#
# FUKUI e VRML
#
#
# moved to patch
#cd $G4SRC/v9.6/geant4.9.6.${G4RELEASE}/source/visualization
export GEANT4_USE_NETWORKVRML=1
export GEANT4_USE_NETWORKDAWN=1
# append to CMakeLists.txt:
# now in patch file
# list(APPEND Geant4_DEFINITIONS -DG4VIS_USE_VRML -DG4VIS_USE_VRMLFILE -DG4VIS_USE_DAWN -DG4VIS_USE_DAWNFILE)
# add_definitions(${Geant4_DEFINITIONS})

#cd VRML
# comment out in sources.cmake:
#if(GEANT4_USE_NETWORKVRML)
#Endif()
# (these 2 lines are enough!!)
# append to CMakeLists.txt:
#list(APPEND Geant4_DEFINITIONS -DG4VIS_USE_VRML -DG4VIS_USE_VRMLFILE -DG4VIS_USE_DAWN -DG4VIS_USE_DAWNFILE)
#add_definitions(${Geant4_DEFINITIONS})

#cd ../FukuiRenderer
# comment out in sources.cmake:
#if(GEANT4_USE_NETWORKDAWN)
# (these 2 lines are enough!!)
# append to CMakeLists.txt:
#list(APPEND Geant4_DEFINITIONS -DG4VIS_USE_VRML -DG4VIS_USE_VRMLFILE -DG4VIS_USE_DAWN -DG4VIS_USE_DAWNFILE)
#add_definitions(${Geant4_DEFINITIONS})

#cd  $G4BUILD
#make

## For private installations:
#yum install gtkglext
#yum install gtkglext-devel
#cd $G4SRC/applications
## download view3dscene from http://castle-engine.sourceforge.net/view3dscene.php
#cp $G4SRC/applications/view3dscene/view3dscene $CUORE_BASE_PATH/bin


#
# Applications
#
# For FARM installations:
#cd $CUOREBASE
#cp geant4/v9.5/G4.9.5exports bin/setup-g4-9.5.sh
#cd $G4SRC/applications/APPLICATION_NAME 
#cd $G4SRC/applications/APPLICATION_NAME 
# edit GNUmakefile and comment out "EXTRALIBS += -l$(G4PL)"
#mkdir build
#cd build
#cmake -DGeant4_DIR=$G4BUILD .. 
#cmake -DCMAKE_BUILD_TYPE:STRING=Debug .. # if debug info is needed
#make
#install exec $G4BIN/
#cd ..
#rm -rf build

# For private installations:
#vrml visualization
#if $OS="fedora"
#	yum install gtkglext
#	yum install gtkglext-devel
#else if $OS="ubuntu"
#	sudo apt-get install gtkglext1-dev
#fi


#cd $G4SRC/applications
# download view3dscene from http://castle-engine.sourceforge.net/view3dscene.php
#cp $G4SRC/applications/view3dscene/view3dscene $CUORE_BASE_DIR/bin


#
# CUORE geometry codes
#
mkdir -vp ${CUORESRC}
mkdir -vp ${CUOREBIN}
for PROGRAM in qshields g2tas #mcuoricino
do
	cd ${CUORESRC}
	svn co ${CUORESVN}/$PROGRAM
done

for PROGRAM in qshields #mcuoricino
do
	mkdir -vp ${CUOREBUILD}/$PROGRAM-build
	cd ${CUOREBUILD}/$PROGRAM-build
	cmake -DGeant4_DIR=${G4BUILD} ${CUORESRC}/$PROGRAM
	make -j"$(nproc)"
	install $PROGRAM ${CUOREBIN}/
	rm -rf ${CUOREBUILD}/${PROGRAM}-build
done


exit 0


#
# g2tas
#
# edit Makefile.gcc and change BINDIR to $CUOREBIN
# check the definition of the dirs containing lib files and eventually change them
cd ${CUOREBUILD}/g2tas
make -f  Makefile.gcc
#make -f  Makefile.gcc install
