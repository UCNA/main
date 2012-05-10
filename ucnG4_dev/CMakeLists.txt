#----------------------------------------------------------------------------
# Usage:
# cmake -DGeant4_DIR=~/geant4.9.5/geant4.9.5-install/lib/Geant4-9.5.0/ ../ucnG4_dev/
# make -j1 VERBOSE=1

#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(ucnG4_prod)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc  ${PROJECT_SOURCE_DIR}/../Physics/*.cc)
LIST(APPEND sources ${PROJECT_SOURCE_DIR}/../IOUtils/strutils.cc)
LIST(APPEND sources ${PROJECT_SOURCE_DIR}/../IOUtils/SMExcept.cc)
LIST(APPEND sources ${PROJECT_SOURCE_DIR}/../IOUtils/QFile.cc)
LIST(APPEND sources ${PROJECT_SOURCE_DIR}/../IOUtils/PathUtils.cc)
LIST(APPEND sources ${PROJECT_SOURCE_DIR}/../BaseTypes/Enums.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh ${PROJECT_SOURCE_DIR}/../Physics/*.hh)
LIST(APPEND sources ${PROJECT_SOURCE_DIR}/../IOUtils/strutils.hh)
LIST(APPEND sources ${PROJECT_SOURCE_DIR}/../IOUtils/SMExcept.hh)
LIST(APPEND sources ${PROJECT_SOURCE_DIR}/../IOUtils/QFile.hh)
LIST(APPEND sources ${PROJECT_SOURCE_DIR}/../IOUtils/PathUtils.hh)
LIST(APPEND headers ${PROJECT_SOURCE_DIR}/../BaseTypes/Enums.hh)

#-------------------------------
# ROOT dependencies
#-------------------------------

######## locate ROOT install, headers, includes
include(FindROOT.cmake)
include_directories(${ROOT_INCLUDE_DIR})
link_directories(${ROOT_LIBRARY_DIR})

######## set up output paths for libraries
set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/lib")
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")

######## set everything needed for the root dictonary and create the dictionary
set(ucnG4_prod_LINKDEF ${PROJECT_SOURCE_DIR}/ucnG4_prod_LinkDef.h)
set(ucnG4_prod_DICTIONARY ${CMAKE_CURRENT_BINARY_DIR}/ucnG4_prodDict.cxx) 
MESSAGE("Generating dictionary ${ucnG4_prod_DICTIONARY}")
set(DICTFILES "${PROJECT_SOURCE_DIR}/include/bmMCEvent.hh;${PROJECT_SOURCE_DIR}/include/bmTrackInfo.hh;${PROJECT_SOURCE_DIR}/include/bmPrimaryInfo.hh")
ROOT_GENERATE_DICTIONARY("${DICTFILES}" "${ucnG4_prod_LINKDEF}" "${ucnG4_prod_DICTIONARY}" "${ROOT_INCLUDE_DIR};${PROJECT_SOURCE_DIR}/include/;${PROJECT_SOURCE_DIR}/../Physics/")
LIST(APPEND sources ${ucnG4_prod_DICTIONARY})
MESSAGE(sources = ${sources})

############### build the library #####################
#add_library(ucnG4_prod_rootlibs SHARED ${sources})
add_library(ucnG4_prod_rootlibs SHARED ${DICTFILES} ${ucnG4_prod_DICTIONARY})
target_link_libraries(ucnG4_prod_rootlibs ${ROOT_LIBRARIES})
SET(ucnG4_prod_MAJOR_VERSION 0)
SET(ucnG4_prod_MINOR_VERSION 0)
SET(ucnG4_prod_PATCH_VERSION 0)
SET(ucnG4_prod_VERSION "${ucnG4_prod_MAJOR_VERSION}.${ucnG4_prod_MINOR_VERSION}.${ucnG4_prod_PATCH_VERSION}")
SET(ucnG4_prod_LIBRARY_PROPERTIES ${ucnG4_prod_LIBRARY_PROPERTIES}
    VERSION "${ucnG4_prod_VERSION}"
    SOVERSION "${ucnG4_prod_MAJOR_VERSION}"
    SUFFIX ".so"
)
set_target_properties(ucnG4_prod_rootlibs PROPERTIES ${ucnG4_prod_LIBRARY_PROPERTIES})
############### install the library ###################
install(TARGETS ucnG4_prod_rootlibs DESTINATION ${CMAKE_BINARY_DIR}/lib)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()


#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/../Physics/ ${PROJECT_SOURCE_DIR}/../BaseTypes/ ${PROJECT_SOURCE_DIR}/../IOUtils/)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(ucnG4_prod ucnG4_prod.cc ${sources} ${headers})
target_link_libraries(ucnG4_prod ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})
add_executable(UCNA_MC_Analyzer UCNA_MC_Analyzer.cc ${sources} ${headers})
target_link_libraries(UCNA_MC_Analyzer ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})
add_executable(SiDet_Analyzer SiDet_Analyzer.cc ${sources} ${headers})
target_link_libraries(SiDet_Analyzer ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS ucnG4_prod UCNA_MC_Analyzer SiDet_Analyzer DESTINATION bin)

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
set(CMAKE_BUILD_TYPE release)
