#----------------------------------------------------------------------------
# Usage:
# cmake -DGeant4_DIR=~/geant4.9.5/geant4.9.5-install/lib/Geant4-9.5.0/ ../ucnG4_dev/
# cmake -DGeant4_DIR=/Users/michael/Applications/geant4/lib/Geant4-9.5.1/ ../UCNA_git/ucnG4_dev/
# make

#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(ucnG4_prod)

# build libUCNA.a dependency
add_custom_target(
	extern_lib
	COMMAND make
	WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/..
)


#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)
include_directories(${PROJECT_SOURCE_DIR}/include
					${PROJECT_SOURCE_DIR}/../Physics
					${PROJECT_SOURCE_DIR}/../IOUtils
					${PROJECT_SOURCE_DIR}/../BaseTypes)
link_directories(${PROJECT_SOURCE_DIR}/..)

#-------------------------------
# ROOT dependencies
#-------------------------------

######## locate ROOT install, headers, includes
include(FindROOT.cmake)
include_directories(${ROOT_INCLUDE_DIR})
link_directories(${ROOT_LIBRARY_DIR})

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
include(${Geant4_USE_FILE})


##############
# Output paths
##############

set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/lib")
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")


###############
# TObject subclasses dictionary & ucnG4 shared library
###############

set(ucnG4_LINKDEF ${PROJECT_SOURCE_DIR}/ucnG4_LinkDef.h)
set(ucnG4_DICTIONARY ${CMAKE_CURRENT_BINARY_DIR}/ucnG4_Dict.cxx) 
MESSAGE("Generating dictionary ${ucnG4_DICTIONARY}")
set(DICTFILES "${PROJECT_SOURCE_DIR}/include/bmMCEvent.hh;${PROJECT_SOURCE_DIR}/include/bmTrackInfo.hh;${PROJECT_SOURCE_DIR}/include/bmPrimaryInfo.hh")
ROOT_GENERATE_DICTIONARY("${DICTFILES}" "${ucnG4_LINKDEF}" "${ucnG4_DICTIONARY}" "")
LIST(APPEND sources ${ucnG4_DICTIONARY})

add_library(ucnG4 SHARED ${sources} ${headers})
add_dependencies(ucnG4 extern_lib)
target_link_libraries(ucnG4 ${ROOT_LIBRARIES} ${Geant4_LIBRARIES} UCNA)
SET(ucnG4_LIBRARY_PROPERTIES
		${ucnG4_LIBRARY_PROPERTIES}
		VERSION "0.0"
		SOVERSION "0"
		SUFFIX ".so"
)
set_target_properties(ucnG4 PROPERTIES ${ucnG4_LIBRARY_PROPERTIES})
install(TARGETS ucnG4 DESTINATION ${LIBRARY_OUTPUT_PATH})
link_directories(${LIBRARY_OUTPUT_PATH})


#################
# Executables
#################

add_executable(ucnG4_prod ucnG4_prod.cc)
target_link_libraries(ucnG4_prod ucnG4 ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} libUCNA.a)
add_dependencies(ucnG4_prod extern_lib ucnG4)

add_executable(UCNA_MC_Analyzer UCNA_MC_Analyzer.cc)
target_link_libraries(UCNA_MC_Analyzer ucnG4 ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} libUCNA.a)
add_dependencies(UCNA_MC_Analyzer extern_lib ucnG4)

add_executable(SiDet_Analyzer SiDet_Analyzer.cc)
target_link_libraries(SiDet_Analyzer ucnG4 ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} libUCNA.a)
add_dependencies(SiDet_Analyzer extern_lib ucnG4)

install(TARGETS ucnG4_prod UCNA_MC_Analyzer SiDet_Analyzer DESTINATION ${EXECUTABLE_OUTPUT_PATH})


#################
# Build settings
#################

#set(CMAKE_BUILD_TYPE debug)
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
set(CMAKE_BUILD_TYPE release)

