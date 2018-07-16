#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(LXe)

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
#
include(${Geant4_USE_FILE})

find_package(cadmesh)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(LXe LXe.cc ${sources} ${headers})
include_directories(${PROJECT_SOURCE_DIR}/include
                    ${PROJECT_SOURCE_DIR}/ds
                    ${Geant4_INCLUDE_DIR}
                    ${ROOT_INCLUDE_DIRS}
                    ${ROOTSYS}/include )
target_link_libraries(LXe ${Geant4_LIBRARIES})
target_link_libraries(LXe ${cadmesh_LIBRARIES})
target_link_libraries(LXe ${ROOT_LIBRARIES})
target_link_libraries(LXe ${PROJECT_SOURCE_DIR}/ds/libEvent.so)

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build LXe. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(LXe_SCRIPTS
    cerenkov.mac defaults.mac LXe.in LXe.out photon.mac reviewEvent.mac vis.mac wls.mac
  )

foreach(_script ${LXe_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS LXe DESTINATION bin)

