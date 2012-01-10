# Copyright (C) 2009,2010 The ESPResSo project
# Copyright (C) 2009,2010 Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany

# This file is part of ESPResSo.

# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

cmake_minimum_required(VERSION 2.6)

project(ESPResSo)

# Cmake modules/macros are in a subdirectory to keep this file cleaner
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE)


########################################################################
# User input options                                                   #
########################################################################
option(WITH_MPI    "Build a parallel (message-passing) version of ESPResSo" OFF)
option(WITH_TK     "Build with tk support" OFF)
set(MYCONFIG "myconfig.h" CACHE STRING "default name of the local config file")

enable_language(C)

find_package(FFTW3)
if (FFTW3_FOUND)
    include_directories(${FFTW3_INCLUDE_DIR})
    list(APPEND LIBRARIES ${FFTW3_LIBRARIES})
    set(FFTW 3)
else(FFTW3_FOUND)
    message(STATUS "FFTW3 not found, diabling it.")
endif(FFTW3_FOUND)

find_package(TCL REQUIRED)
include_directories(${TCL_INCLUDE_PATH})
list(APPEND LIBRARIES ${TCL_LIBRARY})

########################################################################
# Fix stupid flags on MSVC
########################################################################
IF(CMAKE_GENERATOR MATCHES "Visual Studio")
    STRING(REPLACE /MD /MT CMAKE_C_FLAGS_RELEASE ${CMAKE_C_FLAGS_RELEASE})
    SET(CMAKE_C_FLAGS_RELEASE ${CMAKE_C_FLAGS_RELEASE} CACHE STRING "" FORCE)
    STRING(REPLACE /MD /MT CMAKE_C_FLAGS_DEBUG ${CMAKE_C_FLAGS_DEBUG})
    SET(CMAKE_C_FLAGS_DEBUG ${CMAKE_C_FLAGS_DEBUG} CACHE STRING "" FORCE)   
ENDIF(CMAKE_GENERATOR MATCHES "Visual Studio")

########################################################################
#Process MPI settings
########################################################################
if(WITH_MPI)
  find_package(MPI)
  if(MPI_FOUND)
    list(APPEND CMAKE_C_FLAGS ${MPI_COMPILE_FLAGS})
    list(APPEND CMAKE_EXE_LINKER_FLAGS ${MPI_LINK_FLAGS})
    include_directories(${MPI_INCLUDE_PATH})
    list(APPEND LIBRARIES ${MPI_LIBRARIES})
    set(MPI "external")
  else(MPI_FOUND)
    message(FATAL_ERROR "MPI support requested, but no MPI compiler found.")
  endif(MPI_FOUND)
else(WITH_MPI)
  include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/mpifake)  
  set(MPI "fake")
  set(EXTRA_SOURCE mpifake/mpi.h mpifake/mpi.c)
endif(WITH_MPI)

include(EspTestInline)
esp_test_inline(MDINLINE)

include(EspGenerateVersionInfo)
esp_generate_version_info(PACKAGE_VERSION LAST_CHANGE)

########################################################################
# Basic system tests (standard libraries, headers, functions, types)   #
########################################################################
message(STATUS "Using local config file: ${MYCONFIG}")
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/cmake)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/myconfig-default.h.cmakein ${CMAKE_CURRENT_SOURCE_DIR}/cmake/${MYCONFIG})

include(CheckIncludeFiles)
check_include_files("arpa/inet.h;ctype.h;errno.h;fcntl.h;float.h;inttypes.h;limits.h;math.h;netdb.h;netinet/in.h;signal.h;stddef.h;stdio.h;stdlib.h;string.h;strings.h;sys/file.h;sys/resource.h;sys/socket.h;sys/time.h;sys/types.h;time.h;unistd.h" HAVE_HEADERS)
if (NOT HAVE_HEADERS)
  message(FATAL_ERROR "There is a least one header missing.")
endif(NOT HAVE_HEADERS)

include(CheckFunctionExists)
check_function_exists("atexit;floor;getcwd;gethostbyname;memset;pow;rint;select;socket;sqrt"       HAVE_FUNCS)

include(CheckLibraryExists)
check_library_exists(m sqrt "" HAVE_LIBM)

include(CheckTypeSize)
check_type_size(size_t size_t)
if(NOT size_t)
    set(size_t int)
else(NOT size_t)
    set(size_t 0)
endif(NOT size_t)

include(EspSignalType)
ESP_TEST_RETSIGTYPE(RETSIGTYPE)

# Math and thread libraries must often come after all others when linking...
if(HAVE_LIBM)
    list(APPEND	LIBRARIES m)
endif(HAVE_LIBM)

########################################################################
# Specify install locations and which subdirectories to process        #
########################################################################
set(BIN_INSTALL_DIR  ${CMAKE_INSTALL_PREFIX}/bin)

#add_subdirectory(doc)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/acconfig.h.cmakein ${CMAKE_CURRENT_BINARY_DIR}/acconfig.h)
include_directories(${CMAKE_CURRENT_BINARY_DIR})

set(ESPRESSO_BIN "${CMAKE_CURRENT_BINARY_DIR}/Espresso_bin")
set(ESPRESSO_SCRIPTS "${CMAKE_CURRENT_SOURCE_DIR}/scripts")
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/Espresso.cmakein ${CMAKE_CURRENT_BINARY_DIR}/Espresso)

add_executable(Espresso_bin
  main.c config.c config.h initialize.c initialize.h global.c global.h communication.c communication.h binary_file.c binary_file.h 
  interaction_data.c interaction_data.h verlet.c verlet.h grid.c grid.h integrate.c integrate.h cells.c cells.h ghosts.c ghosts.h
  forces.c forces.h rotation.c rotation.h debug.c debug.h particle_data.c particle_data.h thermostat.c thermostat.h dpd.c dpd.h
  statistics.c statistics.h statistics_chain.c statistics_chain.h energy.c energy.h pressure.c pressure.h vmdsock.c vmdsock.h
  imd.c imd.h iccp3m.c iccp3m.h p3m.c p3m.h magnetic_non_p3m_methods.c magnetic_non_p3m_methods.h ewald.c ewald.h fft.c fft.h
  random.c random.h blockfile.c blockfile.h blockfile_tcl.c blockfile_tcl.h polymer.c polymer.h specfunc.c specfunc.h tuning.c tuning.h
  uwerr.c uwerr.h parser.c parser.h domain_decomposition.c domain_decomposition.h nsquare.c nsquare.h layered.c layered.h mmm-common.c 
  mmm-common.h mmm1d.c mmm1d.h mmm2d.c mmm2d.h modes.c modes.h topology.c topology.h nemd.c nemd.h statistics_cluster.c statistics_cluster.h
  elc.c elc.h mdlc_correction.c  mdlc_correction.h statistics_molecule.c statistics_molecule.h errorhandling.c errorhandling.h constraint.c
  constraint.h maggs.c maggs.h mol_cut.h rattle.c rattle.h molforces.c molforces.h virtual_sites.c virtual_sites.h metadynamics.c metadynamics.h
  lb.c lb.h lb-d3q18.h lb-d3q19.h bin.c bin.h lattice.c lattice.h halo.c halo.h statistics_fluid.c statistics_fluid.h lb-boundaries.c 
  lb-boundaries.h utils.c utils.h angle.h pwdist.h angledist.h endangledist.h buckingham.h comfixed.h comforce.h debye_hueckel.h reaction_field.h
  dihedral.h fene.h gb.h harmonic.h imd.h ljcos2.h ljcos.h lj.h ljgen.c ljgen.h steppot.h bmhtf-nacl.h morse.h polynom.h soft_sphere.h subt_lj.h 
  tab.h overlap.h ljangle.h adresso.c adresso.h tunable_slip.h ${MYCONFIG} ${EXTRA_SOURCE})

target_link_libraries(Espresso_bin ${LIBRARIES})
