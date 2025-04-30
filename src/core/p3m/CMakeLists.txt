#
# Copyright (C) 2018-2024 The ESPResSo project
#
# This file is part of ESPResSo.
#
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

target_sources(
  espresso_core PRIVATE common.cpp send_mesh.cpp TuningAlgorithm.cpp
                        FFTBackendLegacy.cpp FFTBuffersLegacy.cpp)

if(ESPRESSO_BUILD_WITH_FFTW)
  add_library(espresso_p3m OBJECT)
  add_library(espresso::p3m ALIAS espresso_p3m)
  if(ESPRESSO_BUILD_WITH_CLANG_TIDY)
    set_target_properties(espresso_p3m PROPERTIES CXX_CLANG_TIDY
                                                  "${HEFFTE_CXX_CLANG_TIDY}")
  endif()
  add_library(espresso_p3m_cpp_flags INTERFACE)
  add_library(espresso::p3m::cpp_flags ALIAS espresso_p3m_cpp_flags)
  target_compile_options(
    espresso_p3m_cpp_flags
    INTERFACE
      -Wno-float-conversion
      $<$<CXX_COMPILER_ID:Clang,AppleClang,IntelLLVM>:-Wno-implicit-float-conversion>
  )
  target_link_libraries(
    espresso_p3m
    PRIVATE
      Heffte espresso::cpp_flags espresso::p3m::cpp_flags espresso::utils
      espresso::config espresso::instrumentation MPI::MPI_CXX Boost::mpi
      $<$<BOOL:${ESPRESSO_BUILD_WITH_SHARED_MEMORY_PARALLELISM}>:Kokkos::kokkos>
      $<$<BOOL:${ESPRESSO_BUILD_WITH_SHARED_MEMORY_PARALLELISM}>:Cabana::Core>)
  set_target_properties(espresso_p3m PROPERTIES POSITION_INDEPENDENT_CODE ON)
  target_link_libraries(espresso_core PUBLIC espresso_p3m)
  install(TARGETS espresso_p3m
          LIBRARY DESTINATION ${ESPRESSO_INSTALL_PYTHON}/espressomd)
endif()
