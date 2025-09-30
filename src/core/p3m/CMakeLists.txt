#
# Copyright (C) 2018-2025 The ESPResSo project
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

if(ESPRESSO_BUILD_WITH_FFTW)
  add_library(espresso_p3m OBJECT)
  add_library(espresso::p3m ALIAS espresso_p3m)
  target_include_directories(espresso_p3m
                             PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/..)
  espresso_set_common_target_properties(espresso_p3m)
  set_target_properties(
    espresso_p3m PROPERTIES CXX_CLANG_TIDY "${ESPRESSO_P3M_CXX_CLANG_TIDY}")
  add_library(espresso_p3m_compiler_flags INTERFACE)
  add_library(espresso::p3m::compiler_flags ALIAS espresso_p3m_compiler_flags)
  target_compile_options(
    espresso_p3m_compiler_flags
    INTERFACE
      $<$<COMPILE_LANGUAGE:CXX>:-Wno-float-conversion>
      $<$<COMPILE_LANG_AND_ID:CXX,Clang,AppleClang,IntelLLVM>:-Wno-implicit-float-conversion>
  )
  target_link_libraries(
    espresso_p3m
    PRIVATE
      espresso::compiler_flags
      espresso::p3m::compiler_flags
      espresso::instrumentation
      espresso::utils
      espresso::config
      Random123
      Boost::mpi
      MPI::MPI_CXX
      Heffte::Heffte
      $<$<BOOL:${ESPRESSO_BUILD_WITH_SHARED_MEMORY_PARALLELISM}>:OpenMP::OpenMP_CXX>
      $<$<BOOL:${ESPRESSO_BUILD_WITH_SHARED_MEMORY_PARALLELISM}>:Kokkos::kokkos>
      $<$<BOOL:${ESPRESSO_BUILD_WITH_SHARED_MEMORY_PARALLELISM}>:Cabana::Core>)
  set_target_properties(espresso_p3m PROPERTIES POSITION_INDEPENDENT_CODE ON)
  target_link_libraries(espresso_core PRIVATE espresso_p3m)
  target_sources(
    espresso_p3m PRIVATE common.cpp send_mesh.cpp TuningAlgorithm.cpp
                         FFTBackendLegacy.cpp FFTBuffersLegacy.cpp)
endif()
