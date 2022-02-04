#
# Copyright (C) 2020 The ESPResSo project
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
set(WalberlaBridgeCodegen_SRC
    generated_kernels/StreamSweepSinglePrecision.cpp
    generated_kernels/StreamSweepDoublePrecision.cpp
    generated_kernels/InitialPDFsSetterSinglePrecision.cpp
    generated_kernels/InitialPDFsSetterDoublePrecision.cpp
    generated_kernels/LeesEdwards_Collision.cpp
    generated_kernels/Dynamic_UBB_single_precision.cpp
    generated_kernels/Dynamic_UBB_double_precision.cpp)
if(WALBERLA_USE_AVX)
  list(APPEND WalberlaBridgeCodegen_SRC
       generated_kernels/CollideSweepSinglePrecisionAVX.cpp
       generated_kernels/CollideSweepDoublePrecisionAVX.cpp
       generated_kernels/CollideSweepSinglePrecisionLeesEdwardsAVX.cpp
       generated_kernels/CollideSweepDoublePrecisionLeesEdwardsAVX.cpp
       generated_kernels/CollideSweepSinglePrecisionThermalizedAVX.cpp
       generated_kernels/CollideSweepDoublePrecisionThermalizedAVX.cpp)
else()
  list(APPEND WalberlaBridgeCodegen_SRC
       generated_kernels/CollideSweepSinglePrecision.cpp
       generated_kernels/CollideSweepDoublePrecision.cpp
       generated_kernels/CollideSweepSinglePrecisionLeesEdwards.cpp
       generated_kernels/CollideSweepDoublePrecisionLeesEdwards.cpp
       generated_kernels/CollideSweepSinglePrecisionThermalized.cpp
       generated_kernels/CollideSweepDoublePrecisionThermalized.cpp)
endif()
set(WalberlaBridge_SRC ${WalberlaBridgeCodegen_SRC} lb_walberla_init.cpp
                       LatticeWalberla.cpp)

add_library(cxx_walberla_interface INTERFACE)
if(WALBERLA_USE_AVX)
  target_link_libraries(cxx_walberla_interface INTERFACE cxx_avx2_interface)
endif(WALBERLA_USE_AVX)

add_library(WalberlaBridge SHARED ${WalberlaBridge_SRC})
target_link_libraries(WalberlaBridge PUBLIC MPI::MPI_CXX EspressoUtils
                      PRIVATE Boost::boost cxx_interface cxx_walberla_interface
                      ${WALBERLA_LIBS})
target_include_directories(
  WalberlaBridge PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${WALBERLA_INCLUDE_DIRS}
  PRIVATE ${walberla_BINARY_DIR}/src)

install(TARGETS WalberlaBridge LIBRARY DESTINATION ${PYTHON_INSTDIR}/espressomd)

if(WITH_TESTS)
  add_subdirectory(tests)
endif(WITH_TESTS)
