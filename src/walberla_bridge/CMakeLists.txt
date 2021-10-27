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
    generated_kernels/StreamSweep.cpp generated_kernels/CollideSweep.cpp
    generated_kernels/CollideSweepThermalized.cpp
    generated_kernels/InitialPDFsSetter.cpp generated_kernels/Dynamic_UBB.cpp)
if(WALBERLA_USE_AVX)
  list(APPEND WalberlaBridgeCodegen_SRC generated_kernels/CollideSweepAVX.cpp
       generated_kernels/CollideSweepThermalizedAVX.cpp)
endif()
set(WalberlaBridge_SRC ${WalberlaBridgeCodegen_SRC} lb_walberla_init.cpp)
add_library(WalberlaBridge SHARED ${WalberlaBridge_SRC})
target_link_libraries(WalberlaBridge PUBLIC MPI::MPI_CXX EspressoUtils
                      PRIVATE Boost::boost cxx_interface ${WALBERLA_LIBS})
target_include_directories(
  WalberlaBridge PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
  PRIVATE ${WALBERLA_INCLUDE_DIRS} ${walberla_BINARY_DIR}/src)

install(TARGETS WalberlaBridge LIBRARY DESTINATION ${PYTHON_INSTDIR}/espressomd)

list(JOIN WalberlaBridgeCodegen_SRC "\",\n    \""
     WalberlaBridgeCodegen_SRC_json)
set(WalberlaBridgeCodegenJson
    "{\n\
  \"EXPECTED_FILES\": [\n\
    \"${WalberlaBridgeCodegen_SRC_json}\"\n\
  ],\n\
  \"CMAKE_VARS\" : {\n\
    \"WALBERLA_OPTIMIZE_FOR_LOCALHOST\": \"${WALBERLA_OPTIMIZE_FOR_LOCALHOST}\",\n\
    \"WALBERLA_DOUBLE_ACCURACY\": \"${WALBERLA_DOUBLE_ACCURACY}\",\n\
    \"CODEGEN_CFG\": \"${codegenCfg}\",\n\
    \"WALBERLA_BUILD_WITH_MPI\": \"${WALBERLA_BUILD_WITH_MPI}\",\n\
    \"WALBERLA_BUILD_WITH_CUDA\": \"${WALBERLA_BUILD_WITH_CUDA}\",\n\
    \"WALBERLA_BUILD_WITH_OPENMP\": \"${WALBERLA_BUILD_WITH_OPENMP}\"\n\
  }\n\
}")
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/codegenCfg.json"
     ${WalberlaBridgeCodegenJson})

if(WITH_TESTS)
  add_subdirectory(tests)
endif(WITH_TESTS)
