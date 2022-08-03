#
# Copyright (C) 2020-2022 The ESPResSo project
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
if(WITH_CUDA AND WALBERLA_BUILD_WITH_CUDA)
  list(
    APPEND
    WalberlaBridgeCodegenCuda_SRC
    generated_kernels/CudaStreamSweepSinglePrecision.cu
    generated_kernels/CudaStreamSweepDoublePrecision.cu
    generated_kernels/CudaInitialPDFsSetterSinglePrecision.cu
    generated_kernels/CudaInitialPDFsSetterDoublePrecision.cu
    generated_kernels/CudaDynamic_UBB_single_precision.cu
    generated_kernels/CudaDynamic_UBB_double_precision.cu
    generated_kernels/CudaUpdateVelocityFromPDFSweepDoublePrecision.cu
    generated_kernels/CudaUpdateVelocityFromPDFSweepSinglePrecision.cu)
endif()

foreach(precision double_precision single_precision)
  list(
    APPEND
    WalberlaBridgeCodegen_SRC
    electrokinetics/generated_kernels/DiffusiveFluxKernel_${precision}.cpp
    electrokinetics/generated_kernels/DiffusiveFluxKernelWithElectrostatic_${precision}.cpp
    electrokinetics/generated_kernels/ContinuityKernel_${precision}.cpp
    electrokinetics/generated_kernels/AdvectiveFluxKernel_${precision}.cpp
    electrokinetics/generated_kernels/FrictionCouplingKernel_${precision}.cpp
    electrokinetics/generated_kernels/FixedFlux_${precision}.cpp
    electrokinetics/generated_kernels/Dirichlet_${precision}.cpp)
endforeach()

# EK reactions
foreach(n_species RANGE 1 5)
  foreach(precision double_precision single_precision)
    list(
      APPEND
      WalberlaBridgeCodegen_SRC
      electrokinetics/reactions/generated_kernels/ReactionKernelBulk_${n_species}_${precision}.cpp
      electrokinetics/reactions/generated_kernels/ReactionKernelIndexed_${n_species}_${precision}.cpp
    )
  endforeach()
endforeach()

add_library(Espresso_walberla_cpp_flags INTERFACE)
add_library(Espresso::walberla::cpp_flags ALIAS Espresso_walberla_cpp_flags)
if(WALBERLA_USE_AVX)
  target_link_libraries(Espresso_walberla_cpp_flags
                        INTERFACE Espresso::avx_flags)
endif(WALBERLA_USE_AVX)

add_library(
  Espresso_walberla SHARED
  lb_walberla_init.cpp electrokinetics/ek_walberla_init.cpp LatticeModel.cpp
  LatticeWalberla.cpp ${WalberlaBridgeCodegen_SRC}
  electrokinetics/reactions/EKReactionImplBulk.cpp
  electrokinetics/reactions/EKReactionIndexed.cpp)
add_library(Espresso::walberla ALIAS Espresso_walberla)

if(WITH_CUDA AND WALBERLA_BUILD_WITH_CUDA)
  add_gpu_library(Espresso_walberla_cuda SHARED
                  ${WalberlaBridgeCodegenCuda_SRC})
  add_library(Espresso::walberla_cuda ALIAS Espresso_walberla_cuda)
  target_link_libraries(Espresso_walberla PUBLIC Espresso::walberla_cuda)
endif()

set_source_files_properties(
  electrokinetics/generated_kernels/Dirichlet_single_precision.cpp
  electrokinetics/ek_walberla_init.cpp
  electrokinetics/reactions/EKReactionIndexed.cpp
  PROPERTIES COMPILE_FLAGS -Wno-unused-variable)

target_link_libraries(
  Espresso_walberla PUBLIC MPI::MPI_CXX Espresso::utils
  PRIVATE Boost::boost Espresso::cpp_flags Espresso::walberla::cpp_flags
          ${WALBERLA_LIBS})
target_include_directories(
  Espresso_walberla PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${WALBERLA_INCLUDE_DIRS}
  PRIVATE ${walberla_BINARY_DIR}/src)

install(TARGETS Espresso_walberla
        LIBRARY DESTINATION ${ESPRESSO_INSTALL_PYTHON}/espressomd)

if(WITH_TESTS)
  add_subdirectory(tests)
endif(WITH_TESTS)
