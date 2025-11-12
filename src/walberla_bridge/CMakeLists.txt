#
# Copyright (C) 2020-2024 The ESPResSo project
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

set(SKIP_CLANG_TIDY_CHECKS "")
set(SKIP_CLANG_TIDY_CHECKS_CXX "")
set(SKIP_CLANG_TIDY_CHECKS_CUDA "")
# silence waLBerla diagnostics
list(APPEND SKIP_CLANG_TIDY_CHECKS "-readability-avoid-const-params-in-decls")
list(APPEND SKIP_CLANG_TIDY_CHECKS "-readability-else-after-return")
list(APPEND SKIP_CLANG_TIDY_CHECKS "-readability-simplify-boolean-expr")
list(APPEND SKIP_CLANG_TIDY_CHECKS "-modernize-pass-by-value")
list(APPEND SKIP_CLANG_TIDY_CHECKS "-bugprone-crtp-constructor-accessibility")
list(APPEND SKIP_CLANG_TIDY_CHECKS "-bugprone-narrowing-conversions")
list(APPEND SKIP_CLANG_TIDY_CHECKS "-bugprone-exception-escape")
list(APPEND SKIP_CLANG_TIDY_CHECKS "-bugprone-branch-clone")
if(WALBERLA_BUILD_WITH_CUDA)
  # silence diagnostics from cuda header files
  list(APPEND SKIP_CLANG_TIDY_CHECKS "-bugprone-casting-through-void")
  list(APPEND SKIP_CLANG_TIDY_CHECKS "-modernize-redundant-void-arg")
  list(APPEND SKIP_CLANG_TIDY_CHECKS "-modernize-use-nullptr")
  list(APPEND SKIP_CLANG_TIDY_CHECKS "-readability-non-const-parameter")
  # silence nullptr dereference in cuda::thrust
  list(APPEND SKIP_CLANG_TIDY_CHECKS_CUDA
       "-clang-analyzer-core.NonNullParamChecker")
endif()
# codegen-specific Clang-Tidy overrides
list(APPEND SKIP_CLANG_TIDY_CHECKS "-modernize-use-auto")
list(APPEND SKIP_CLANG_TIDY_CHECKS "-readability-non-const-parameter")
list(APPEND SKIP_CLANG_TIDY_CHECKS_CXX "-clang-analyzer-deadcode.DeadStores")
list(APPEND SKIP_CLANG_TIDY_CHECKS_CUDA
     "-bugprone-multi-level-implicit-pointer-conversion")

foreach(LANG CXX CUDA)
  set(ESPRESSO_WALBERLA_${LANG}_CLANG_TIDY_CODEGEN
      "${ESPRESSO_${LANG}_CLANG_TIDY}")
  espresso_override_clang_tidy_checks(
    ESPRESSO_WALBERLA_${LANG}_CLANG_TIDY_CODEGEN "${SKIP_CLANG_TIDY_CHECKS}"
    "${SKIP_CLANG_TIDY_CHECKS_${LANG}}")
  unset(SKIP_CLANG_TIDY_CHECKS_${LANG})
endforeach()
unset(SKIP_CLANG_TIDY_CHECKS)

function(espresso_configure_walberla_target)
  set(TARGET_NAME ${ARGV0})
  espresso_set_common_target_properties(${TARGET_NAME})
  target_link_libraries(${TARGET_NAME} PRIVATE ${WALBERLA_LIBS})
  target_include_directories(
    ${TARGET_NAME} PUBLIC include PRIVATE ${WALBERLA_INCLUDE_DIRS}
                                          ${walberla_BINARY_DIR}/src)
  if(ESPRESSO_BUILD_WITH_CLANG_TIDY)
    foreach(LANG CXX CUDA)
      set(TARGET_CLANG_TIDY "${ESPRESSO_${LANG}_CLANG_TIDY}")
      if(TARGET_NAME MATCHES "_codegen")
        set(TARGET_CLANG_TIDY "${ESPRESSO_WALBERLA_${LANG}_CLANG_TIDY_CODEGEN}")
      endif()
      set_target_properties(${TARGET_NAME} PROPERTIES ${LANG}_CLANG_TIDY
                                                      "${TARGET_CLANG_TIDY}")
    endforeach()
  endif()
  get_target_property(TARGET_TYPE ${TARGET_NAME} TYPE)
  if(TARGET_TYPE STREQUAL "OBJECT_LIBRARY")
    set_target_properties(${TARGET_NAME} PROPERTIES POSITION_INDEPENDENT_CODE
                                                    ON)
  else()
    install(TARGETS ${TARGET_NAME}
            LIBRARY DESTINATION ${ESPRESSO_INSTALL_PYTHON}/espressomd)
  endif()
endfunction()

add_library(espresso_walberla_compiler_flags INTERFACE)
add_library(espresso::walberla::compiler_flags ALIAS
            espresso_walberla_compiler_flags)
target_link_libraries(
  espresso_walberla_compiler_flags
  INTERFACE
    espresso::compiler_flags
    $<$<BOOL:${ESPRESSO_BUILD_WITH_WALBERLA_AVX}>:espresso::avx_flags>
    $<$<BOOL:${ESPRESSO_BUILD_WITH_SHARED_MEMORY_PARALLELISM}>:OpenMP::OpenMP_CXX>
)
add_library(espresso_walberla_codegen_compiler_flags INTERFACE)
add_library(espresso::walberla_codegen::compiler_flags ALIAS
            espresso_walberla_codegen_compiler_flags)
# cmake-format: off
target_compile_options(
  espresso_walberla_codegen_compiler_flags
  INTERFACE
  $<$<COMPILE_LANGUAGE:CXX>:-Wno-cast-align>
  $<$<COMPILE_LANG_AND_ID:CXX,GNU>:-Wno-shadow=compatible-local>
  # silence diagnostics from MPI vendors
  $<$<AND:$<STREQUAL:${ESPRESSO_MPIEXEC_VENDOR},OpenMPI>,$<COMPILE_LANG_AND_ID:CXX,Clang,IntelLLVM,AppleClang>>:-Wno-cast-function-type-mismatch>
  $<$<AND:$<STREQUAL:${ESPRESSO_MPIEXEC_VENDOR},OpenMPI>,$<COMPILE_LANG_AND_ID:CUDA,Clang,IntelLLVM>>:-Wno-cast-function-type-mismatch>
  $<$<AND:$<STREQUAL:${ESPRESSO_MPIEXEC_VENDOR},OpenMPI>,$<COMPILE_LANG_AND_ID:CXX,GNU>>:-Wno-cast-function-type>
  $<$<AND:$<STREQUAL:${ESPRESSO_MPIEXEC_VENDOR},OpenMPI>,$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>>:-Xcompiler=-Wno-cast-function-type>
  $<$<AND:$<STREQUAL:${ESPRESSO_MPIEXEC_VENDOR},MPICH>,$<COMPILE_LANGUAGE:CXX>>:-Wno-cast-qual>
)
# cmake-format: on

add_library(espresso_walberla SHARED)
add_library(espresso_walberla_codegen OBJECT)
add_library(espresso::walberla ALIAS espresso_walberla)
add_library(espresso::walberla_codegen ALIAS espresso_walberla_codegen)

espresso_configure_walberla_target(espresso_walberla)
espresso_configure_walberla_target(espresso_walberla_codegen)

target_link_libraries(
  espresso_walberla PUBLIC MPI::MPI_CXX espresso::utils
  PRIVATE espresso::walberla::compiler_flags espresso::walberla_codegen
          $<$<BOOL:${ESPRESSO_BUILD_WITH_FFTW}>:Heffte::Heffte>)
target_link_libraries(
  espresso_walberla_codegen PRIVATE espresso::walberla::compiler_flags
                                    espresso::walberla_codegen::compiler_flags)

function(espresso_set_walberla_codegen_kernel_source_properties source_file)
  # disable sanitizers on kernels due to excessive runtime and VRAM usage
  set_source_files_properties(
    ${source_file} TARGET_DIRECTORY espresso_walberla_codegen
    PROPERTIES
      COMPILE_FLAGS
      "-fno-sanitize=undefined $<$<CONFIG:Debug>:$<IF:$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>,-Xptxas=-O0 -Xcompiler=-Og,-Og>>"
  )
endfunction()

if(WALBERLA_BUILD_WITH_CUDA)
  target_link_libraries(espresso_walberla PRIVATE CUDA::cuda_driver
                                                  CUDA::cudart)
  target_link_libraries(espresso_walberla_codegen PRIVATE CUDA::cudart)
endif()

add_subdirectory(src)
if(ESPRESSO_BUILD_TESTS)
  add_subdirectory(tests)
endif()
