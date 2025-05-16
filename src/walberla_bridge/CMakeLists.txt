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

set(WALBERLA_CXX_CLANG_TIDY "${ESPRESSO_CXX_CLANG_TIDY}")
set(WALBERLA_CUDA_CLANG_TIDY "${ESPRESSO_CUDA_CLANG_TIDY}")
set(WALBERLA_CXX_CLANG_TIDY_CODEGEN "${ESPRESSO_CXX_CLANG_TIDY}")
set(WALBERLA_CUDA_CLANG_TIDY_CODEGEN "${ESPRESSO_CUDA_CLANG_TIDY}")
set(SKIP_CLANG_TIDY_CHECKS "")
set(SKIP_CLANG_TIDY_CHECKS_CXX "")
set(SKIP_CLANG_TIDY_CHECKS_CUDA "")
# silence waLBerla diagnostics
list(APPEND SKIP_CLANG_TIDY_CHECKS "-clang-analyzer-deadcode.DeadStores")
list(APPEND SKIP_CLANG_TIDY_CHECKS "-readability-non-const-parameter")
list(APPEND SKIP_CLANG_TIDY_CHECKS "-readability-avoid-const-params-in-decls")
list(APPEND SKIP_CLANG_TIDY_CHECKS "-readability-else-after-return")
list(APPEND SKIP_CLANG_TIDY_CHECKS "-readability-simplify-boolean-expr")
list(APPEND SKIP_CLANG_TIDY_CHECKS "-modernize-use-auto")
list(APPEND SKIP_CLANG_TIDY_CHECKS "-modernize-loop-convert")
list(APPEND SKIP_CLANG_TIDY_CHECKS "-modernize-pass-by-value")
list(APPEND SKIP_CLANG_TIDY_CHECKS "-modernize-use-equals-delete")
list(APPEND SKIP_CLANG_TIDY_CHECKS "-modernize-return-braced-init-list")
list(APPEND SKIP_CLANG_TIDY_CHECKS "-bugprone-crtp-constructor-accessibility")
list(APPEND SKIP_CLANG_TIDY_CHECKS "-bugprone-narrowing-conversions")
list(APPEND SKIP_CLANG_TIDY_CHECKS "-bugprone-exception-escape")
list(APPEND SKIP_CLANG_TIDY_CHECKS "-bugprone-branch-clone")
if(WALBERLA_BUILD_WITH_CUDA)
  # silence diagnostics from cuda header files
  list(APPEND SKIP_CLANG_TIDY_CHECKS "-bugprone-casting-through-void")
  list(APPEND SKIP_CLANG_TIDY_CHECKS "-modernize-redundant-void-arg")
  list(APPEND SKIP_CLANG_TIDY_CHECKS "-modernize-use-nullptr")
  # silence nullptr dereference in cuda::thrust
  list(APPEND SKIP_CLANG_TIDY_CHECKS_CUDA
       "-clang-analyzer-core.NonNullParamChecker")
endif()

if(ESPRESSO_BUILD_WITH_CLANG_TIDY)
  espresso_override_clang_tidy_checks(
    WALBERLA_CXX_CLANG_TIDY "${SKIP_CLANG_TIDY_CHECKS}"
    "${SKIP_CLANG_TIDY_CHECKS_CXX}")
  espresso_override_clang_tidy_checks(
    WALBERLA_CUDA_CLANG_TIDY "${SKIP_CLANG_TIDY_CHECKS}"
    "${SKIP_CLANG_TIDY_CHECKS_CUDA}")
endif()

# codegen-specific Clang-Tidy overrides
list(APPEND SKIP_CLANG_TIDY_CHECKS_CUDA
     "-bugprone-multi-level-implicit-pointer-conversion")
if(ESPRESSO_BUILD_WITH_CLANG_TIDY)
  espresso_override_clang_tidy_checks(
    WALBERLA_CXX_CLANG_TIDY_CODEGEN "${SKIP_CLANG_TIDY_CHECKS}"
    "${SKIP_CLANG_TIDY_CHECKS_CXX}")
  espresso_override_clang_tidy_checks(
    WALBERLA_CUDA_CLANG_TIDY_CODEGEN "${SKIP_CLANG_TIDY_CHECKS}"
    "${SKIP_CLANG_TIDY_CHECKS_CUDA}")
endif()

function(espresso_configure_walberla_target_clangtidy)
  set(TARGET_NAME ${ARGV0})
  if(ESPRESSO_BUILD_WITH_CLANG_TIDY)
    set(TARGET_LANG "CXX")
    if(TARGET_NAME MATCHES "_cuda")
      set(TARGET_LANG "CUDA")
    endif()
    set(TARGET_SUFFIX "")
    if(TARGET_NAME MATCHES "_codegen")
      set(TARGET_SUFFIX "_CODEGEN")
    endif()
    set_target_properties(
      ${TARGET_NAME}
      PROPERTIES ${TARGET_LANG}_CLANG_TIDY
                 "${WALBERLA_${TARGET_LANG}_CLANG_TIDY${TARGET_SUFFIX}}")
  endif()
endfunction()

function(espresso_configure_walberla_target)
  set(TARGET_NAME ${ARGV0})
  espresso_configure_walberla_target_clangtidy(${TARGET_NAME})
  target_link_libraries(${TARGET_NAME} PRIVATE ${WALBERLA_LIBS})
  target_include_directories(
    ${TARGET_NAME} PUBLIC include PRIVATE ${WALBERLA_INCLUDE_DIRS}
                                          ${walberla_BINARY_DIR}/src)
  get_target_property(TARGET_TYPE ${TARGET_NAME} TYPE)
  if(TARGET_TYPE STREQUAL "OBJECT_LIBRARY")
    set_target_properties(${TARGET_NAME} PROPERTIES POSITION_INDEPENDENT_CODE
                                                    ON)
  else()
    install(TARGETS ${TARGET_NAME}
            LIBRARY DESTINATION ${ESPRESSO_INSTALL_PYTHON}/espressomd)
  endif()
endfunction()

add_library(espresso_walberla_cpp_flags INTERFACE)
add_library(espresso::walberla::cpp_flags ALIAS espresso_walberla_cpp_flags)
target_link_libraries(
  espresso_walberla_cpp_flags
  INTERFACE
    espresso::cpp_flags
    $<$<BOOL:${ESPRESSO_BUILD_WITH_WALBERLA_AVX}>:espresso::avx_flags>
    $<$<BOOL:${ESPRESSO_BUILD_WITH_SHARED_MEMORY_PARALLELISM}>:OpenMP::OpenMP_CXX>
)
add_library(espresso_walberla_cuda_flags INTERFACE)
add_library(espresso::walberla::cuda_flags ALIAS espresso_walberla_cuda_flags)
target_link_libraries(
  espresso_walberla_cuda_flags
  INTERFACE espresso::cuda_flags
            $<$<BOOL:${ESPRESSO_BUILD_WITH_WALBERLA_AVX}>:espresso::avx_flags>)

add_library(espresso_walberla SHARED)
add_library(espresso_walberla_codegen OBJECT)
add_library(espresso::walberla ALIAS espresso_walberla)
add_library(espresso::walberla_codegen ALIAS espresso_walberla_codegen)

espresso_configure_walberla_target(espresso_walberla)
espresso_configure_walberla_target(espresso_walberla_codegen)

target_link_libraries(
  espresso_walberla PUBLIC MPI::MPI_CXX espresso::utils
  PRIVATE espresso::walberla::cpp_flags espresso::walberla_codegen)
target_link_libraries(espresso_walberla_codegen
                      PRIVATE espresso::walberla::cpp_flags)

if(WALBERLA_BUILD_WITH_CUDA)
  espresso_add_gpu_library(espresso_walberla_cuda SHARED)
  espresso_add_gpu_library(espresso_walberla_codegen_cuda OBJECT)
  add_library(espresso::walberla_cuda ALIAS espresso_walberla_cuda)
  add_library(espresso::walberla_codegen_cuda ALIAS
              espresso_walberla_codegen_cuda)
  espresso_configure_walberla_target(espresso_walberla_cuda)
  espresso_configure_walberla_target(espresso_walberla_codegen_cuda)
  target_link_libraries(
    espresso_walberla_cuda PUBLIC espresso::utils
    PRIVATE CUDA::cuda_driver CUDA::cudart espresso::walberla_codegen_cuda)
  target_link_libraries(espresso_walberla_codegen_cuda PRIVATE CUDA::cuda_driver
                                                               CUDA::cudart)
endif()

add_subdirectory(src)
if(ESPRESSO_BUILD_TESTS)
  add_subdirectory(tests)
endif()
