# Find a CUDA-capable compiler, include its libraries and declare a custom
# `add_library()` wrapper function named `add_gpu_library()`.

if(EXISTS "$ENV{NVCC}" AND NOT CUDA_NVCC_EXECUTABLE)
  set(CUDA_NVCC_EXECUTABLE $ENV{NVCC} CACHE FILEPATH "Path to CUDA compiler.")
endif()

# Option 1: the cuda compiler is the same as the C++ compiler (e.g. Clang)
if(CUDA_NVCC_EXECUTABLE STREQUAL CMAKE_CXX_COMPILER)
  message(STATUS "Found CUDA-capable host compiler: ${CUDA_NVCC_EXECUTABLE}")
  set(CUDA 1)
  set(CUDA_COMPILER_EXE ${CUDA_NVCC_EXECUTABLE})
  if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang"
     OR CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    execute_process(COMMAND ${CUDA_NVCC_EXECUTABLE} ${CMAKE_CXX_FLAGS}
                            --verbose
                    ERROR_VARIABLE CUDA_DIR_STRING)
    string(REGEX
           REPLACE "^.*Found CUDA installation: ([^,]+).*\$"
                   "\\1"
                   CUDA_DIR
                   "${CUDA_DIR_STRING}")
    string(REGEX
           REPLACE "^.*Found CUDA installation: .* version ([0-9\.]+|unknown).*\$"
                   "\\1"
                   CUDA_VERSION
                   "${CUDA_DIR_STRING}")
    if(NOT CUDA_DIR_STRING MATCHES "Found CUDA installation" OR CUDA_VERSION STREQUAL "unknown")
      message(FATAL_ERROR "Clang found no compatible CUDA library.")
    endif()
    message(STATUS "Found CUDA version: ${CUDA_VERSION}")
    message(STATUS "Found CUDA installation: ${CUDA_DIR}")
    if(CUDA_VERSION VERSION_LESS 7.0)
      message(FATAL_ERROR "CUDA version does not match requirements.")
    endif()
  else()
      set(CUDA_DIR "/usr/local/cuda")
  endif()
  find_library(CUDART_LIBRARY NAMES cudart PATHS ${CUDA_DIR}/lib64 ${CUDA_DIR}/lib /usr/local/nvidia/lib /usr/lib/x86_64-linux-gnu NO_DEFAULT_PATH)
  find_library(CUFFT_LIBRARY NAMES cufft PATHS ${CUDA_DIR}/lib64 ${CUDA_DIR}/lib /usr/local/nvidia/lib /usr/lib/x86_64-linux-gnu NO_DEFAULT_PATH)

  function(add_gpu_library)
    set(options STATIC SHARED MODULE EXCLUDE_FROM_ALL)
    set(oneValueArgs)
    set(multiValueArgs)
    cmake_parse_arguments(ARG "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )
    list(REMOVE_AT ARG_UNPARSED_ARGUMENTS 0)
    set_source_files_properties(${ARG_UNPARSED_ARGUMENTS} PROPERTIES LANGUAGE "CXX" COMPILE_FLAGS "${CUDA_NVCC_FLAGS}")
    add_library(${ARGV})
    set_target_properties(${ARGV0} PROPERTIES LINKER_LANGUAGE "CXX")
    target_link_libraries(${ARGV0} PRIVATE ${CUDA_LIBRARY} ${CUDART_LIBRARY})
    target_link_libraries(${ARGV0} PRIVATE ${CUFFT_LIBRARY})

    foreach(file ${ARG_UNPARSED_ARGUMENTS})
      if(${file} MATCHES "\\.cu$")
              if (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 3.8.9)
                set_source_files_properties(${file} PROPERTIES COMPILE_FLAGS "--cuda-gpu-arch=sm_30 --cuda-gpu-arch=sm_52")
              else()
                set_source_files_properties(${file} PROPERTIES COMPILE_FLAGS "--cuda-gpu-arch=sm_30")
              endif()
      endif()
    endforeach()
  endfunction()
else()
  # Option 2: the cuda compiler is HIP
  set(ROCM_HOME "/opt/rocm" CACHE STRING "Path to AMD ROCm")
  list(APPEND CMAKE_MODULE_PATH "${ROCM_HOME}/hip/cmake")
  find_package(HIP 1.5.18494 QUIET MODULE)
  if(HIP_FOUND)
    if(HIP_VERSION VERSION_LESS "3.1")
      set(HCC_PATH "${HIP_ROOT_DIR}")
    else()
      set(HCC_PATH "${ROCM_HOME}/hcc")
    endif()

    find_package(HIP MODULE)
    message(STATUS "Found HIP compiler: ${HIP_HIPCC_EXECUTABLE}")
    set(CUDA 1)
    set(HIP 1)
    set(CUDA_COMPILER_EXE ${HIP_HIPCC_EXECUTABLE})
    list(APPEND HIP_HCC_FLAGS "-I${HIP_ROOT_DIR}/include -I${ROCM_HOME}/include -Wno-c99-designator -Wno-macro-redefined -Wno-duplicate-decl-specifier -std=c++${CMAKE_CXX_STANDARD}")
    list(APPEND HIP_HCC_FLAGS "-pedantic -Wall -Wextra -Wno-sign-compare -Wno-unused-function -Wno-unused-variable -Wno-unused-parameter -Wno-missing-braces -Wno-gnu-anonymous-struct -Wno-nested-anon-types -Wno-gnu-zero-variadic-macro-arguments")
    if(NOT HIP_VERSION VERSION_LESS "3.3")
      list(APPEND HIP_HCC_FLAGS "-Wno-deprecated-copy")
    endif()
    if(WARNINGS_ARE_ERRORS)
      list(APPEND HIP_HCC_FLAGS "-Werror")
    endif()

    find_library(ROCFFT_LIB name "rocfft" PATHS "${ROCM_HOME}/lib")

    function(add_gpu_library)
      hip_add_library(${ARGV})
      set_target_properties(${ARGV0} PROPERTIES LINKER_LANGUAGE HIP)
      target_link_libraries(${ARGV0} PRIVATE "${ROCFFT_LIB}")
    endfunction()
  else()
    # Option 3: the cuda compiler is NVCC
    find_package(CUDA 9.0)
    if(CUDA_FOUND)
      if(NOT CUDA_NVCC_EXECUTABLE STREQUAL "${CUDA_TOOLKIT_ROOT_DIR}/bin/nvcc")
        get_filename_component(NVCC_EXECUTABLE_DIRNAME "${CUDA_NVCC_EXECUTABLE}" DIRECTORY)
        get_filename_component(NVCC_EXECUTABLE_DIRNAME "${NVCC_EXECUTABLE_DIRNAME}" DIRECTORY)
        message(WARNING "Your nvcc (${CUDA_NVCC_EXECUTABLE}) does not appear to match your CUDA libraries (in ${CUDA_TOOLKIT_ROOT_DIR}). While Espresso will still compile, you might get unexpected crashes. Please point CUDA_TOOLKIT_ROOT_DIR to your CUDA toolkit path, e.g. by adding -DCUDA_TOOLKIT_ROOT_DIR=${NVCC_EXECUTABLE_DIRNAME} to your cmake command.")
      endif()
      set(CUDA 1)
      set(CUDA_COMPILER_EXE ${CUDA_NVCC_EXECUTABLE})

      set(CUDA_LINK_LIBRARIES_KEYWORD PUBLIC)

      set(CUDA_NVCC_FLAGS_DEBUG "${CUDA_NVCC_FLAGS_DEBUG} -g")
      set(CUDA_NVCC_FLAGS_RELEASE "${CUDA_NVCC_FLAGS_RELEASE} -O3 -DNDEBUG")
      set(CUDA_NVCC_FLAGS_MINSIZEREL "${CUDA_NVCC_FLAGS_MINSIZEREL} -Os -DNDEBUG")
      set(CUDA_NVCC_FLAGS_RELWITHDEBINFO "${CUDA_NVCC_FLAGS_RELWITHDEBINFO} -g -O2")
      set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -gencode=arch=compute_30,code=sm_30 -gencode=arch=compute_52,code=sm_52 -gencode=arch=compute_52,code=compute_52 -std=c++${CMAKE_CXX_STANDARD}")
      SET(CUDA_PROPAGATE_HOST_FLAGS OFF)

      if(WARNINGS_ARE_ERRORS)
        set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -Xcompiler -Werror")
      endif()

      if (CMAKE_OSX_SYSROOT)
        set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -Xcompiler -isysroot -Xcompiler ${CMAKE_OSX_SYSROOT}")
      endif()

      function(add_gpu_library)
        cuda_add_library(${ARGV})
        set_property(TARGET ${ARGV0} PROPERTY CUDA_SEPARABLE_COMPILATION ON)
        target_link_libraries(${ARGV0} PRIVATE ${CUDA_CUFFT_LIBRARIES})
      endfunction()

    endif(CUDA_FOUND)
  endif()
endif()

include( FindPackageHandleStandardArgs )
FIND_PACKAGE_HANDLE_STANDARD_ARGS( CudaCompiler REQUIRED_VARS CUDA_COMPILER_EXE )
