file(GLOB EspressoCore_SRC
          "*.cpp" "io/mpiio/mpiio.cpp")


if( WITH_COVERAGE )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Og --coverage -fprofile-arcs -ftest-coverage")
  link_libraries(gcov)
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_BINARY_DIR})

add_custom_command(OUTPUT ${CMAKE_BINARY_DIR}/src/core/config-features.hpp ${CMAKE_BINARY_DIR}/src/core/config-features.cpp
  COMMAND
  ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/src/core/gen_featureconfig.py
  ${CMAKE_SOURCE_DIR}/src/features.def
  ${CMAKE_BINARY_DIR}/src/core/config-features.hpp
  ${CMAKE_BINARY_DIR}/src/core/config-features.cpp
  DEPENDS ${CMAKE_SOURCE_DIR}/src/features.def ${CMAKE_SOURCE_DIR}/src/core/gen_featureconfig.py
  )

add_custom_command(OUTPUT config-version.cpp
                   COMMAND ${CMAKE_SOURCE_DIR}/config/genversion.sh -c > ${CMAKE_BINARY_DIR}/src/core/config-version.cpp
)

# Custom target for everything that needs to be triggered when the myconfig changes
add_library(EspressoConfig SHARED config-features.cpp)
add_dependencies(EspressoConfig myconfig check_myconfig)
install(TARGETS EspressoConfig LIBRARY DESTINATION ${LIBDIR})
set_target_properties(EspressoConfig PROPERTIES SOVERSION ${SOVERSION})

add_library(EspressoCore SHARED ${EspressoCore_SRC} config-version.cpp)
install(TARGETS EspressoCore LIBRARY DESTINATION ${LIBDIR})
set_target_properties(EspressoCore PROPERTIES SOVERSION ${SOVERSION})

target_link_libraries(EspressoCore EspressoConfig ${LIBRARIES} Actor ObjectInFluid ImmersedBoundary Shapes Constraints Observables Accumulators cluster_analysis VirtualSites)

if(SCAFACOS)
  target_link_libraries(EspressoCore Scafacos)
endif(SCAFACOS)
if(GSL)
  include_directories(${GSL_INCLUDE_DIRS})
endif(GSL)


if(GSL)
  target_link_libraries(EspressoCore gsl gslcblas)
endif(GSL)

if(GSL)
  target_link_libraries(EspressoCore gsl gslcblas)
endif(GSL)

if(CUDA)
  file(GLOB EspressoCuda_SRC
          "*.cu"
          )

  cuda_include_directories(${CMAKE_CURRENT_SOURCE_DIR})
  cuda_include_directories(${CMAKE_CURRENT_BINARY_DIR})

  cuda_add_library(EspressoCuda SHARED ${EspressoCuda_SRC})
  install(TARGETS EspressoCuda DESTINATION ${LIBDIR})
  set_target_properties(EspressoCuda PROPERTIES SOVERSION ${SOVERSION})
  add_dependencies(EspressoCuda EspressoConfig)

  CUDA_ADD_CUFFT_TO_TARGET(EspressoCuda)
  target_link_libraries(EspressoCore EspressoCuda)
endif(CUDA)

# Subdirectories
add_subdirectory(observables)
add_subdirectory(accumulators)
add_subdirectory(io)
if(HDF5_FOUND)
    target_link_libraries(EspressoCore H5mdCore)
endif(HDF5_FOUND)

if(WITH_UNIT_TESTS)
  add_subdirectory(unit_tests)
endif(WITH_UNIT_TESTS)

add_subdirectory(actor)
add_subdirectory(immersed_boundary)
add_subdirectory(virtual_sites)
add_subdirectory(lbboundaries)
add_subdirectory(object-in-fluid)
add_subdirectory(shapes)
add_subdirectory(constraints)
add_subdirectory(cluster_analysis)

if(SCAFACOS)
  add_subdirectory(scafacos)
endif(SCAFACOS)
