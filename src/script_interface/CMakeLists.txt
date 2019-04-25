set(EspressoScriptInterface_SRC
    "${CMAKE_CURRENT_SOURCE_DIR}/accumulators/initialize.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/initialize.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/collision_detection/initialize.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/constraints/initialize.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/cluster_analysis/initialize.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/lbboundaries/initialize.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/virtual_sites/initialize.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/observables/initialize.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/pair_criteria/initialize.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/mpiio/initialize.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/ScriptInterfaceBase.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/shapes/initialize.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/ParallelScriptInterface.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/ScriptInterfaceBase.cpp"
)

if(H5MD)
    list(APPEND EspressoScriptInterface_SRC
        "${CMAKE_CURRENT_SOURCE_DIR}/h5md/initialize.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/h5md/h5md.cpp"
        )
endif(H5MD)

if( WITH_COVERAGE )
  if( CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-instr-generate -fcoverage-mapping")
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage -fprofile-arcs -ftest-coverage")
    link_libraries(gcov)
  endif()
endif()

add_library(EspressoScriptInterface SHARED ${EspressoScriptInterface_SRC})
install(TARGETS EspressoScriptInterface LIBRARY DESTINATION ${PYTHON_INSTDIR}/espressomd)
set_target_properties(EspressoScriptInterface PROPERTIES SOVERSION ${SOVERSION})
target_link_libraries(EspressoScriptInterface PRIVATE EspressoConfig PUBLIC EspressoCore Shapes LBBoundaries cluster_analysis mpiio pdbreader)
if(HDF5_FOUND)
 target_link_libraries(EspressoScriptInterface PRIVATE H5mdCore)
endif()
target_include_directories(EspressoScriptInterface PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_include_directories(EspressoScriptInterface PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
