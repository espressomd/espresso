set(EspressoScriptInterface_SRC
    "${CMAKE_CURRENT_SOURCE_DIR}/initialize.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/constraints/initialize.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/correlators/initialize.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/lbboundaries/initialize.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/observables/initialize.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/ScriptInterfaceBase.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/shapes/Cylinder.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/shapes/HollowCone.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/shapes/initialize.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/shapes/Maze.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/shapes/Pore.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/shapes/Rhomboid.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/shapes/Slitpore.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/shapes/Sphere.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/shapes/SpheroCylinder.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/shapes/Stomatocyte.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/shapes/Wall.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/Variant.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/ParallelScriptInterface.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/ParallelScriptInterfaceSlave.cpp"
)

if(H5MD)
    list(APPEND EspressoScriptInterface_SRC 
        "${CMAKE_CURRENT_SOURCE_DIR}/h5md/initialize.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/h5md/h5md.cpp"
        )
endif(H5MD)
add_library(EspressoScriptInterface SHARED ${EspressoScriptInterface_SRC})
add_dependencies(EspressoScriptInterface EspressoConfig EspressoConfig)
set_target_properties(EspressoScriptInterface PROPERTIES MACOSX_RPATH TRUE)
target_link_libraries(EspressoScriptInterface EspressoCore)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

include_directories(${CMAKE_SOURCE_DIR}/src)
include_directories(${CMAKE_SOURCE_DIR}/src/core)
include_directories(${CMAKE_BINARY_DIR}/src/core)
if (HDF5_FOUND)
    include_directories(${CMAKE_SOURCE_DIR}/libs/h5xx)
    include_directories(${CMAKE_SOURCE_DIR}/src/core/io/writer/h5md)
    include_directories(${CMAKE_BINARY_DIR}/src/core/io/writer/h5md)
endif(HDF5_FOUND)
