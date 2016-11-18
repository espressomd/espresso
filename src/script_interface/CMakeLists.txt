file(GLOB_RECURSE EspressoScriptInterface_SRC
          "*.?pp"
          )
      #set(EspressoScriptInterface_SRC
      #    "${CMAKE_CURRENT_SOURCE_DIR}/ParallelScriptInterface.hpp"
      #    "${CMAKE_CURRENT_SOURCE_DIR}/ParallelScriptInterfaceSlave.hpp"
      #    "${CMAKE_CURRENT_SOURCE_DIR}/Parameter.hpp"
      #    "${CMAKE_CURRENT_SOURCE_DIR}/ScriptInterface.hpp"
      #    "${CMAKE_CURRENT_SOURCE_DIR}/ScriptInterfaceBase.cpp"
      #    "${CMAKE_CURRENT_SOURCE_DIR}/ScriptInterfaceBase.hpp"
      #    "${CMAKE_CURRENT_SOURCE_DIR}/ScriptObjectRegistry.hpp"
      #    "${CMAKE_CURRENT_SOURCE_DIR}/initialize.cpp"
      #    "${CMAKE_CURRENT_SOURCE_DIR}/initialize.hpp"
      #    "${CMAKE_CURRENT_SOURCE_DIR}/constraints/Constraint.hpp"
      #    "${CMAKE_CURRENT_SOURCE_DIR}/constraints/Constraints.hpp"
      #    "${CMAKE_CURRENT_SOURCE_DIR}/constraints/initialize.hpp"
      #    "${CMAKE_CURRENT_SOURCE_DIR}/constraints/initialize.cpp"
      #    "${CMAKE_CURRENT_SOURCE_DIR}/correlators/AutoUpdateCorrelators.hpp"
      #    "${CMAKE_CURRENT_SOURCE_DIR}/correlators/Correlator.hpp"
      #    "${CMAKE_CURRENT_SOURCE_DIR}/correlators/initialize.hpp"
      #    "${CMAKE_CURRENT_SOURCE_DIR}/correlators/initialize.cpp"
      #    "${CMAKE_CURRENT_SOURCE_DIR}/lbboundaries/LBBoundaries.hpp"
      #    "${CMAKE_CURRENT_SOURCE_DIR}/lbboundaries/LBBoundary.hpp"
      #    "${CMAKE_CURRENT_SOURCE_DIR}/lbboundaries/initialize.hpp"
      #    "${CMAKE_CURRENT_SOURCE_DIR}/lbboundaries/initialize.cpp"
      #    "${CMAKE_CURRENT_SOURCE_DIR}/observables/initialize.hpp"
      #    "${CMAKE_CURRENT_SOURCE_DIR}/observables/initialize.cpp"
      #    "${CMAKE_CURRENT_SOURCE_DIR}/observables/AutoUpdateObservables.hpp"
      #    "${CMAKE_CURRENT_SOURCE_DIR}/observables/PidObservable.hpp"
      #    "${CMAKE_CURRENT_SOURCE_DIR}/observables/Observable.hpp"
      #    "${CMAKE_CURRENT_SOURCE_DIR}/observables/ParamlessObservable.hpp"
      #    "${CMAKE_CURRENT_SOURCE_DIR}/observables/ProfileObservable.hpp"
      #    "${CMAKE_CURRENT_SOURCE_DIR}/observables/ProfileOb"
      #    "${CMAKE_CURRENT_SOURCE_DIR}/shapes/Cylinder.cpp"
      #    "${CMAKE_CURRENT_SOURCE_DIR}/shapes/Cylinder.hpp"
      #    "${CMAKE_CURRENT_SOURCE_DIR}/shapes/HollowCone.cpp"
      #    "${CMAKE_CURRENT_SOURCE_DIR}/shapes/HollowCone.hpp"
      #    "${CMAKE_CURRENT_SOURCE_DIR}/shapes/Maze.cpp"
      #    "${CMAKE_CURRENT_SOURCE_DIR}/shapes/Maze.hpp"
      #    "${CMAKE_CURRENT_SOURCE_DIR}/shapes/NoWhere.hpp"
      #    "${CMAKE_CURRENT_SOURCE_DIR}/shapes/Pore.cpp"
      #    "${CMAKE_CURRENT_SOURCE_DIR}/shapes/Pore.hpp"
      #    "${CMAKE_CURRENT_SOURCE_DIR}/shapes/Rhomboid.cpp"
      #    "${CMAKE_CURRENT_SOURCE_DIR}/shapes/Rhomboid.hpp"
      #    "${CMAKE_CURRENT_SOURCE_DIR}/shapes/Shape.hpp"
      #    "${CMAKE_CURRENT_SOURCE_DIR}/shapes/Slitpore.cpp"
      #    "${CMAKE_CURRENT_SOURCE_DIR}/shapes/Slitpore.hpp"
      #    "${CMAKE_CURRENT_SOURCE_DIR}/shapes/Sphere.cpp"
      #    "${CMAKE_CURRENT_SOURCE_DIR}/shapes/Sphere.hpp"
      #    "${CMAKE_CURRENT_SOURCE_DIR}/shapes/SpheroCylinder.cpp"
      #    "${CMAKE_CURRENT_SOURCE_DIR}/shapes/SpheroCylinder.hpp"
      #    "${CMAKE_CURRENT_SOURCE_DIR}/shapes/Stomatocyte.cpp"
      #    "${CMAKE_CURRENT_SOURCE_DIR}/shapes/Stomatocyte.hpp"
      #    "${CMAKE_CURRENT_SOURCE_DIR}/shapes/Wall.cpp"
      #    "${CMAKE_CURRENT_SOURCE_DIR}/shapes/Wall.hpp"
      #    "${CMAKE_CURRENT_SOURCE_DIR}/shapes/initialize.cpp"
      #    "${CMAKE_CURRENT_SOURCE_DIR}/shapes/initialize.hpp"
      #)

      #if(H5MD)
      #    list(APPEND EspressoScriptInterface_SRC 
      #        "${CMAKE_CURRENT_SOURCE_DIR}/h5md/initialize.hpp"
      #        "${CMAKE_CURRENT_SOURCE_DIR}/h5md/initialize.cpp"
      #        "${CMAKE_CURRENT_SOURCE_DIR}/h5md/h5md.hpp"
      #        "${CMAKE_CURRENT_SOURCE_DIR}/h5md/h5md.cpp"
      #        )
      #endif(H5MD)
message("${EspressoScriptInterface_SRC}")
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
