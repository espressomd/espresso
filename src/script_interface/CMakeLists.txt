file(GLOB_RECURSE EspressoScriptInterface_SRC
          "*.cpp"
          )

add_library(EspressoScriptInterface SHARED ${EspressoScriptInterface_SRC})
target_link_libraries(EspressoScriptInterface EspressoCore)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

include_directories(${CMAKE_SOURCE_DIR}/src)
include_directories(${CMAKE_SOURCE_DIR}/src/core)
include_directories(${CMAKE_BINARY_DIR}/src/core)

