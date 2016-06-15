file(GLOB_RECURSE EspressoScriptInterface_SRC
          "*.cpp"
          )

add_library(EspressoScriptInterface SHARED ${EspressoScriptInterface_SRC})

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../core)
