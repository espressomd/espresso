file(GLOB cython_SRC *.pyx)
file(GLOB cython_HEADER *.pxd)
file(GLOB cython_AUX *.py)
file(COPY ${cython_AUX} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
#add_custom_target(gen_pxiconfig_cpp)

add_custom_command(OUTPUT gen_pxiconfig.cpp
                   COMMAND ${PYTHON_EXECUTABLE}
                   ${CMAKE_CURRENT_SOURCE_DIR}/gen_pxiconfig.py
                   ${CMAKE_SOURCE_DIR}/src/features.def ${CMAKE_CURRENT_BINARY_DIR}/gen_pxiconfig.cpp
                   DEPENDS ${CMAKE_SOURCE_DIR}/src/features.def
                  )

include_directories(${CMAKE_BINARY_DIR}/src)
include_directories(${CMAKE_SOURCE_DIR}/src/core)
include_directories(${CMAKE_BINARY_DIR}/src/core)

add_executable(gen_pxiconfig gen_pxiconfig.cpp)
#add_dependencies(gen_pxiconfig gen_pxiconfig.cpp)

#dd_custom_target(myconfig-pxi)
#add_dependencies(myconfig-pxi gen_pxiconfig)

add_custom_command(OUTPUT _espresso/myconfig.pxi
		   COMMAND mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/_espresso
                   COMMAND ${CMAKE_CURRENT_BINARY_DIR}/gen_pxiconfig > ${CMAKE_CURRENT_BINARY_DIR}/_espresso/myconfig.pxi
                   DEPENDS gen_pxiconfig
                   )

#add_custom_target(cython-generated-files)
#add_dependencies(cython-generated-files myconfig-pxi)
include_directories(${PYTHON_INCLUDE_DIRS})
include_directories(${NUMPY_INCLUDE_DIR})

foreach(cython_file ${cython_SRC})
  get_filename_component(basename ${cython_file} NAME_WE)
  list(APPEND cython_generated_SRC "${basename}.cpp")
  add_custom_command(OUTPUT ${basename}.cpp
                     COMMAND ${CYTHON_EXECUTABLE}  --cplus
                     -I ${CMAKE_CURRENT_BINARY_DIR}/_espresso
                     -I ${CMAKE_CURRENT_SOURCE_DIR}
                     -I ${CMAKE_CURRENT_SOURCE_DIR}/_espresso
                     ${cython_file}
                     -o ${CMAKE_CURRENT_BINARY_DIR}/${basename}.cpp
                     DEPENDS _espresso/myconfig.pxi 
                     )

  add_library(${basename} SHARED ${basename}.cpp)
  set_target_properties(${basename} PROPERTIES PREFIX "")
  target_link_libraries(${basename} EspressoCore ${LIBRARIES})
endforeach()

#include_directories(${PYTHON_INCLUDE_DIRS})
#include_directories(${NUMPY_INCLUDE_DIR})

#add_library(EspressoCython SHARED ${cython_generated_SRC})
#add_dependencies(EspressoCython myconfig-pxi)
#target_link_libraries(EspressoCython EspressoCore ${LIBRARIES})
