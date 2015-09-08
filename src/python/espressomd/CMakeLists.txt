file(GLOB cython_SRC *.pyx)
file(GLOB cython_HEADER *.pxd)

add_custom_target(gen_pxiconfig_cpp)

add_custom_command(TARGET gen_pxiconfig_cpp
                   COMMAND ${PYTHON_EXECUTABLE}
                   ${CMAKE_CURRENT_SOURCE_DIR}/gen_pxiconfig.py
                   ${CMAKE_SOURCE_DIR}/src/features.def ${CMAKE_CURRENT_BINARY_DIR}/gen_pxiconfig.cpp
                   DEPENDS ${CMAKE_SOURCE_DIR}/src/features.def
                  )

include_directories(${CMAKE_BINARY_DIR}/src)
include_directories(${CMAKE_SOURCE_DIR}/src/core)
include_directories(${CMAKE_BINARY_DIR}/src/core)

add_executable(gen_pxiconfig gen_pxiconfig.cpp)
add_dependencies(gen_pxiconfig gen_pxi_config_cpp)

add_custom_target(myconfig-pxi)
add_dependencies(myconfig-pxi gen_pxiconfig)

add_custom_command(TARGET myconfig-pxi
                   COMMAND ${CMAKE_CURRENT_BINARY_DIR}/gen_pxiconfig > ${CMAKE_CURRENT_BINARY_DIR}/_espresso/myconfig.pxi
                  )

add_custom_target(cython-generated-files)
add_dependencies(cython-generated-files myconfig-pxi)

foreach(cython_file ${cython_SRC})
  get_filename_component(basename ${cython_file} NAME_WE)
  list(APPEND cython_generated_SRC "${basename}.cpp")
  add_custom_command(TARGET cython-generated-files
                     COMMAND ${CYTHON_EXECUTABLE}  --cplus
                     -I ${CMAKE_CURRENT_BINARY_DIR}/_espresso
                     -I ${CMAKE_CURRENT_SOURCE_DIR}
                     ${cython_file}
                     -o ${CMAKE_CURRENT_BINARY_DIR}/${basename}.cpp
                    )
endforeach()
