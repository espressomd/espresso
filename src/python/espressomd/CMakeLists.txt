file(GLOB cython_SRC *.pyx)
file(GLOB cython_HEADER *.pxd)
file(GLOB cython_AUX *.py)
file(COPY ${cython_AUX} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

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

add_custom_command(OUTPUT _espresso/myconfig.pxi
		   COMMAND mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/_espresso
                   COMMAND ${CMAKE_CURRENT_BINARY_DIR}/gen_pxiconfig > ${CMAKE_CURRENT_BINARY_DIR}/_espresso/myconfig.pxi
                   DEPENDS gen_pxiconfig
                   )

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

  if (${basename} STREQUAL "_tcl")
    file(GLOB_RECURSE Espresso_tcl_SRC
            "${CMAKE_CURRENT_SOURCE_DIR}/../../tcl/*.cpp"
            )
    add_library(${basename} SHARED ${basename}.cpp ${Espresso_tcl_SRC})
    target_include_directories(${basename} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../../tcl)
    add_definitions(-DESPRESSO_SCRIPTS_DEFAULT="${CMAKE_BINARY_DIR}/scripts")
  else()
    add_library(${basename} SHARED ${basename}.cpp)
  endif()
  set_target_properties(${basename} PROPERTIES PREFIX "")
  set_target_properties(${basename} PROPERTIES MACOSX_RPATH TRUE)
  if (APPLE)
     set_target_properties(${basename} PROPERTIES SUFFIX ".so")
  endif()
  target_link_libraries(${basename} EspressoCore ${LIBRARIES})
endforeach()

