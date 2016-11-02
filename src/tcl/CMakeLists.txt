file(GLOB_RECURSE Espresso_SRC
          "*.cpp"
)

list(REMOVE_ITEM Espresso_SRC ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)
list(REMOVE_ITEM Espresso_SRC ${CMAKE_CURRENT_SOURCE_DIR}/scriptsdir.cpp)
#list(REMOVE_ITEM Espresso_SRC ${CMAKE_CURRENT_SOURCE_DIR}/initialize_interpreter.cpp)

include_directories(${CMAKE_SOURCE_DIR}/src/core)
include_directories(${CMAKE_BINARY_DIR}/src/core)
include_directories(${CMAKE_SOURCE_DIR}/src/tcl)

file(GLOB espresso_scripts ${CMAKE_SOURCE_DIR}/scripts/*tcl)
file(COPY ${espresso_scripts} DESTINATION ${CMAKE_BINARY_DIR}/scripts)
install(DIRECTORY ${CMAKE_BINARY_DIR}/scripts DESTINATION ${DATA})

add_library(EspressoTcl SHARED ${Espresso_SRC})
set_target_properties(EspressoTcl PROPERTIES SOVERSION ${SOVERSION})
set_target_properties(EspressoTcl PROPERTIES MACOSX_RPATH TRUE)
add_dependencies(EspressoTcl EspressoConfig)
install(TARGETS EspressoTcl DESTINATION ${LIBDIR})

add_executable(Espresso main.cpp scriptsdir.cpp ${MPI_SOURCE})
set_target_properties(Espresso PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
target_compile_definitions(Espresso PRIVATE
  -DESPRESSO_SCRIPTS_DEFAULT="${CMAKE_BINARY_DIR}/scripts")
add_dependencies(Espresso EspressoConfig)
target_link_libraries(Espresso EspressoTcl EspressoCore ${LIBRARIES})

add_executable(Espresso_install main.cpp scriptsdir.cpp ${MPI_SOURCE})
set_target_properties(Espresso_install PROPERTIES
  OUTPUT_NAME Espresso)
target_compile_definitions(Espresso_install PRIVATE
  -DESPRESSO_SCRIPTS_DEFAULT="${CMAKE_INSTALL_PREFIX}/${DATA}/scripts")
add_dependencies(Espresso_install EspressoConfig)
target_link_libraries(Espresso_install EspressoTcl EspressoCore ${LIBRARIES})
install(TARGETS Espresso_install DESTINATION ${BINDIR})
