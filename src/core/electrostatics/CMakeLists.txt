
file(GLOB Electrostatics_SRC *.cpp)
add_library(Electrostatics SHARED ${Electrostatics_SRC})
set_target_properties(Electrostatics PROPERTIES MACOSX_RPATH TRUE)
set_target_properties(Electrostatics PROPERTIES LINKER_LANGUAGE CXX)

if(WITH_SCAFACOS)
  add_dependencies(Electrostatics Scafacos)
  set_target_properties(Electrostatics PROPERTIES MACOSX_RPATH TRUE)
  target_link_libraries(Electrostatics Scafacos)

  add_subdirectory(scafacos)
endif(WITH_SCAFACOS)
