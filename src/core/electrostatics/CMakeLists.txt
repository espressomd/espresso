file(GLOB Electrostatics_SRC *.cpp)
add_library(Electrostatics SHARED ${Electrostatics_SRC})

set_target_properties(Electrostatics PROPERTIES LINKER_LANGUAGE CXX)

if(SCAFACOS)
  add_dependencies(Electrostatics Scafacos)
  target_link_libraries(Electrostatics Scafacos)

  add_subdirectory(scafacos)
endif(SCAFACOS)
