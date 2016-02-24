add_library(Electrostatics INTERFACE)

if(SCAFACOS)
  add_dependencies(Electrostatics Scafacos)
  target_link_libraries(Electrostatics Scafacos)

  add_subdirectory(scafacos)
endif(SCAFACOS)
