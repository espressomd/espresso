add_library(Electrostatics INTERFACE)

if(SCAFACOS)
  target_link_libraries(Electrostatics Scafacos)

  add_subdirectory(scafacos)
endif(SCAFACOS)
