if(WITH_WALBERLA)
  target_sources(
    espresso_script_interface
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/initialize.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/FluidWalberla.cpp)
  target_link_libraries(espresso_script_interface PRIVATE espresso::walberla)
endif(WITH_WALBERLA)
