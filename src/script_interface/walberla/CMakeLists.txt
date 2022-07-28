if(WITH_WALBERLA)
  target_sources(
    Espresso_script_interface
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/initialize.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/FluidWalberla.cpp)
  target_link_libraries(Espresso_script_interface PRIVATE Espresso::walberla)
endif(WITH_WALBERLA)
