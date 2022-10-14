add_subdirectory(lattice_boltzmann)
if(WITH_WALBERLA_ELECTROKINETICS)
  add_subdirectory(electrokinetics)
endif()

target_sources(
        espresso_walberla
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/LatticeModel.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/LatticeWalberla.cpp)
