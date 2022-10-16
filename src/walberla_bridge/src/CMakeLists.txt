add_subdirectory(electrokinetics)
add_subdirectory(lattice_boltzmann)

target_sources(
        espresso_walberla
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/LatticeModel.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/LatticeWalberla.cpp)
