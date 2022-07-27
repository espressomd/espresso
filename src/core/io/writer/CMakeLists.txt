if(H5MD)
  target_link_libraries(Espresso_core PUBLIC ${HDF5_LIBRARIES}
                                             Boost::filesystem h5xx)
  target_include_directories(Espresso_core PUBLIC ${CMAKE_CURRRENT_SOURCE_DIR}
                                                  ${HDF5_INCLUDE_DIRS})

  target_compile_definitions(Espresso_core PUBLIC H5XX_USE_MPI)

  target_sources(
    Espresso_core PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/h5md_core.cpp
                          ${CMAKE_CURRENT_SOURCE_DIR}/h5md_specification.cpp)
endif(H5MD)
