if(H5MD)
  target_sources(
    Espresso_core PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/h5md_core.cpp
                          ${CMAKE_CURRENT_SOURCE_DIR}/h5md_specification.cpp)
endif(H5MD)
