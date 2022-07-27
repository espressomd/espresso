target_sources(
  Espresso_script_interface
  PRIVATE $<$<BOOL:H5MD>:${CMAKE_CURRENT_SOURCE_DIR}/initialize.cpp>
          $<$<BOOL:H5MD>:${CMAKE_CURRENT_SOURCE_DIR}/h5md.cpp>)
