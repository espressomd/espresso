target_sources(
  ScriptInterface
  PRIVATE $<$<BOOL:H5MD>:${CMAKE_CURRENT_SOURCE_DIR}/initialize.cpp>
          $<$<BOOL:H5MD>:${CMAKE_CURRENT_SOURCE_DIR}/h5md.cpp>)
