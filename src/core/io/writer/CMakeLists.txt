target_sources(
  EspressoCore
  PRIVATE "$<$<BOOL:${H5MD}>:${CMAKE_CURRENT_SOURCE_DIR}/h5md_core.cpp>")
