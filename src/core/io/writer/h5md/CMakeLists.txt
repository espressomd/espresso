add_library(H5mdCore SHARED h5md_core.cpp)
add_dependencies(H5mdCore EspressoConfig)
set_target_properties(H5mdCore PROPERTIES MACOSX_RPATH TRUE)
