file(GLOB Correlators_SRC *.?pp)
add_library(Correlators SHARED ${Correlators_SRC})
add_dependencies(Correlators EspressoConfig)
set_target_properties(Correlators PROPERTIES MACOSX_RPATH TRUE)
