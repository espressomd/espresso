file(GLOB ObjectInFluid_SRC *.cpp)
add_library(ObjectInFluid SHARED ${ObjectInFluid_SRC})
add_dependencies(ObjectInFluid EspressoConfig)
set_target_properties(ObjectInFluid PROPERTIES MACOSX_RPATH TRUE)
