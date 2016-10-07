file(GLOB EspressoUtils_SRC *.cpp)
add_library(EspressoUtils SHARED ${EspressoUtils_SRC})
set_target_properties(EspressoUtils PROPERTIES MACOSX_RPATH TRUE)
