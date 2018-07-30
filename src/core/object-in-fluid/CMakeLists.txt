file(GLOB ObjectInFluid_SRC *.cpp)
add_library(ObjectInFluid SHARED ${ObjectInFluid_SRC})
install(TARGETS ObjectInFluid LIBRARY DESTINATION ${LIBDIR} ARCHIVE DESTINATION ${LIBDIR})
set_target_properties(ObjectInFluid PROPERTIES SOVERSION ${SOVERSION})
add_dependencies(ObjectInFluid EspressoConfig)

