file(GLOB ImmersedBoundary_SRC *.cpp)
add_library(ImmersedBoundary SHARED ${ImmersedBoundary_SRC})
install(TARGETS ImmersedBoundary LIBRARY DESTINATION ${LIBDIR} ARCHIVE DESTINATION ${LIBDIR})
set_target_properties(ImmersedBoundary PROPERTIES SOVERSION ${SOVERSION})
add_dependencies(ImmersedBoundary EspressoConfig)


set_target_properties(ImmersedBoundary PROPERTIES LINKER_LANGUAGE CXX)
