set(Accumulators_SRC
    Correlator.hpp
    Correlator.cpp
    MeanVarianceCalculator.hpp
    MeanVarianceCalculator.cpp
)
add_library(Accumulators SHARED ${Accumulators_SRC})
install(TARGETS Accumulators LIBRARY DESTINATION ${PYTHON_INSTDIR}/espressomd)
set_target_properties(Accumulators PROPERTIES SOVERSION ${SOVERSION})
add_dependencies(Accumulators EspressoConfig)
set_target_properties(Accumulators PROPERTIES MACOSX_RPATH TRUE)
target_include_directories(Accumulators PRIVATE ${CMAKE_SOURCE_DIR}/src/core ${CMAKE_BINARY_DIR}/src/core)
