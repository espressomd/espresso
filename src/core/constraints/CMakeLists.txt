set(Constraints_SRC
    HomogeneousMagneticField.cpp
    ShapeBasedConstraint.cpp
)
add_library(Constraints SHARED ${Constraints_SRC})
install(TARGETS Constraints LIBRARY DESTINATION ${PYTHON_INSTDIR}/espressomd)
set_target_properties(Constraints PROPERTIES SOVERSION ${SOVERSION})
target_link_libraries(Constraints PRIVATE EspressoConfig)
target_include_directories(Constraints PRIVATE ${CMAKE_SOURCE_DIR}/src/core ${CMAKE_BINARY_DIR}/src/core)
