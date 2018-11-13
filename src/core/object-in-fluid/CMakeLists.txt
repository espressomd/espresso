set(ObjectInFluid_SRC
    affinity.cpp
    membrane_collision.cpp
    oif_global_forces.cpp
    oif_local_forces.cpp
    out_direction.cpp
)
add_library(ObjectInFluid SHARED ${ObjectInFluid_SRC})
install(TARGETS ObjectInFluid LIBRARY DESTINATION ${PYTHON_INSTDIR}/espressomd ARCHIVE DESTINATION ${PYTHON_INSTDIR}/espressomd)
set_target_properties(ObjectInFluid PROPERTIES SOVERSION ${SOVERSION})
target_include_directories(ObjectInFluid PRIVATE ${CMAKE_SOURCE_DIR}/src/core ${CMAKE_BINARY_DIR}/src/core)
add_dependencies(ObjectInFluid EspressoConfig)

