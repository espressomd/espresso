set(ObjectInFluid_SRC
    affinity.cpp
    membrane_collision.cpp
    oif_global_forces.cpp
    oif_local_forces.cpp
    out_direction.cpp
)
add_library(ObjectInFluid SHARED ${ObjectInFluid_SRC})
install(TARGETS ObjectInFluid LIBRARY DESTINATION ${LIBDIR} ARCHIVE DESTINATION ${LIBDIR})
set_target_properties(ObjectInFluid PROPERTIES SOVERSION ${SOVERSION})
add_dependencies(ObjectInFluid EspressoConfig)

