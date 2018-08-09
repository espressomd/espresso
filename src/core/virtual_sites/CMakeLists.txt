set(VirtualSites_SRC
    lb_inertialess_tracers.cpp
    lb_inertialess_tracers_cuda_interface.cpp
    virtual_sites_com.cpp
    VirtualSitesInertialessTracers.cpp
    VirtualSitesRelative.cpp
)
add_library(VirtualSites SHARED ${VirtualSites_SRC})
install(TARGETS VirtualSites LIBRARY DESTINATION ${LIBDIR} ARCHIVE DESTINATION ${LIBDIR})
set_target_properties(VirtualSites PROPERTIES SOVERSION ${SOVERSION})
add_dependencies(VirtualSites EspressoConfig)

if(CUDA)
  set(VirtualSitesCuda_SRC
    lb_inertialess_tracers_cuda.cu
  )
  cuda_add_library(VirtualSitesCuda SHARED ${VirtualSitesCuda_SRC})
  install(TARGETS VirtualSitesCuda DESTINATION ${LIBDIR})
  set_target_properties(VirtualSitesCuda PROPERTIES SOVERSION ${SOVERSION})
  add_dependencies(VirtualSitesCuda EspressoConfig)

  CUDA_ADD_CUFFT_TO_TARGET(VirtualSitesCuda)
  target_link_libraries(VirtualSites VirtualSitesCuda)
endif(CUDA)

set_target_properties(VirtualSites PROPERTIES LINKER_LANGUAGE CXX)
target_include_directories(VirtualSites PRIVATE ${CMAKE_SOURCE_DIR}/src/core ${CMAKE_BINARY_DIR}/src/core)
