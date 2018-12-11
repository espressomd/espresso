if(CUDA)
  set(VirtualSitesCuda_SRC
    lb_inertialess_tracers_cuda.cu
  )
  cuda_add_library(VirtualSitesCuda SHARED ${VirtualSitesCuda_SRC})
  install(TARGETS VirtualSitesCuda DESTINATION ${PYTHON_INSTDIR}/espressomd)
  set_target_properties(VirtualSitesCuda PROPERTIES SOVERSION ${SOVERSION})
  target_link_libraries(VirtualSitesCuda EspressoConfig)

  CUDA_ADD_CUFFT_TO_TARGET(VirtualSitesCuda)
endif(CUDA)
