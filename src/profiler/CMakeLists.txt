add_library(Profiler INTERFACE)
target_include_directories(Profiler INTERFACE "include")

if(WITH_PROFILER)
  find_package(caliper REQUIRED)

  target_link_libraries(Profiler INTERFACE caliper-mpi)
  target_compile_definitions(Profiler INTERFACE HAVE_CALIPER)
endif()
