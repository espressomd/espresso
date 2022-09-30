add_library(espresso_profiler INTERFACE)
add_library(espresso::profiler ALIAS espresso_profiler)
target_include_directories(espresso_profiler INTERFACE "include")

if(WITH_PROFILER)
  find_package(caliper REQUIRED)

  target_link_libraries(espresso_profiler INTERFACE caliper-mpi)
  target_compile_definitions(espresso_profiler INTERFACE HAVE_CALIPER)
endif()
