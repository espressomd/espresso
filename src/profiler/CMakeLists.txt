add_library(Espresso_profiler INTERFACE)
add_library(Espresso::profiler ALIAS Espresso_profiler)
target_include_directories(Espresso_profiler INTERFACE "include")

if(WITH_PROFILER)
  find_package(caliper REQUIRED)

  target_link_libraries(Espresso_profiler INTERFACE caliper-mpi)
  target_compile_definitions(Espresso_profiler INTERFACE HAVE_CALIPER)
endif()
