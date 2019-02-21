
add_library(Profiler SHARED src/profiler.cpp)
target_include_directories(Profiler PUBLIC "include")

if(WITH_PROFILER)
find_package(caliper REQUIRED)

target_link_libraries(Profiler PUBLIC caliper-mpi)
target_compile_definitions(Profiler PUBLIC HAVE_CALIPER)
endif()

