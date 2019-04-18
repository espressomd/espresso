add_library(utils INTERFACE)
target_include_directories(utils INTERFACE include)

if(WITH_TESTS)
  add_subdirectory(tests)
endif()

