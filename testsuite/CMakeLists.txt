include(unit_test)

if(WITH_COVERAGE_PYTHON)
  set(PYPRESSO_OPTIONS --coverage)
endif(WITH_COVERAGE_PYTHON)

add_subdirectory(python)
add_subdirectory(cmake)
add_subdirectory(scripts)
