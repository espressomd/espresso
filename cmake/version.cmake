if(GIT_EXECUTABLE)
  # Get the name of the working branch
  execute_process(
    COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} OUTPUT_VARIABLE GIT_BRANCH
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  # Get the latest abbreviated commit hash of the working branch
  execute_process(
    COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_COMMIT_HASH OUTPUT_STRIP_TRAILING_WHITESPACE)

  # Get branch status
  execute_process(
    COMMAND ${GIT_EXECUTABLE} diff-index --quiet HEAD --
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    RESULT_VARIABLE GIT_DIFF_INDEX_RESULT
    OUTPUT_VARIABLE GIT_DIFF_INDEX_OUTPUT OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(GIT_DIFF_INDEX_RESULT EQUAL 0)
    set(GIT_STATE "CLEAN")
  else()
    set(GIT_STATE "DIRTY")
  endif()

endif(GIT_EXECUTABLE)

configure_file(${PROJECT_SOURCE_DIR}/src/config/version.hpp.in version.hpp.tmp)
execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different version.hpp.tmp
                        version.hpp)
