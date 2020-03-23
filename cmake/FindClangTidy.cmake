# Find Clang-tidy

# get Clang version
string(REGEX
        REPLACE "^([1-9]+)\\.[0-9]+.*$"
        "\\1"
        CLANG_MAJOR_VERSION
        "${CMAKE_CXX_COMPILER_VERSION}")
string(REGEX
       REPLACE "^[1-9]+\\.([0-9]+).*$"
               "\\1"
               CLANG_MINOR_VERSION
               "${CMAKE_CXX_COMPILER_VERSION}")
# find Clang-tidy
find_program(CLANG_TIDY_EXE
             NAMES "clang-tidy-${CLANG_MAJOR_VERSION}.${CLANG_MINOR_VERSION}" "clang-tidy-${CLANG_MAJOR_VERSION}" "clang-tidy"
             DOC "Path to clang-tidy executable")
if(CLANG_TIDY_EXE)
  execute_process(COMMAND ${CLANG_TIDY_EXE} --version
                OUTPUT_VARIABLE CLANG_TIDY_OUTPUT RESULT_VARIABLE CLANG_TIDY_STATUS ERROR_QUIET)
  if(CLANG_TIDY_STATUS EQUAL 0)
    string(REGEX MATCH "LLVM version ([0-9]+\\.[0-9\\.]+)" CLANG_TIDY_VERSION "${CLANG_TIDY_OUTPUT}")
    string(REGEX MATCH "([0-9\\.]+)" CLANG_TIDY_VERSION "${CLANG_TIDY_VERSION}")
  endif()
endif()

include( FindPackageHandleStandardArgs )
FIND_PACKAGE_HANDLE_STANDARD_ARGS( ClangTidy REQUIRED_VARS CLANG_TIDY_EXE
                                   VERSION_VAR CLANG_TIDY_VERSION)
