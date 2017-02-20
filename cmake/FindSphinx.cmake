include(FindPackageHandleStandardArgs)

find_program(SPHINX_EXECUTABLE NAMES sphinx-build sphinx-build-${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}
    HINTS
    $ENV{SPHINX_DIR}
    PATHS /opt/local/ /usr/local/ $ENV{HOME}/Library/Python/2.7/
    PATH_SUFFIXES bin
    DOC "Sphinx documentation generator."
)
 
find_program(SPHINX_API_DOC_EXE NAMES sphinx-apidoc sphinx-apidoc-${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}
    HINTS
    PATH_SUFFIXES bin
    PATHS /opt/local/ /usr/local/ $ENV{HOME}/Library/Python/2.7/
    DOC "Sphinx api-doc executable."
)

if( Sphinx_FIND_VERSION )
  execute_process(COMMAND "${SPHINX_EXECUTABLE}" --version
                  OUTPUT_VARIABLE QUERY_VERSION
                  RESULT_VARIABLE QUERY_VERSION_RESULT
                  ERROR_QUIET)
  if( NOT QUERY_VERSION_RESULT )
    string(REGEX MATCH "[0-9.]+" AVAILABLE_VERSION "${QUERY_VERSION}")
    if( "${Sphinx_FIND_VERSION}" VERSION_GREATER "${AVAILABLE_VERSION}" )
      message(FATAL_ERROR "Sphinx version ${Sphinx_FIND_VERSION} requested but only ${AVAILABLE_VERSION} found!")
    endif( )
  endif( )
endif( )

find_package_handle_standard_args(Sphinx DEFAULT_MSG
    SPHINX_EXECUTABLE
)
 
mark_as_advanced(SPHINX_EXECUTABLE)
