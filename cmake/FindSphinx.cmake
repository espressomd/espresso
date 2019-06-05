include(FindPackageHandleStandardArgs)

find_program(SPHINX_EXECUTABLE
             NAMES sphinx-build
                   sphinx-build-${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}
             HINTS $ENV{SPHINX_DIR}
             PATHS /opt/local/ /usr/local/ $ENV{HOME}/Library/Python/2.7/
             PATH_SUFFIXES bin
             DOC "Sphinx documentation generator.")

find_program(SPHINX_API_DOC_EXE
             NAMES sphinx-apidoc
                   sphinx-apidoc-${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}
             HINTS
             PATH_SUFFIXES bin
             PATHS /opt/local/ /usr/local/ $ENV{HOME}/Library/Python/2.7/
             DOC "Sphinx api-doc executable.")

if(Sphinx_FIND_VERSION)
  execute_process(COMMAND "${SPHINX_EXECUTABLE}" --version
                  OUTPUT_VARIABLE QUERY_VERSION_OUT
                  ERROR_VARIABLE QUERY_VERSION_ERR
                  RESULT_VARIABLE QUERY_VERSION_RESULT)

  # Sphinx switched at some point from returning ther version on stdout to
  # printing it at stderr. Since we do not know ther version yet, we use stdout
  # if it matches a version regex, or stderr otherwise.
  if(QUERY_VERSION_OUT MATCHES "[0-9.]+")
    set(QUERY_VERSION "${QUERY_VERSION_OUT}")
  else()
    set(QUERY_VERSION "${QUERY_VERSION_ERR}")
  endif()

  if(NOT QUERY_VERSION_RESULT)
    string(REGEX MATCH
                 "[0-9.]+"
                 SPHINX_VERSION
                 "${QUERY_VERSION}")
  endif()

  set(SPHINX_VERSION_COMPATIBLE TRUE)
  # Blacklist broken version
  if("${SPHINX_VERSION}" VERSION_EQUAL "2.1.0")
    message(WARNING "Sphinx version 2.1.0 is not compatible.")
    set(SPHINX_VERSION_COMPATIBLE FALSE)
  endif()
endif()

find_package_handle_standard_args(Sphinx
                                  REQUIRED_VARS
                                  SPHINX_EXECUTABLE
                                  SPHINX_API_DOC_EXE
                                  SPHINX_VERSION_COMPATIBLE
                                  VERSION_VAR
                                  SPHINX_VERSION)

mark_as_advanced(SPHINX_EXECUTABLE)
