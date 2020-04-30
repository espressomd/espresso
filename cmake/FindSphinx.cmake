include(FindPackageHandleStandardArgs)

set(SPHINX_EXECUTABLE ${PYTHON_EXECUTABLE} -m sphinx)

execute_process(
  COMMAND ${SPHINX_EXECUTABLE} --version OUTPUT_VARIABLE QUERY_VERSION_OUT
  ERROR_VARIABLE QUERY_VERSION_ERR RESULT_VARIABLE QUERY_VERSION_RESULT)

if(NOT QUERY_VERSION_RESULT)
  # Sphinx switched at some point from returning their version on stdout to
  # printing it at stderr. Since we do not know their version yet, we use stdout
  # if it matches a version regex, or stderr otherwise.
  if(QUERY_VERSION_OUT MATCHES "[0-9]+\.[0-9.]+")
    set(QUERY_VERSION "${QUERY_VERSION_OUT}")
  else()
    set(QUERY_VERSION "${QUERY_VERSION_ERR}")
  endif()

  string(REGEX MATCH "[0-9]+\.[0-9.]+" SPHINX_VERSION "${QUERY_VERSION}")

  if("${SPHINX_VERSION}" VERSION_LESS "1.7")
    set(SPHINX_API_DOC_EXE ${PYTHON_EXECUTABLE} -m sphinx.apidoc)
  else()
    set(SPHINX_API_DOC_EXE ${PYTHON_EXECUTABLE} -m sphinx.ext.apidoc)
  endif()
endif()

set(SPHINX_VERSION_COMPATIBLE TRUE)
# Blacklist broken versions
if("${SPHINX_VERSION}" VERSION_EQUAL "2.1.0" OR "${SPHINX_VERSION}" VERSION_EQUAL "3.0.0")
  message(WARNING "Sphinx version ${SPHINX_VERSION} is not compatible.")
  set(SPHINX_VERSION_COMPATIBLE FALSE)
endif()

find_package_handle_standard_args(
  Sphinx REQUIRED_VARS SPHINX_EXECUTABLE SPHINX_API_DOC_EXE
  SPHINX_VERSION_COMPATIBLE VERSION_VAR SPHINX_VERSION)

mark_as_advanced(SPHINX_EXECUTABLE)
mark_as_advanced(SPHINX_API_DOC_EXE)
