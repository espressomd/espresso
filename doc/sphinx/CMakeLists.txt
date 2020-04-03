find_package(Sphinx 1.6.6)
if(SPHINX_FOUND)
  # configured documentation tools and intermediate build results
  set(SPHINX_BASE "${CMAKE_CURRENT_BINARY_DIR}")
  set(BINARY_BUILD_DIR "${SPHINX_BASE}/_build")
  set(SPHINX_PYTHON_DIR "${CMAKE_SOURCE_DIR}/src/python/espressomd")
  # Sphinx cache with pickled ReST documents
  set(SPHINX_CACHE_DIR "${SPHINX_BASE}/_doctrees")

  # HTML output directory
  set(SPHINX_HTML_DIR "${SPHINX_BASE}/html")

  # Sphinx configuration file
  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/conf.py.in"
                 "${SPHINX_BASE}/conf.py" @ONLY)

  # Copy the figure directory to the build dir
  file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/figures
       DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

  # Copy the _static directory to the build dir
  file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/_static
       DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

  configure_file(samples.py ${CMAKE_CURRENT_BINARY_DIR}/samples.py)

  # Files to be copied to the build directory
  set(FILE_LIST
      "${CMAKE_SOURCE_DIR}/doc/logo/logo_48x48.png"
      "${CMAKE_CURRENT_SOURCE_DIR}/advanced_methods.rst"
      "${CMAKE_CURRENT_SOURCE_DIR}/analysis.rst"
      "${CMAKE_CURRENT_SOURCE_DIR}/appendix.rst"
      "${CMAKE_CURRENT_SOURCE_DIR}/constraints.rst"
      "${CMAKE_CURRENT_SOURCE_DIR}/contributing.rst"
      "${CMAKE_CURRENT_SOURCE_DIR}/electrostatics.rst"
      "${CMAKE_CURRENT_SOURCE_DIR}/index.rst"
      "${CMAKE_CURRENT_SOURCE_DIR}/installation.rst"
      "${CMAKE_CURRENT_SOURCE_DIR}/inter_bonded.rst"
      "${CMAKE_CURRENT_SOURCE_DIR}/inter_non-bonded.rst"
      "${CMAKE_CURRENT_SOURCE_DIR}/introduction.rst"
      "${CMAKE_CURRENT_SOURCE_DIR}/io.rst"
      "${CMAKE_CURRENT_SOURCE_DIR}/lb.rst"
      "${CMAKE_CURRENT_SOURCE_DIR}/magnetostatics.rst"
      "${CMAKE_CURRENT_SOURCE_DIR}/particles.rst"
      "${CMAKE_CURRENT_SOURCE_DIR}/running.rst"
      "${CMAKE_CURRENT_SOURCE_DIR}/system_manipulation.rst"
      "${CMAKE_CURRENT_SOURCE_DIR}/system_setup.rst"
      "${CMAKE_CURRENT_SOURCE_DIR}/ug.rst"
      "${CMAKE_CURRENT_SOURCE_DIR}/under_the_hood.rst"
      "${CMAKE_CURRENT_SOURCE_DIR}/visualization.rst"
      "${CMAKE_CURRENT_SOURCE_DIR}/zrefs.bib"
      "${CMAKE_CURRENT_SOURCE_DIR}/zreferences.rst")

  foreach(file ${FILE_LIST})
    get_filename_component(basename ${file} NAME)
    configure_file(${file} ${SPHINX_BASE}/${basename} COPYONLY)
  endforeach()

  # Place files to be excluded from SPHINX documentation
  set(EXCLUDE "${SPHINX_PYTHON_DIR}/gen_code_info.py"
              "${SPHINX_PYTHON_DIR}/gen_pxiconfig*")

  add_custom_target(
    sphinx
    ${PYTHON_EXECUTABLE} "${CMAKE_CURRENT_BINARY_DIR}/samples.py"
    COMMAND ${SPHINX_API_DOC_EXE} -f -o ${SPHINX_BASE} ${SPHINX_PYTHON_DIR}
            ${EXCLUDE}
    COMMAND
      ${SPHINX_EXECUTABLE} -q -W -b html -c "${SPHINX_BASE}" -d
      "${SPHINX_CACHE_DIR}" "${CMAKE_CURRENT_BINARY_DIR}" "${SPHINX_HTML_DIR}"
    DEPENDS espressomd pypresso
    COMMENT "Building HTML documentation with Sphinx.")
  add_custom_target(
    doctest ${SPHINX_EXECUTABLE} -b doctest -d "${SPHINX_CACHE_DIR}"
            "${CMAKE_CURRENT_BINARY_DIR}" "${SPHINX_HTML_DIR}" DEPENDS sphinx
    COMMENT "Performing doctests within the Sphinx documentation.")
endif(SPHINX_FOUND)
