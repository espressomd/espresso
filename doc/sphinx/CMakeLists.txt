find_package(Sphinx)
if(SPHINX_FOUND) 
    # configured documentation tools and
    # intermediate build results
    set(SPHINX_BASE "${CMAKE_CURRENT_BINARY_DIR}")
    set(BINARY_BUILD_DIR "${SPHINX_BASE}/_build")
    set(SPHINX_PYTHON_DIR "${CMAKE_SOURCE_DIR}/src/python/espressomd") 
    # Sphinx cache with pickled ReST documents
    set(SPHINX_CACHE_DIR "${SPHINX_BASE}/_doctrees")
     
    # HTML output directory
    set(SPHINX_HTML_DIR "${SPHINX_BASE}/html")
     
    configure_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/conf.py.in"
        "${SPHINX_BASE}/conf.py"
        @ONLY)

    configure_file(
        "${CMAKE_SOURCE_DIR}/doc/logo/logo_48x48.png"
        "${SPHINX_BASE}"
        COPYONLY)
    
    configure_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/index.rst"
        "${SPHINX_BASE}/index.rst"
        COPYONLY)
    
    configure_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/introduction.rst"
        "${SPHINX_BASE}/introduction.rst"
        COPYONLY)

    add_custom_target(sphinx ALL
        ${SPHINX_API_DOC_EXE}
            -f
            -o ${SPHINX_BASE}
            ${SPHINX_PYTHON_DIR}
        COMMAND
        ${SPHINX_EXECUTABLE}
            -Q -b html
            -c "${SPHINX_BASE}"
            -d "${SPHINX_CACHE_DIR}"
            "${CMAKE_CURRENT_BINARY_DIR}"
            "${SPHINX_HTML_DIR}"
        DEPENDS espressomd pypresso 
        COMMENT "Building HTML documentation with Sphinx")

endif(SPHINX_FOUND) 
