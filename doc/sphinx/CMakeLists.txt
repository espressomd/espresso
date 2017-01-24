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
     
    # Sphinx configuration file
    configure_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/conf.py.in"
        "${SPHINX_BASE}/conf.py"
        @ONLY)
    
    # Files to be copied to the build directory
    set(FILE_LIST 
        "${CMAKE_SOURCE_DIR}/doc/logo/logo_48x48.png"
        "${CMAKE_CURRENT_SOURCE_DIR}/index.rst"
        "${CMAKE_CURRENT_SOURCE_DIR}/defs.rst"
        "${CMAKE_CURRENT_SOURCE_DIR}/introduction.rst"
        "${CMAKE_CURRENT_SOURCE_DIR}/analysis-core.tex.rst"
        "${CMAKE_CURRENT_SOURCE_DIR}/analysis.tex.rst"
        "${CMAKE_CURRENT_SOURCE_DIR}/aux.tex.rst"
        "${CMAKE_CURRENT_SOURCE_DIR}/contributing.tex.rst"
        "${CMAKE_CURRENT_SOURCE_DIR}/electrokinetics.tex.rst"
        "${CMAKE_CURRENT_SOURCE_DIR}/examples.tex.rst"
        "${CMAKE_CURRENT_SOURCE_DIR}/features.tex.rst"
        "${CMAKE_CURRENT_SOURCE_DIR}/ibm.tex.rst"
        "${CMAKE_CURRENT_SOURCE_DIR}/installation.tex.rst"
        "${CMAKE_CURRENT_SOURCE_DIR}/inter.tex.rst"
        "${CMAKE_CURRENT_SOURCE_DIR}/internal.tex.rst"
        "${CMAKE_CURRENT_SOURCE_DIR}/lb.tex.rst"
        "${CMAKE_CURRENT_SOURCE_DIR}/mbtools.tex.rst"
        "${CMAKE_CURRENT_SOURCE_DIR}/memd.tex.rst"
        "${CMAKE_CURRENT_SOURCE_DIR}/mmm.tex.rst"
        "${CMAKE_CURRENT_SOURCE_DIR}/run.tex.rst"
        )

    foreach(file ${FILE_LIST})
            get_filename_component(basename ${file} NAME) 
            configure_file(
                ${file}
                ${SPHINX_BASE}/${basename}
                COPYONLY
                )
    endforeach()

    add_custom_target(sphinx ALL
        ${SPHINX_API_DOC_EXE}
            -f
            -o ${SPHINX_BASE}
            ${SPHINX_PYTHON_DIR}
        COMMAND
        ${SPHINX_EXECUTABLE}
            -b html
            -c "${SPHINX_BASE}"
            -d "${SPHINX_CACHE_DIR}"
            "${CMAKE_CURRENT_BINARY_DIR}"
            "${SPHINX_HTML_DIR}"
        DEPENDS espressomd pypresso 
        COMMENT "Building HTML documentation with Sphinx")
endif(SPHINX_FOUND) 
