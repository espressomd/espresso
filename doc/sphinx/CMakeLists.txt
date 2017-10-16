find_package(Sphinx 1.5.2)
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

    # Copy the figure directory to the build dir
    file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/figures DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

    # Files to be copied to the build directory
    set(FILE_LIST 
        "${CMAKE_SOURCE_DIR}/doc/logo/logo_48x48.png"
        "${CMAKE_CURRENT_SOURCE_DIR}/index.rst"
        "${CMAKE_CURRENT_SOURCE_DIR}/dg.rst"
        "${CMAKE_CURRENT_SOURCE_DIR}/ug.rst"
        "${CMAKE_CURRENT_SOURCE_DIR}/introduction.rst"
        "${CMAKE_CURRENT_SOURCE_DIR}/installation.rst"
        "${CMAKE_CURRENT_SOURCE_DIR}/part.rst"
        "${CMAKE_CURRENT_SOURCE_DIR}/inter.rst"
        "${CMAKE_CURRENT_SOURCE_DIR}/setup.rst"
        "${CMAKE_CURRENT_SOURCE_DIR}/run.rst"
        "${CMAKE_CURRENT_SOURCE_DIR}/analysis.rst"
        "${CMAKE_CURRENT_SOURCE_DIR}/analysis-core.rst"
        "${CMAKE_CURRENT_SOURCE_DIR}/io.rst"
        "${CMAKE_CURRENT_SOURCE_DIR}/lb.rst"
        "${CMAKE_CURRENT_SOURCE_DIR}/electrokinetics.rst"
        "${CMAKE_CURRENT_SOURCE_DIR}/oif.rst"
        "${CMAKE_CURRENT_SOURCE_DIR}/ibm.rst"
        "${CMAKE_CURRENT_SOURCE_DIR}/internal.rst"
        "${CMAKE_CURRENT_SOURCE_DIR}/contributing.rst"
        "${CMAKE_CURRENT_SOURCE_DIR}/features.rst"
        "${CMAKE_CURRENT_SOURCE_DIR}/examples.rst"
        "${CMAKE_CURRENT_SOURCE_DIR}/memd.rst"
        "${CMAKE_CURRENT_SOURCE_DIR}/mmm.rst"
        "${CMAKE_CURRENT_SOURCE_DIR}/cuda.rst"
        "${CMAKE_CURRENT_SOURCE_DIR}/shapes.rst"
        "${CMAKE_CURRENT_SOURCE_DIR}/zrefs.bib"
        )

    foreach(file ${FILE_LIST})
            get_filename_component(basename ${file} NAME) 
            configure_file(
                ${file}
                ${SPHINX_BASE}/${basename}
                COPYONLY
                )
    endforeach()
    
    # Place files to be excluded from SPHINX documentation
    set(EXCLUDE
        "${SPHINX_PYTHON_DIR}/gen_code_info.py"
        "${SPHINX_PYTHON_DIR}/gen_pxiconfig*"
    )

    add_custom_target(sphinx
        ${SPHINX_API_DOC_EXE}
            -f
            -P
            -o ${SPHINX_BASE}
            ${SPHINX_PYTHON_DIR}
            ${EXCLUDE}
        COMMAND
        ${SPHINX_EXECUTABLE}
            -q
            -b html
            -c "${SPHINX_BASE}"
            -d "${SPHINX_CACHE_DIR}"
            "${CMAKE_CURRENT_BINARY_DIR}"
            "${SPHINX_HTML_DIR}"
        DEPENDS espressomd pypresso 
        COMMENT "Building HTML documentation with Sphinx")
endif(SPHINX_FOUND) 
