include(FindPackageHandleStandardArgs)

find_program(SPHINX_EXECUTABLE NAMES sphinx-build
    HINTS
    $ENV{SPHINX_DIR}
    PATHS /opt/local/ /usr/local/ $ENV{HOME}/Library/Python/2.7/
    PATH_SUFFIXES bin
    DOC "Sphinx documentation generator."
)
 
find_program(SPHINX_API_DOC_EXE NAMES sphinx-apidoc
    HINTS
    PATH_SUFFIXES bin
    PATHS /opt/local/ /usr/local/ $ENV{HOME}/Library/Python/2.7/
    DOC "Sphinx api-doc executable."
)

find_package_handle_standard_args(Sphinx DEFAULT_MSG
    SPHINX_EXECUTABLE
)
 
mark_as_advanced(SPHINX_EXECUTABLE)
