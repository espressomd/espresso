# find the Python C++ headers
execute_process(
  COMMAND ${PYTHON_EXECUTABLE} -c
          "import distutils.sysconfig as cg; print(cg.get_python_inc())"
  OUTPUT_VARIABLE PYTHON_INCLUDE_DIRS
  OUTPUT_STRIP_TRAILING_WHITESPACE)
# find Python installation directory
if(NOT PYTHON_INSTDIR)
  execute_process(
    COMMAND
      ${PYTHON_EXECUTABLE} -c
      "import distutils.sysconfig as cg; print(cg.get_python_lib(prefix='${CMAKE_INSTALL_PREFIX}', plat_specific=True, standard_lib=False).replace('${CMAKE_INSTALL_PREFIX}/', '', 1))"
    OUTPUT_VARIABLE PYTHON_INSTDIR
    OUTPUT_STRIP_TRAILING_WHITESPACE)
endif(NOT PYTHON_INSTDIR)
include( FindPackageHandleStandardArgs )
FIND_PACKAGE_HANDLE_STANDARD_ARGS( PythonHeaders REQUIRED_VARS PYTHON_INCLUDE_DIRS PYTHON_INSTDIR )
