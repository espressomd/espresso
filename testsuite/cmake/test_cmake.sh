#!/bin/bash

# comment out un-necessary tests
TEST_DOXYGEN=1
TEST_SPHINX=1
TEST_PYTHON=1


function set_up() {

if [ ! -f myconfig.hpp ]
then
cp ../../src/core/myconfig-default.hpp myconfig.hpp
fi

cmake ..
make

# run installer
make install DESTDIR="install"

if [ ! -z "${TEST_DOXYGEN}" ]
then
# run Doxygen
make doc
fi

if [ ! -z "${TEST_SPHINX}" ]
then
# run Sphinx
make sphinx
fi

}


# load bash unit testing library
source testsuite/test_cmake/BashUnitTests.sh


# test binaries and libraries are correclty installed
function test_install() {
  local OUTDIR="install"
  for filepath in "${OUTDIR}/usr/local/bin/ipypresso" \
                  "${OUTDIR}/usr/local/bin/pypresso" \
                  "${OUTDIR}/usr/local/lib/libEspressoCore.so" \
                  "${OUTDIR}/usr/local/lib/libEspressoCore.so" \
                  "${OUTDIR}/usr/local/lib/python2.7/site-packages/espressomd/_init.so" \
                  "${OUTDIR}/usr/local/lib/python2.7/site-packages/espressomd/__init__.py" \
                  "${OUTDIR}/usr/local/lib/python3.?/site-packages/espressomd/_init.so" \
                  "${OUTDIR}/usr/local/lib/python3.?/site-packages/espressomd/__init__.py" \
                  ;
  do
    assert_file_exists "${filepath}"
  done
}


if [ ! -z "${TEST_PYTHON}" ]
then
# test if Python module can be imported
function test_Python() {
  local OUTDIR="install"
  assert_return_code python2 -c "import sys;sys.path.insert(0, '${OUTDIR}/usr/local/lib/python2.7/site-packages');import espressomd"
  assert_return_code python3 -c "import sys;sys.path.insert(0, '${OUTDIR}/usr/local/lib/python3.5/site-packages');import espressomd"
}
fi


if [ ! -z "${TEST_DOXYGEN}" ]
then
# test Doxygen documentation
function test_Doxygen() {
  assert_non_zero "$(ls doc/doxygen/html/*_8hpp_source.html    | wc -l)" "No output for .hpp sources"
  assert_non_zero "$(ls doc/doxygen/html/*_8hpp.html           | wc -l)" "No output for C++ headers"
  assert_non_zero "$(ls doc/doxygen/html/*_8cpp_source.html    | wc -l)" "No output for .cpp sources"
  assert_non_zero "$(ls doc/doxygen/html/*_8cpp.html           | wc -l)" "No output for C++ code"
  assert_non_zero "$(ls doc/doxygen/html/*_8cu_source.html     | wc -l)" "No output for .cu sources"
  assert_non_zero "$(ls doc/doxygen/html/*_8cu.html            | wc -l)" "No output for C++ CUDA"
  assert_non_zero "$(ls doc/doxygen/html/struct*.html          | wc -l)" "No output for C++ structs"
  assert_non_zero "$(ls doc/doxygen/html/class*.html           | wc -l)" "No output for C++ classes"
  assert_non_zero "$(ls doc/doxygen/html/namespace*.html       | wc -l)" "No output for C++ namespaces"
  assert_non_zero "$(ls doc/doxygen/html/functions_func_*.html | wc -l)" "No output for C++ functions/methods"
  assert_non_zero "$(ls doc/doxygen/html/functions_vars_*.html | wc -l)" "No output for C++ variables/members"
}
fi


if [ ! -z "${TEST_SPHINX}" ]
then
# test Sphinx documentation
function test_Sphinx() {
  local OUTDIR="install"
  for filepath in "doc/sphinx/modules.rst" \
                  "doc/sphinx/html/modules.html" \
                  "doc/sphinx/html/index.html" \
                  ;
  do
    assert_file_exists "${filepath}"
  done
}
fi



# run tests
run_test_suite

