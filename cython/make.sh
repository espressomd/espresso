
if [ -e test_bin ]; then
  echo "Yeah you have an Espresso binary in the right place"
else
  echo "please compile espresso into a subdirectory test_bin of cython"
  echo "please add the Compiler Flag -fPIC"
fi

cd test_bin
make -j4
cd ..
cython espresso.pyx  && \
gcc -std=c99 -fPIC -shared -c ./errexit.c && \
gcc -std=c99 -fPIC -shared -c espresso.c  -I/usr/local/epd-7.0-2_64/include/python2.7 -I../src -Itest_bin/src\
&& gcc -std=c99 -shared -o espresso.so errexit.o espresso.o  test_bin/Espresso-scriptsdir.o  \
    -lpython2.7 -L/usr/local/epd-7.0-2_64/lib \
    -ltcl8.5 \
    -lEspresso -L./test_bin/src \
    -lmpi \
    -lfftw3
