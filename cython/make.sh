
cd test_bin
make -j4
cd ..
cython espresso.pyx  && \
gcc -m32 -fPIC -shared -c ./errexit.c && \
gcc -std=c99 -m32 -fPIC -shared -c espresso.c  -I/Library/Frameworks/Python.framework/Versions/Current/Headers -I../src -Itest_bin/src\
&& gcc -std=c99 -m32 -shared -o espresso.so errexit.o espresso.o  test_bin/Espresso-scriptsdir.o  \
    -lpython2.7 -L/Library/Frameworks/Python.framework/Versions/Current/lib  \
    -ltcl \
    -lEspresso -L./test_bin/src \
