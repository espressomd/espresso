
all: libespresso_main.a ./espresso/espresso.so

./espresso/espresso.so:
	cd espresso ; python setup.py build_ext --inplace ; cd ..

libespresso_main.a: espresso_bin/src/libEspresso.a errexit.o 
	ar rcs libespresso_main.a ./espresso_bin/src/*.o ./espresso_bin/Espresso-scriptsdir.o ./errexit.o

errexit.o:
	gcc -c errexit.c

espresso_bin:
	mkdir espresso_bin

espresso_bin/src/libEspresso.a: ../configure
	@cd espresso_bin ; ../../configure && make -j  ; cd ..

../configure:
	@cd .. &&  ./bootstrap.sh && cd cython
