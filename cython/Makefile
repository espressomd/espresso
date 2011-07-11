

PYXFILES=$(wildcard espresso/*.pyx)
PYXSOS=$(patsubst %.pyx,%.so, ${PYXFILES})

all: libespresso_main.so cython

test: 
	@echo ${PYXFILES}
	@echo ${PYXSOS}

.PHONY: 
cython:
	cd espresso ; python setup.py build_ext --inplace ; cd ..


libespresso_main.so: espresso_bin/src/libEspresso.a errexit.o 
	gcc -shared -o libespresso_main.so ./espresso_bin/src/*.o ./espresso_bin/Espresso-scriptsdir.o ./errexit.o

errexit.o:
	gcc -c errexit.c

espresso_bin/src/libEspresso.a: ../configure
	@if [ -e espresso_bin ]; then echo "espresso_bin already exists"; else mkdir espresso_bin; fi; cd espresso_bin ; ../../configure CFLAGS="-fPIC" && make -j  ; cd ..

../configure:
	@cd .. &&  ./bootstrap.sh && cd cython
