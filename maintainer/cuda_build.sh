cd ${CI_PROJECT_DIR};
cp maintainer/configs/maxset.hpp myconfig.hpp;
mkdir build && cd build;
cmake .. -DCUDA_NVCC_FLAGS="-D_FORCE_INLINES";
make && make check_python; cat testsuite/python/Testing/Temporary/LastTest.log
