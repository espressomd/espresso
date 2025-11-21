export PATH=/usr/local/cuda-12.8/bin:$PATH
export LD_LIBRARY_PATH=/usr/local/cuda-12.8/lib64:$LD_LIBRARY_PATH
export CUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda-12.8

ESPRESSO_DIR=/tikhome/weeber/es
BUILD_DIR=/tikhome/weeber/es/build_mt
VENV_DIR=/work/weeber/es_mt-env

rm -rf ${BUILD_DIR}
mkdir -p ${BUILD_DIR}
cd ${BUILD_DIR}

source ${VENV_DIR}/bin/activate

export TORCH_PREFIX=$(python -c "import torch; print(torch.utils.cmake_prefix_path)")
export MTS_PREFIX=$(python -c "import metatensor; print(metatensor.utils.cmake_prefix_path)")
export MTS_TORCH_PREFIX=$(python -c "import metatensor.torch; print(metatensor.torch.utils.cmake_prefix_path)")
export MTA_TORCH_PREFIX=$(python -c "import metatomic.torch; print(metatomic.torch.utils.cmake_prefix_path)")
export CMAKE_PREFIX_PATH="$TORCH_PREFIX;$MTS_PREFIX;$MTS_TORCH_PREFIX;$MTA_TORCH_PREFIX"

cd ${BUILD_DIR}
cmake ../ \
  -D CMAKE_BUILD_TYPE=Debug \
  -D ESPRESSO_BUILD_WITH_CUDA=OFF \
  -D ESPRESSO_BUILD_WITH_CCACHE=OFF \
  -D ESPRESSO_BUILD_WITH_WALBERLA=OFF \
  -D ESPRESSO_BUILD_WITH_WALBERLA_AVX=OFF \
  -D ESPRESSO_BUILD_WITH_GSL=OFF \
  -D ESPRESSO_BUILD_WITH_METATENSOR=ON \
  -D CMAKE_PREFIX_PATH="$CMAKE_PREFIX_PATH"

make -j$(nproc)
