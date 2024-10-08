#include "config/config.hpp"

#ifdef METATENSOR
#undef CUDA

#include <torch/cuda.h>
#include <torch/script.h>
#include <torch/version.h>
#endif
