#include "config/config.hpp"

#ifdef METATENSOR
#undef CUDA

#include <torch/version.h>
#include <torch/script.h>
#include <torch/cuda.h>
#endif
