#include "config/config.hpp"

#ifdef ESPRESSO_METATENSOR
#undef CUDA

#include <torch/cuda.h>
#include <torch/script.h>
#include <torch/version.h>

std::string load_metadata(const std::string &path);

#endif
