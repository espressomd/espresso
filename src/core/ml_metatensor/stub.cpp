#include "config/config.hpp"

#ifdef METATENSOR
#undef CUDA
#include <torch/cuda.h>
#include <torch/script.h>
#include <torch/version.h>

#if TORCH_VERSION_MAJOR >= 2
#include <torch/mps.h>
#endif

#include <metatensor/torch.hpp>
#include <metatensor/torch/atomistic.hpp>
#endif
