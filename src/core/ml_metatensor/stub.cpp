#include  "config/config.hpp"

#ifdef METATENSOR
#undef CUDA
#include <torch/version.h>
#include <torch/script.h>
#include <torch/cuda.h>

#if TORCH_VERSION_MAJOR >= 2
    #include <torch/mps.h>
#endif

#include <metatensor/torch.hpp>
#include <metatensor/torch/atomistic.hpp>
#endif
