#include "config/config.hpp"

#ifdef ESPRESSO_METATENSOR
#undef CUDA
#include <torch/cuda.h>
#include <torch/script.h>
#include <torch/version.h>

#if TORCH_VERSION_MAJOR >= 2
#include <torch/mps.h>
#endif

#include <metatensor/torch.hpp>
#include <metatomic/torch.hpp>

#include "model.hpp"
#include "stub.hpp"

std::string load_metadata(const std::string &path) {
	MetatensorModel model(path);
	return model.print_metadata();
}
#endif
