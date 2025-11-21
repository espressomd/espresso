#include "config/config.hpp"

#ifdef ESPRESSO_METATENSOR
#include <memory>
#include <metatensor/torch.hpp>
#include <metatomic/torch.hpp>

#include <torch/torch.h>

struct MetatensorModel {
  MetatensorModel(const std::string &path);
  std::string print_metadata();

  std::unique_ptr<torch::jit::Module> model;
  metatomic_torch::ModelCapabilities capabilities;
  metatomic_torch::ModelEvaluationOptions evaluation_options;
};
#endif
