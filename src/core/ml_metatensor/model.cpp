#include "config/config.hpp"

#ifdef ESPRESSO_METATENSOR
#include "model.hpp"
#include <string>

MetatensorModel::MetatensorModel(const std::string &path) {
  this->model = std::make_unique<torch::jit::Module>(
      metatomic_torch::load_atomistic_model(
          path)); // TODO: Capture c10::Error and add custom exception?

  auto capabilities_ivalue = this->model->run_method("capabilities");
  this->capabilities =
      capabilities_ivalue
          .toCustomClass<metatomic_torch::ModelCapabilitiesHolder>();

  if (!this->capabilities->outputs().contains("energy")) {
    throw std::runtime_error("Metatensor model at " + path +
                             " does not have energy as output!");
  }
}

std::string MetatensorModel::print_metadata() {
  if (!this->model) {
    throw std::runtime_error("No model loaded!");
  }
  auto metadata_ivalue = this->model->run_method("metadata");
  auto metadata =
      metadata_ivalue.toCustomClass<metatomic_torch::ModelMetadataHolder>();
  return metadata->print();
}
#endif
