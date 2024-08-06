using ModelPtr = std::unique_ptr<torch::jit::Module>;
using NeighborListRequest = 
   std::pair<double, metatensor_torch::NeighborListOptions>;


ModelPtr load_model(const std::string& path, const std::string& extensions_directory, torch::device device) {

 return std::make_unique<torch::jit::Module>(
            metatensor_torch::load_atomistic_model(path, extensions)
        );
}


metatensor_torch::ModelCapabilitiesHolder 
get_model_capabilites(const ModelPtr& model) {
  auto capabilities_ivalue = model->run_method("capabilities");
  return capabilities_ivalue.toCustomClass<metatensor_torch::ModelCapabilitiesHolder>();
};

bool modle_provides_energy(const torch_metatensor::ModelCapabilitiesHolder& capabilities) {
  return (capabilities->outputs().contains("energy"));
}


metatensor_toch::ModelMetadataHolder get_model_metadata(ModelPtr& model) {
    auto metadata_ivalue = model->run_method("metadata");
  return metadata_ivalue.toCustomClass<metatensor_torch::ModelMetadataHolder>();
}



double required_range(const metatensor_torch::ModelCapabilities& capabilities, const metatensor_torch::ModelEvaluatoinOptoins& evaluation_optoins) {
    return range = mts_data->capabilities->engine_interaction_range(evaluation_options->length_unit());
}

std::vector<NeighborListRequest> get_requested_neighbor_lists(ModelPtr& model) {

    std::vector<NeighborListRequest> res;
    auto requested_nl = mts_data->model->run_method("requested_neighbor_lists");
    for (const auto& ivalue: requested_nl.toList()) {
        auto options = ivalue.get().toCustomClass<metatensor_torch::NeighborListOptionsHolder>();
        auto cutoff = options->engine_cutoff(mts_data->evaluation_options->length_unit());

        res.push_back({cutoff, options});
    }
   return res;
}
