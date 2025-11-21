struct PairInfo {
  int part_id_1,
  int part_id_2,
  Utils::Vector3d distance;
}

using Sample = std::array<int_32_t,5>;
using Distances =
    std::variant<std::vector<std::array<double,3>>, std::vector<std::array<float,3>>>;


template <typename PairIterable>
TorchTensorBlock neighbor_list_from_pairs(const metatensor_torch::System& system, const PairIterable& pairs) {
    auto dtype = system->positions().scalar_type();
    auto device = system->positions().device();
    std::vector<Sample> samples;
   Distances distances; 
   if (dtype == torch::kFloat64) {
     distances = {std::vector<std::array<double,3>>()};
   }
   else if (dtype == torch::kFloat32) {
     distances = {std::vector<std::array<float,3>>()};
   }
   else {
    throw std::runtime_error("Unsupported floating poitn data type");
   }

   for (auto const& pair: pairs) {
       auto sample = Sample{
         pair.particle_id_1, pair.particle_id_2, 0, 0, 0};
       samples.push_back(sample);
       (*distances).push_back(pair.distance);
    }


    int64_t n_pairs = samples.size();
    auto samples_tensor = torch::from_blob(
            reinterpret_cast<int32_t*>(samples.data()),
            {n_pairs, 5},
            torch::TensorOptions().dtype(torch::kInt32).device(torch::kCPU)
        );

        auto samples = torch::make_intrusive<metatensor_torch::LabelsHolder>(
            std::vector<std::string>{"first_atom", "second_atom", "cell_shift_a", "cell_shift_b", "cell_shift_c"},
            samples_values
        );

       distances_vectors = torch::from_blob(
                (*distances).data(),
                {n_pairs, 3, 1},
                torch::TensorOptions().dtype(dtype).device(torch::kCPU)
            );
        return neighbors = torch::make_intrusive<metatensor_torch::TensorBlockHolder>(
            distances_vectors.to(dtype).to(device),
            samples->to(device),
            std::vector<metatensor_torch::TorchLabels>{
                metatensor_torch::LabelsHolder::create({"xyz"}, {{0}, {1}, {2}})->to(device),
            },
            metatensor_torch::LabelsHolder::create({"distance"}, {{0}})->to(device)
        );

}

void add_neighbor_list_to_system(MetatensorTorch::system& system, 
    const TorchTensorBlock& neighbors,
    const NeighborListOptions& options) {
        metatensor_torch::register_autograd_neighbors(system, neighbors, options_.check_consistency);
        system->add_neighbor_list(options, neighbors);
}


