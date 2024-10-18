#include "metatensor/torch/atomistic/system.hpp"
#include "utils/Vector.hpp"
#include <variant>

struct PairInfo {
  int part_id_1;
  int part_id_2;
  Utils::Vector3d distance;
};

using Sample = std::array<int32_t, 5>;
using Distances = std::variant<std::vector<std::array<double, 3>>,
                               std::vector<std::array<float, 3>>>;

template <typename PairIterable>
metatensor_torch::TorchTensorBlock
neighbor_list_from_pairs(const metatensor_torch::System &system,
                         const PairIterable &pairs) {
  auto dtype = system->positions().scalar_type();
  auto device = system->positions().device();
  std::vector<Sample> samples;
  Distances distances;

  if (dtype == torch::kFloat64) {
    distances = {std::vector<std::array<double, 3>>()};
  } else if (dtype == torch::kFloat32) {
    distances = {std::vector<std::array<float, 3>>()};
  } else {
    throw std::runtime_error("Unsupported floating point data type");
  }

  for (auto const &pair : pairs) {
    samples.emplace_back(pair.particle_id_1, pair.particle_id_2, 0, 0, 0);
    std::visit([&pair](auto &vec) { vec.push_back(pair.distance); }, distances);
  }

  auto n_pairs = static_cast<int64_t>(samples.size());

  auto samples_tensor = torch::from_blob(
      samples.data(), {n_pairs, 5},
      torch::TensorOptions().dtype(torch::kInt32).device(torch::kCPU));

  auto samples_ptr = torch::make_intrusive<metatensor_torch::LabelsHolder>(
      std::vector<std::string>{"first_atom", "second_atom", "cell_shift_a",
                               "cell_shift_b", "cell_shift_c"},
      samples);

  auto distances_vectors = torch::from_blob(
      std::visit([](auto &vec) { return vec.data(); }, distances),
      {n_pairs, 3, 1}, torch::TensorOptions().dtype(dtype).device(torch::kCPU));

  auto neighbors = torch::make_intrusive<metatensor_torch::TensorBlockHolder>(
      distances_vectors.to(dtype).to(device), samples_ptr->to(device),
      std::vector<metatensor_torch::TorchLabels>{
          metatensor_torch::LabelsHolder::create({"xyz"}, {{0}, {1}, {2}})
              ->to(device),
      },
      metatensor_torch::LabelsHolder::create({"distance"}, {{0}})->to(device));

  return neighbors;
}

void add_neighbor_list_to_system(
    metatensor_torch::System &system,
    const metatensor_torch::TorchTensorBlock &neighbors,
    const metatensor_torch::NeighborListOptions &options,
    bool check_consistency) {
  metatensor_torch::register_autograd_neighbors(system, neighbors,
                                                check_consistency);
  system->add_neighbor_list(options, neighbors);
}
