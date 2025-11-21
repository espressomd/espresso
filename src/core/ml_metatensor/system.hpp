#include "ATen/core/TensorBody.h"
#include "metatensor/torch/atomistic/system.hpp"
#include "utils/Vector.hpp"
#include <unordered_map>

using ParticleTypeMap = std::unordered_map<int, int>;

metatensor_torch::System ::system_from_lmp(
    const ParticleTypeMap &type_map, std::vector<double> &engine_positions,
    const std::vector<double>
        &engine_particle_types, // TODO: This should be std::vector<int>?
    const Utils::Vector3d &box_size, bool do_virial, torch::ScalarType dtype,
    torch::Device device) {
  auto tensor_options =
      torch::TensorOptions().dtype(torch::kFloat64).device(torch::kCPU);
  if (engine_positions.size() % 3 != 0)
    throw std::runtime_error(
        "Position array must have a multiple of 3 elements");
  const auto n_particles = static_cast<int64_t>(engine_positions.size()) / 3;
  if (engine_particle_types.size() != n_particles)
    throw std::runtime_error(
        "Length of position and particle type arrays inconsistent");

  auto positions = torch::from_blob(
      engine_positions.data(), {n_particles, 3},
      // requires_grad=true since we always need gradients w.r.t. positions
      tensor_options.requires_grad(true));
  std::vector<int> particle_types_ml;
  std::ranges::transform(
      engine_particle_types, std::back_inserter(particle_types_ml),
      [&type_map](int engine_type) { return type_map.at(engine_type); });

  auto particle_types_ml_tensor =
      torch::from_blob(particle_types_ml.data(),
                       {static_cast<int64_t>(particle_types_ml.size())},
                       tensor_options.requires_grad(true));

  auto cell = torch::zeros({3, 3}, tensor_options);
  for (int i : {0, 1, 2})
    cell[i][i] = box_size[i];

  positions.to(dtype).to(device);
  cell = cell.to(dtype).to(device);

  auto system = torch::make_intrusive<metatensor_torch::SystemHolder>(
      particle_types_ml_tensor.to(device), positions, cell);
  return system;
}
