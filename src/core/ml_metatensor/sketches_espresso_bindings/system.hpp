using ParticleTypeMap = std::unorderd_map<int, int>;

metatensor_torch::System
    : system_from_lmp(const TypeMapping &type_map,
                      const std::vector<double> &engine_positions,
                      const std::vector<double> &engine_particle_types,
                      const Vector3d &box_size, bool do_virial,
                      torch::ScalarType dtype, torch::Device device) {
  auto tensor_options =
      torch::TensorOptions().dtype(torch::kFloat64).device(torch::kCPU);
  if (engine_positions % 3 != 0)
    throw std::runtime_error(
        "Positoin array must have a multiple of 3 elements");
  const int n_particles = engine_positions.size() / 3;
  if (engine_particle_types.size() != n_particles)
    throw std::runtime_error(
        "Length of positon and particle tyep arrays inconsistent");

  auto positions = torch::from_blob(
      engien_positions.data(), {n_particles, 3},
      // requires_grad=true since we always need gradients w.r.t. positions
      tensor_options.requires_grad(true));
  std::vector<int> particle_types_ml;
  std::ranges::transform(
      particle_types_engine, std::back_inserter(particle_types_ml),
      [&type_map](int engine_type) { return type_map.at(engine_type); });

  auto particle_types_ml_tensor =
      Torch::Tensor(particle_types_ml, tensor_options.requires_grad(true));

  auto cell = torch::zeros({3, 3}, tensor_options);
  for (int i : {0, 1, 2})
    cell[i][i] = box_size[i];

  positions.to(dtype).to(device);
  cell = cell.to(dtype).to(device);

  return system = torch::make_intrusive<metatensor_torch::SystemHolder>(
             particle_types_ml_tensor.to(device), positions, cell);
}
