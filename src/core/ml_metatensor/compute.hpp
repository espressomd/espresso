#include "metatensor/torch/atomistic/system.hpp"
#include <metatensor/torch/atomistic/model.hpp>
#include <metatensor/torch/tensor.hpp>

metatensor_torch::TensorMapHolder
run_model(metatensor_torch::System &system, int64_t n_particles,
          const metatensor_torch::ModelEvaluationOptions evaluation_options,
          torch::Dtype dtype, torch::Device device, bool check_consistency) {

  // only run the calculation for atoms actually in the current domain
  auto options = torch::TensorOptions().dtype(torch::kInt32);
  auto selected_atoms_values = torch::zeros({n_particles, 2}, options);

  for (int i = 0; i < n_particles; i++) {
    selected_atoms_values[i][0] = 0;
    selected_atoms_values[i][1] = i;
  }
  auto selected_atoms = torch::make_intrusive<metatensor_torch::LabelsHolder>(
      std::vector<std::string>{"system", "atom"}, selected_atoms_values);
  evaluation_options->set_selected_atoms(selected_atoms->to(device));

  torch::IValue result_ivalue;
  model->forward({std::vector<metatensor_torch::System>{system},
                  evaluation_options, check_consistency});

  auto result = result_ivalue.toGenericDict();
  auto energy =
      result.at("energy").toCustomClass<metatensor_torch::TensorMapHolder>();
  auto energy_tensor =
      metatensor_torch::TensorMapHolder::block_by_id(energy, 0);
}

double get_energy(metatensor_torch::TensorMapHolder &energy,
                  bool energy_is_per_atom) {
  auto energy_block = metatensor_torch::TensorMapHolder::block_by_id(energy, 0);
  auto energy_tensor = energy_block->values();
  auto energy_detached =
      energy_tensor.detach().to(torch::kCPU).to(torch::kFloat64);
  auto energy_samples = energy_block->samples();

  // store the energy returned by the model
  torch::Tensor global_energy;
  if (energy_is_per_atom) {
    assert(energy_samples->size() == 2);
    assert(energy_samples->names()[0] == "system");
    assert(energy_samples->names()[1] == "atom");

    auto samples_values = energy_samples->values().to(torch::kCPU);
    auto samples = samples_values.accessor<int32_t, 2>();

    //        int n_atoms = selected_atoms_values.sizes();
    //        assert(samples_values.sizes() == selected_atoms_values.sizes());

    auto energies = energy_detached.accessor<double, 2>();
    global_energy = energy_detached.sum(0);
    assert(energy_detached.sizes() == std::vector<int64_t>({1}));
  } else {
    assert(energy_samples->size() == 1);
    assert(energy_samples->names()[0] == "system");

    assert(energy_detached.sizes() == std::vector<int64_t>({1, 1}));
    global_energy = energy_detached.reshape({1});
  }

  return global_energy.item<double>();
}

torch::Tensor get_forces(metatensor::TensorMap &energy,
                         metatensor_torch::System &system) {
  // reset gradients to zero before calling backward
  system->positions().mutable_grad() = torch::Tensor();

  auto energy_block = metatensor_torch::TensorMapHolder::block_by_id(energy, 0);
  auto energy_tensor = energy_block->values();

  // compute forces/virial with backward propagation
  energy_tensor.backward(-torch::ones_like(energy_tensor));
  auto forces_tensor = system->positions().grad();
  assert(forces_tensor.is_cpu() &&
         forces_tensor.scalar_type() == torch::kFloat64);
  return forces_tensor;
}
