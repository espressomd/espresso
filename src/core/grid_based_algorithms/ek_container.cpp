#include "ek_container.hpp"

#include "lb_interface.hpp"

#include "errorhandling.hpp"

#include "ek_reactions.hpp"

namespace EK {
EKContainer<EKinWalberlaBase<double>> ek_container;

double get_tau() { return ek_container.get_tau(); }

int get_steps_per_md_step(double md_timestep) {
  return static_cast<int>(std::round(get_tau() / md_timestep));
}

void propagate() {
  // first calculate the charge for the potential, for that get all the
  // field-ids from the ekspecies pass the potential-field-id to the
  // flux-kernels of the eks for this the integrate function has to be split
  // with a public interface to diffusive and advective-flux this should also
  // allow the back-coupling to the LB with a field-id

  if (ek_container.empty()) {
    return;
  }

  ek_container.reset_charge();
  std::for_each(ek_container.begin(), ek_container.end(), [](auto const &ek) {
    ek_container.add_charge(ek->get_density_id(), ek->get_valency());
  });

  ek_container.solve_poisson();

  // TODO: find a way for a proper interface
  const std::size_t velocity_field_id = []() -> std::size_t {
    try {
      return Walberla::get_velocity_field_id();
    } catch (const std::runtime_error &e) {
      return {};
    }
  }();
  const std::size_t force_field_id = []() -> std::size_t {
    try {
      return Walberla::get_force_field_id();
    } catch (const std::runtime_error &e) {
      return {};
    }
  }();

  std::for_each(ek_container.begin(), ek_container.end(),
                [velocity_field_id, force_field_id](auto const &ek) {
                  try {
                    ek->integrate(ek_container.get_potential_field_id(),
                                  velocity_field_id, force_field_id);
                  } catch (const std::runtime_error &e) {
                    runtimeErrorMsg() << e.what();
                  };
                });

  EK::perform_reactions();
}
} // namespace EK
