/*
 * Copyright (C) 2022-2023 The ESPResSo project
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "config/config.hpp"

#include "ek_container.hpp"
#include "ek_reactions.hpp"
#include "errorhandling.hpp"
#include "lb_interface.hpp"
#include "lb_walberla_instance.hpp"

#ifdef WALBERLA
#include <walberla_bridge/electrokinetics/EKContainer.hpp>
#endif // WALBERLA

#include <cmath>

#ifdef WALBERLA
#include <algorithm>
#include <cstddef>
#include <stdexcept>
#endif // WALBERLA

namespace EK {

#ifdef WALBERLA
EKContainer<EKinWalberlaBase> ek_container;
#endif // WALBERLA

double get_tau() {
#ifdef WALBERLA
  return ek_container.get_tau();
#else
  throw NoEKActive();
#endif // WALBERLA
}

int get_steps_per_md_step(double md_timestep) {
  return static_cast<int>(std::round(get_tau() / md_timestep));
}

void propagate() {
#ifdef WALBERLA
  // first calculate the charge for the potential, for that get all the
  // field-ids from the ekspecies pass the potential-field-id to the
  // flux-kernels of the eks for this the integrate function has to be split
  // with a public interface to diffusive and advective-flux this should also
  // allow the back-coupling to the LB with a field-id

  if (ek_container.empty()) {
    return;
  }

  if (!ek_container.is_poisson_solver_set()) {
    runtimeErrorMsg() << "EK requires a Poisson solver.";
    return;
  }

  ek_container.reset_charge();
  std::for_each(ek_container.begin(), ek_container.end(), [](auto const &ek) {
    ek_container.add_charge(ek->get_density_id(), ek->get_valency(),
                            ek->is_double_precision());
  });

  ek_container.solve_poisson();

  auto velocity_field_id = std::size_t{};
  auto force_field_id = std::size_t{};
  try {
    auto const lbf = ::lb_walberla();
    velocity_field_id = lbf->get_velocity_field_id();
    force_field_id = lbf->get_force_field_id();
  } catch (std::runtime_error const &) {
  }

  std::for_each(ek_container.begin(), ek_container.end(),
                [velocity_field_id, force_field_id](auto const &ek) {
                  try {
                    ek->integrate(ek_container.get_potential_field_id(),
                                  velocity_field_id, force_field_id);
                  } catch (std::runtime_error const &e) {
                    runtimeErrorMsg() << e.what();
                  }
                });

  EK::perform_reactions();

  for (auto const &species : ek_container) {
    species->ghost_communication();
  }
#endif // WALBERLA
}

} // namespace EK
