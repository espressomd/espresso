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

#ifdef WALBERLA

#include "ek/EKReactions.hpp"
#include "ek/EKWalberla.hpp"
#include "errorhandling.hpp"
#include "grid.hpp"
#include "integrate.hpp"
#include "lb/Implementation.hpp"
#include "lb/LBWalberla.hpp"
#include "system/System.hpp"

#include <walberla_bridge/electrokinetics/EKContainer.hpp>
#include <walberla_bridge/electrokinetics/EKinWalberlaBase.hpp>
#include <walberla_bridge/electrokinetics/reactions/EKReactionBase.hpp>
#include <walberla_bridge/lattice_boltzmann/LBWalberlaBase.hpp>

#include <cstddef>
#include <variant>

namespace EK {

double EKWalberla::get_tau() const { return ek_container->get_tau(); }

bool EKWalberla::is_ready_for_propagation() const noexcept {
  return not ek_container->empty();
}

struct FieldsConnector {
  std::size_t velocity_field_id{};
  std::size_t force_field_id{};
  void operator()(LB::Solver::Implementation const &impl) {
    using lb_value_type = std::shared_ptr<LB::LBWalberla>;
    if (impl.solver.has_value()) {
      if (auto const *ptr = std::get_if<lb_value_type>(&(*impl.solver))) {
        auto const &instance = **ptr;
        velocity_field_id = instance.lb_fluid->get_velocity_field_id();
        force_field_id = instance.lb_fluid->get_force_field_id();
      }
    }
  }
};

void EKWalberla::propagate() {
  // first calculate the charge for the potential, for that get all the
  // field-ids from the ekspecies pass the potential-field-id to the
  // flux-kernels of the eks for this the integrate function has to be split
  // with a public interface to diffusive and advective-flux this should also
  // allow the back-coupling to the LB with a field-id

  if (ek_container->empty()) {
    return;
  }

  ek_container->reset_charge();
  for (auto const &ek_species : *ek_container) {
    ek_container->add_charge(ek_species->get_density_id(),
                             ek_species->get_valency(),
                             ek_species->is_double_precision());
  }
  ek_container->solve_poisson();

  FieldsConnector connector{};
  System::get_system().lb.connect(connector);
  for (auto const &ek_species : *ek_container) {
    try {
      ek_species->integrate(ek_container->get_potential_field_id(),
                            connector.velocity_field_id,
                            connector.force_field_id);
    } catch (std::runtime_error const &e) {
      runtimeErrorMsg() << e.what();
    }
  }

  perform_reactions();

  for (auto const &ek_species : *ek_container) {
    ek_species->ghost_communication();
  }
}

void EKWalberla::perform_reactions() {
  for (auto const &ek_reaction : *ek_reactions) {
    ek_reaction->perform_reaction();
  }
}

void EKWalberla::veto_time_step(double time_step) const {
  walberla_tau_sanity_checks("EK", ek_container->get_tau(), time_step);
}

void EKWalberla::sanity_checks() const {
  auto const &lattice = ek_container->get_lattice();
  auto const agrid = ::box_geo.length()[0] / lattice.get_grid_dimensions()[0];
  auto [my_left, my_right] = lattice.get_local_domain();
  my_left *= agrid;
  my_right *= agrid;
  walberla_agrid_sanity_checks("EK", my_left, my_right, agrid);
  // EK time step and MD time step must agree
  walberla_tau_sanity_checks("EK", ek_container->get_tau());
}

} // namespace EK

#endif // WALBERLA
