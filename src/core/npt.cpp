/*
 * Copyright (C) 2010-2022 The ESPResSo project
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
#include "npt.hpp"

#ifdef ESPRESSO_NPT

#include "BoxGeometry.hpp"
#include "PropagationMode.hpp"
#include "cell_system/CellStructure.hpp"
#include "communication.hpp"
#include "config/config.hpp"
#include "electrostatics/coulomb.hpp"
#include "errorhandling.hpp"
#include "integrators/Propagation.hpp"
#include "magnetostatics/dipoles.hpp"
#include "system/System.hpp"

#include <utils/Vector.hpp>
#include <utils/mpi/gather_buffer.hpp>

#include <boost/mpi/collectives.hpp>

#include <algorithm>
#include <cmath>
#include <functional>
#include <stdexcept>

void System::System::synchronize_npt_state() {
  boost::mpi::broadcast(comm_cart, nptiso->p_epsilon, 0);
  boost::mpi::broadcast(comm_cart, nptiso->volume, 0);
  boost::mpi::broadcast(comm_cart, npt_inst_pressure->p_inst, 0);
}

void NptIsoParameters::coulomb_dipole_sanity_checks(
    System::System const &system) const {
#ifdef ESPRESSO_ELECTROSTATICS
  if (dimension < 3 and not cubic_box and system.coulomb.impl->solver) {
    throw std::runtime_error("If electrostatics is being used you must "
                             "use the cubic box NpT.");
  }
#endif
#ifdef ESPRESSO_DIPOLES
  if (dimension < 3 and not cubic_box and system.dipoles.impl->solver) {
    throw std::runtime_error("If magnetostatics is being used you must "
                             "use the cubic box NpT.");
  }
#endif
}

NptIsoParameters::NptIsoParameters(double ext_pressure, double piston,
                                   Utils::Vector<bool, 3> const &rescale,
                                   bool cubic_box)
    : piston{piston}, inv_piston{piston}, p_ext{ext_pressure},
      cubic_box{cubic_box} {

  if (ext_pressure < 0.0) {
    throw std::runtime_error("The external pressure must be positive");
  }
  if (piston <= 0.0) {
    throw std::runtime_error("The piston mass must be positive");
  }

  /* set the NpT geometry */
  for (auto const i : {0u, 1u, 2u}) {
    if (rescale[i]) {
      geometry |= nptgeom_dir[i];
      dimension += 1;
      non_const_dim = static_cast<int>(i);
    }
  }

  if (dimension == 0 or non_const_dim == -1) {
    throw std::runtime_error(
        "You must enable at least one of the x y z components "
        "as fluctuating dimension(s) for box length motion");
  }
}

Utils::Vector<bool, 3> NptIsoParameters::get_direction() const {
  return {static_cast<bool>(geometry & nptgeom_dir[0]),
          static_cast<bool>(geometry & nptgeom_dir[1]),
          static_cast<bool>(geometry & nptgeom_dir[2])};
}

void System::System::npt_ensemble_init(bool recalc_forces) {
  nptiso->inv_piston = 1. / nptiso->piston;
  nptiso->volume =
      std::pow(box_geo->length()[nptiso->non_const_dim], nptiso->dimension);
  if (recalc_forces) {
    npt_inst_pressure->p_inst = Utils::Vector2d{};
    npt_inst_pressure->p_vir = Utils::Vector3d{};
    npt_inst_pressure->p_vel = Utils::Vector3d{};
  }

  auto const particle_number = cell_structure->local_particles().size();
  nptiso->particle_number =
      boost::mpi::all_reduce(::comm_cart, particle_number, std::plus<>());

  auto const dt = get_time_step();
  nptiso->half_dt_inv_piston = 0.5 * dt * nptiso->inv_piston;
  nptiso->half_dt_inv_piston_and_Nf = -nptiso->half_dt_inv_piston;
  if (particle_number > 1) {
    nptiso->half_dt_inv_piston_and_Nf *=
        (1. + 1. / static_cast<double>(particle_number - 1));
  }

  auto &mass_list = nptiso->mass_list;
  mass_list.clear();
  for (auto &p : cell_structure->local_particles()) {
    mass_list.emplace_back(p.mass());
  }
  mass_list.erase(std::ranges::unique(mass_list).begin(), mass_list.end());
  Utils::Mpi::gather_buffer(mass_list, ::comm_cart);
  if (::this_node == 0) {
    mass_list.erase(std::ranges::unique(mass_list).begin(), mass_list.end());
    std::ranges::sort(mass_list);
  }
  boost::mpi::broadcast(::comm_cart, mass_list, 0);
}

void System::System::npt_add_virial_contribution(double energy) {
  if (has_npt_enabled()) {
    npt_inst_pressure->p_vir[0] += energy;
  }
}
#endif // ESPRESSO_NPT
