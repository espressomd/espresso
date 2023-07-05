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

#ifdef NPT

#include "communication.hpp"
#include "config/config.hpp"
#include "electrostatics/coulomb.hpp"
#include "errorhandling.hpp"
#include "event.hpp"
#include "integrate.hpp"
#include "magnetostatics/dipoles.hpp"
#include "system/System.hpp"

#include <utils/Vector.hpp>

#include <boost/mpi/collectives/broadcast.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>

static constexpr Utils::Vector3i nptgeom_dir{{1, 2, 4}};

NptIsoParameters nptiso = {};

void synchronize_npt_state() {
  boost::mpi::broadcast(comm_cart, nptiso.p_inst, 0);
  boost::mpi::broadcast(comm_cart, nptiso.p_diff, 0);
  boost::mpi::broadcast(comm_cart, nptiso.volume, 0);
}

void NptIsoParameters::coulomb_dipole_sanity_checks() const {
#ifdef ELECTROSTATICS
  if (dimension < 3 and not cubic_box and System::get_system().coulomb.solver) {
    throw std::runtime_error("If electrostatics is being used you must "
                             "use the cubic box NpT.");
  }
#endif

#ifdef DIPOLES
  if (dimension < 3 and not cubic_box and System::get_system().dipoles.solver) {
    throw std::runtime_error("If magnetostatics is being used you must "
                             "use the cubic box NpT.");
  }
#endif
}

NptIsoParameters::NptIsoParameters(double ext_pressure, double piston,
                                   Utils::Vector<bool, 3> const &rescale,
                                   bool cubic_box)
    : piston{piston}, p_ext{ext_pressure}, cubic_box{cubic_box} {

  if (ext_pressure < 0.0) {
    throw std::runtime_error("The external pressure must be positive");
  }
  if (piston <= 0.0) {
    throw std::runtime_error("The piston mass must be positive");
  }

  inv_piston = ::nptiso.inv_piston;
  volume = ::nptiso.volume;
  p_inst = ::nptiso.p_inst;
  p_diff = ::nptiso.p_diff;
  p_vir = ::nptiso.p_vir;
  p_vel = ::nptiso.p_vel;

  /* set the NpT geometry */
  for (auto const i : {0, 1, 2}) {
    if (rescale[i]) {
      geometry |= ::nptgeom_dir[i];
      dimension += 1;
      non_const_dim = i;
    }
  }

  if (dimension == 0 or non_const_dim == -1) {
    throw std::runtime_error(
        "You must enable at least one of the x y z components "
        "as fluctuating dimension(s) for box length motion");
  }

  coulomb_dipole_sanity_checks();
}

Utils::Vector<bool, 3> NptIsoParameters::get_direction() const {
  return {static_cast<bool>(geometry & ::nptgeom_dir[0]),
          static_cast<bool>(geometry & ::nptgeom_dir[1]),
          static_cast<bool>(geometry & ::nptgeom_dir[2])};
}

void npt_ensemble_init(const BoxGeometry &box) {
  if (integ_switch == INTEG_METHOD_NPT_ISO) {
    /* prepare NpT-integration */
    nptiso.inv_piston = 1. / nptiso.piston;
    nptiso.volume = pow(box.length()[nptiso.non_const_dim], nptiso.dimension);
    if (recalc_forces) {
      nptiso.p_inst = 0.0;
      nptiso.p_vir = Utils::Vector3d{};
      nptiso.p_vel = Utils::Vector3d{};
    }
  }
}

void integrator_npt_sanity_checks() {
  if (integ_switch == INTEG_METHOD_NPT_ISO) {
    try {
      nptiso.coulomb_dipole_sanity_checks();
    } catch (std::runtime_error const &err) {
      runtimeErrorMsg() << err.what();
    }
  }
}

/** reset virial part of instantaneous pressure */
void npt_reset_instantaneous_virials() {
  if (integ_switch == INTEG_METHOD_NPT_ISO)
    nptiso.p_vir = Utils::Vector3d{};
}

void npt_add_virial_contribution(double energy) {
  if (integ_switch == INTEG_METHOD_NPT_ISO) {
    nptiso.p_vir[0] += energy;
  }
}

void npt_add_virial_contribution(const Utils::Vector3d &force,
                                 const Utils::Vector3d &d) {
  if (integ_switch == INTEG_METHOD_NPT_ISO) {
    nptiso.p_vir += hadamard_product(force, d);
  }
}
#endif // NPT
