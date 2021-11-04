/*
 * Copyright (C) 2010-2019 The ESPResSo project
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
#include "config.hpp"
#include "electrostatics_magnetostatics/coulomb.hpp"
#include "electrostatics_magnetostatics/dipole.hpp"
#include "errorhandling.hpp"
#include "event.hpp"
#include "integrate.hpp"

#include <utils/Vector.hpp>

#include <boost/mpi/collectives/broadcast.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>

NptIsoParameters nptiso = {0.0,
                           0.0,
                           0.0,
                           0.0,
                           0.0,
                           0.0,
                           {0.0, 0.0, 0.0},
                           {0.0, 0.0, 0.0},
                           0,
                           {NPTGEOM_XDIR, NPTGEOM_YDIR, NPTGEOM_ZDIR},
                           0,
                           false,
                           0};

void synchronize_npt_state() {
  boost::mpi::broadcast(comm_cart, nptiso.p_inst, 0);
  boost::mpi::broadcast(comm_cart, nptiso.p_diff, 0);
  boost::mpi::broadcast(comm_cart, nptiso.volume, 0);
}

void mpi_bcast_nptiso_geom_barostat_local() {
  boost::mpi::broadcast(comm_cart, nptiso.geometry, 0);
  boost::mpi::broadcast(comm_cart, nptiso.dimension, 0);
  boost::mpi::broadcast(comm_cart, nptiso.cubic_box, 0);
  boost::mpi::broadcast(comm_cart, nptiso.non_const_dim, 0);
  boost::mpi::broadcast(comm_cart, nptiso.p_ext, 0);
  boost::mpi::broadcast(comm_cart, nptiso.piston, 0);
  on_thermostat_param_change();
}

REGISTER_CALLBACK(mpi_bcast_nptiso_geom_barostat_local)

/** Broadcast nptiso geometry and barostat parameters to all nodes. */
void mpi_bcast_nptiso_geom_barostat() {
  mpi_call_all(mpi_bcast_nptiso_geom_barostat_local);
}

void integrator_npt_coulomb_dipole_sanity_checks(
    NptIsoParameters const &params) {
#ifdef ELECTROSTATICS
  if (params.dimension < 3 && !params.cubic_box && coulomb.prefactor > 0) {
    throw std::runtime_error("If electrostatics is being used you must "
                             "use the cubic box NpT.");
  }
#endif

#ifdef DIPOLES
  if (params.dimension < 3 && !params.cubic_box && dipole.prefactor > 0) {
    throw std::runtime_error("If magnetostatics is being used you must "
                             "use the cubic box NpT.");
  }
#endif

#if defined(ELECTROSTATICS) && defined(CUDA)
  if (coulomb.method == COULOMB_P3M_GPU) {
    throw std::runtime_error("NpT virial cannot be calculated on P3M GPU");
  }
#endif
}

void nptiso_init(double ext_pressure, double piston, bool xdir_rescale,
                 bool ydir_rescale, bool zdir_rescale, bool cubic_box) {

  if (ext_pressure < 0.0) {
    throw std::runtime_error("The external pressure must be positive.");
  }
  if (piston <= 0.0) {
    throw std::runtime_error("The piston mass must be positive.");
  }

  NptIsoParameters new_nptiso = {piston,
                                 nptiso.inv_piston,
                                 nptiso.volume,
                                 ext_pressure,
                                 nptiso.p_inst,
                                 nptiso.p_diff,
                                 nptiso.p_vir,
                                 nptiso.p_vel,
                                 0,
                                 nptiso.nptgeom_dir,
                                 0,
                                 cubic_box,
                                 -1};

  /* set the NpT geometry */
  if (xdir_rescale) {
    new_nptiso.geometry |= NPTGEOM_XDIR;
    new_nptiso.dimension += 1;
    new_nptiso.non_const_dim = 0;
  }
  if (ydir_rescale) {
    new_nptiso.geometry |= NPTGEOM_YDIR;
    new_nptiso.dimension += 1;
    new_nptiso.non_const_dim = 1;
  }
  if (zdir_rescale) {
    new_nptiso.geometry |= NPTGEOM_ZDIR;
    new_nptiso.dimension += 1;
    new_nptiso.non_const_dim = 2;
  }

  if (new_nptiso.dimension == 0 || new_nptiso.non_const_dim == -1) {
    throw std::runtime_error(
        "You must enable at least one of the x y z components "
        "as fluctuating dimension(s) for box length motion!");
  }

  integrator_npt_coulomb_dipole_sanity_checks(new_nptiso);

  nptiso = new_nptiso;

  mpi_bcast_nptiso_geom_barostat();
}

void npt_ensemble_init(const BoxGeometry &box) {
  if (integ_switch == INTEG_METHOD_NPT_ISO) {
    /* prepare NpT-integration */
    nptiso.inv_piston = 1 / nptiso.piston;
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
      integrator_npt_coulomb_dipole_sanity_checks(nptiso);
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
