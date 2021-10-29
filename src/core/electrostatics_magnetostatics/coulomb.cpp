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

#include "config.hpp"

#ifdef ELECTROSTATICS

#include "electrostatics_magnetostatics/coulomb.hpp"

#include "ParticleRange.hpp"
#include "cells.hpp"
#include "communication.hpp"
#include "electrostatics_magnetostatics/common.hpp"
#include "electrostatics_magnetostatics/debye_hueckel.hpp"
#include "electrostatics_magnetostatics/elc.hpp"
#include "electrostatics_magnetostatics/icc.hpp"
#include "electrostatics_magnetostatics/mmm1d.hpp"
#include "electrostatics_magnetostatics/p3m-common.hpp"
#include "electrostatics_magnetostatics/p3m.hpp"
#include "electrostatics_magnetostatics/p3m_gpu.hpp"
#include "electrostatics_magnetostatics/reaction_field.hpp"
#include "electrostatics_magnetostatics/scafacos.hpp"
#include "errorhandling.hpp"
#include "grid_based_algorithms/electrokinetics.hpp"
#include "integrate.hpp"
#include "npt.hpp"

#include <utils/Vector.hpp>
#include <utils/constants.hpp>

#include <mpi.h>

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <limits>
#include <stdexcept>

Coulomb_parameters coulomb;

namespace Coulomb {
Utils::Vector9d calc_pressure_long_range(const ParticleRange &particles) {
  switch (coulomb.method) {
#ifdef P3M
  case COULOMB_ELC_P3M:
    fprintf(stderr,
            "WARNING: pressure calculated, but ELC pressure not implemented\n");
    break;
  case COULOMB_P3M_GPU:
    fprintf(
        stderr,
        "WARNING: pressure calculated, but GPU P3M pressure not implemented\n");
    break;
  case COULOMB_P3M: {
    p3m_charge_assign(particles);
    return p3m_calc_kspace_pressure_tensor();
  }
#endif
  case COULOMB_MMM1D:
  case COULOMB_MMM1D_GPU:
    fprintf(
        stderr,
        "WARNING: pressure calculated, but MMM1D pressure not implemented\n");
    break;
  default:
    break;
  }
  return {};
}

bool sanity_checks() {
  bool failed = false;
  switch (coulomb.method) {
  case COULOMB_MMM1D:
    if (MMM1D_sanity_checks())
      failed = true;
    break;
#ifdef P3M
  case COULOMB_ELC_P3M:
    try {
      ELC_sanity_checks(elc_params);
    } catch (std::runtime_error const &err) {
      runtimeErrorMsg() << err.what();
      failed = true;
    }
    // fall through
  case COULOMB_P3M_GPU:
  case COULOMB_P3M:
    if (p3m_sanity_checks())
      failed = true;
    break;
#endif
  default:
    break;
  }
  return failed;
}

double cutoff(const Utils::Vector3d &box_l) {
  switch (coulomb.method) {
  case COULOMB_MMM1D:
    return std::numeric_limits<double>::infinity();
#ifdef P3M
  case COULOMB_ELC_P3M:
    return std::max(elc_params.space_layer, p3m.params.r_cut_iL * box_l[0]);
  case COULOMB_P3M_GPU:
  case COULOMB_P3M:
    /* do not use precalculated r_cut here, might not be set yet */
    return p3m.params.r_cut_iL * box_l[0];
#endif
  case COULOMB_DH:
    return dh_params.r_cut;
  case COULOMB_RF:
    return rf_params.r_cut;
#ifdef SCAFACOS
  case COULOMB_SCAFACOS:
    return Scafacos::fcs_coulomb()->get_r_cut();
#endif
  default:
    return -1.0;
  }
}

void deactivate() {
  switch (coulomb.method) {
#ifdef P3M
  case COULOMB_ELC_P3M:
  case COULOMB_P3M_GPU:
  case COULOMB_P3M:
    break;
#endif
  case COULOMB_DH:
    dh_params = {};
    break;
  case COULOMB_RF:
    rf_params = {};
    break;
  case COULOMB_MMM1D:
    mmm1d_params.maxPWerror = 1e40;
    break;
  default:
    break;
  }
}

void update_dependent_particles() {
  icc_iteration(cell_structure.local_particles(),
                cell_structure.ghost_particles());
}

void on_observable_calc() {
  switch (coulomb.method) {
#ifdef P3M
  case COULOMB_ELC_P3M:
  case COULOMB_P3M_GPU:
  case COULOMB_P3M:
    p3m_count_charged_particles();
    break;
#endif
  default:
    break;
  }
}

void on_coulomb_change() {
  switch (coulomb.method) {
  case COULOMB_DH:
    break;
#ifdef P3M
#ifdef CUDA
  case COULOMB_P3M_GPU:
    if (this_node == 0) {
      try {
        p3m_gpu_init(p3m.params.cao, p3m.params.mesh, p3m.params.alpha);
      } catch (std::runtime_error const &err) {
        runtimeErrorMsg() << err.what();
      }
    }
    break;
#endif
  case COULOMB_ELC_P3M:
    ELC_init();
    // fall through
  case COULOMB_P3M:
    p3m_init();
    break;
#endif
  case COULOMB_MMM1D:
    MMM1D_init();
    break;
  default:
    break;
  }
}

void on_boxl_change() {
  switch (coulomb.method) {
#ifdef P3M
  case COULOMB_ELC_P3M:
    ELC_init();
    // fall through
  case COULOMB_P3M_GPU:
  case COULOMB_P3M:
    p3m_scaleby_box_l();
    break;
#endif
  case COULOMB_MMM1D:
    MMM1D_init();
    break;
#ifdef SCAFACOS
  case COULOMB_SCAFACOS:
    Scafacos::fcs_coulomb()->update_system_params();
    break;
#endif
  default:
    break;
  }
}

void init() {
  switch (coulomb.method) {
  case COULOMB_DH:
    break;
#ifdef P3M
  case COULOMB_ELC_P3M:
    ELC_init();
    // fall through
  case COULOMB_P3M:
    p3m_init();
    break;
  case COULOMB_P3M_GPU:
    break;
#endif
  case COULOMB_MMM1D:
    MMM1D_init();
    break;
  default:
    break;
  }
}

void calc_long_range_force(const ParticleRange &particles) {
  switch (coulomb.method) {
#ifdef P3M
  case COULOMB_ELC_P3M:
    if (elc_params.dielectric_contrast_on) {
      ELC_P3M_modify_p3m_sums_both(particles);
      ELC_p3m_charge_assign_both(particles);
      ELC_P3M_self_forces(particles);
    } else
      p3m_charge_assign(particles);

    p3m_calc_kspace_forces(true, false, particles);

    if (elc_params.dielectric_contrast_on)
      ELC_P3M_restore_p3m_sums(particles);

    ELC_add_force(particles);

    break;
#endif
#ifdef CUDA
  case COULOMB_P3M_GPU:
    if (this_node == 0) {
      p3m_gpu_add_farfield_force();
    }
    /* there is no NPT handling here as long as we cannot compute energies.
       This is checked in integrator_npt_sanity_checks() when integration
       starts. */
    break;
#endif
#ifdef P3M
  case COULOMB_P3M:
    p3m_charge_assign(particles);
#ifdef NPT
    if (integ_switch == INTEG_METHOD_NPT_ISO) {
      auto const energy = p3m_calc_kspace_forces(true, true, particles);
      npt_add_virial_contribution(energy);
    } else
#endif
      p3m_calc_kspace_forces(true, false, particles);
    break;
#endif
#ifdef SCAFACOS
  case COULOMB_SCAFACOS:
    Scafacos::fcs_coulomb()->add_long_range_force();
    break;
#endif
  default:
    break;
  }

  /* Add fields from EK if enabled */
#ifdef ELECTROKINETICS
  if (this_node == 0) {
    ek_calculate_electrostatic_coupling();
  }
#endif
}

double calc_energy_long_range(const ParticleRange &particles) {
  double energy = 0.0;
  switch (coulomb.method) {
#ifdef P3M
  case COULOMB_P3M_GPU:
    runtimeWarningMsg()
        << "long range energy calculation not implemented for GPU P3M";
    break;
  case COULOMB_P3M:
    p3m_charge_assign(particles);
    energy = p3m_calc_kspace_forces(false, true, particles);
    break;
  case COULOMB_ELC_P3M:
    // assign the original charges first
    // they may not have been assigned yet
    p3m_charge_assign(particles);
    if (!elc_params.dielectric_contrast_on)
      energy = p3m_calc_kspace_forces(false, true, particles);
    else {
      energy = 0.5 * p3m_calc_kspace_forces(false, true, particles);
      energy += 0.5 * coulomb.prefactor *
                ELC_P3M_dielectric_layers_energy_self(particles);

      // assign both original and image charges now
      ELC_p3m_charge_assign_both(particles);
      ELC_P3M_modify_p3m_sums_both(particles);

      energy += 0.5 * p3m_calc_kspace_forces(false, true, particles);

      // assign only the image charges now
      ELC_p3m_charge_assign_image(particles);
      ELC_P3M_modify_p3m_sums_image(particles);

      energy -= 0.5 * p3m_calc_kspace_forces(false, true, particles);

      // restore modified sums
      ELC_P3M_restore_p3m_sums(particles);
    }
    energy += ELC_energy(particles);
    break;
#endif
#ifdef SCAFACOS
  case COULOMB_SCAFACOS:
    energy += Scafacos::fcs_coulomb()->long_range_energy();
    break;
#endif
  default:
    break;
  }
  return energy;
}

int icc_sanity_check() {
  switch (coulomb.method) {
#ifdef P3M
  case COULOMB_ELC_P3M: {
    if (elc_params.dielectric_contrast_on) {
      runtimeErrorMsg() << "ICC conflicts with ELC dielectric contrast";
      return 1;
    }
    break;
  }
#endif
  case COULOMB_DH: {
    runtimeErrorMsg() << "ICC does not work with Debye-Hueckel.";
    return 1;
  }
  case COULOMB_RF: {
    runtimeErrorMsg() << "ICC does not work with COULOMB_RF.";
    return 1;
  }
  default:
    break;
  }

#ifdef NPT
  if (integ_switch == INTEG_METHOD_NPT_ISO) {
    runtimeErrorMsg() << "ICC does not work in the NPT ensemble";
    return 1;
  }
#endif

  return 0;
}

void bcast_coulomb_params() {
  switch (coulomb.method) {
  case COULOMB_NONE:
  case COULOMB_SCAFACOS:
    break;
#ifdef P3M
  case COULOMB_ELC_P3M:
    MPI_Bcast(&elc_params, sizeof(ELCParameters), MPI_BYTE, 0, comm_cart);
    // fall through
  case COULOMB_P3M_GPU:
  case COULOMB_P3M:
    MPI_Bcast(&p3m.params, sizeof(P3MParameters), MPI_BYTE, 0, comm_cart);
    break;
#endif
  case COULOMB_DH:
    MPI_Bcast(&dh_params, sizeof(DebyeHueckelParameters), MPI_BYTE, 0,
              comm_cart);
    break;
  case COULOMB_MMM1D:
  case COULOMB_MMM1D_GPU:
    MPI_Bcast(&mmm1d_params, sizeof(MMM1DParameters), MPI_BYTE, 0, comm_cart);
    break;
  case COULOMB_RF:
    MPI_Bcast(&rf_params, sizeof(ReactionFieldParameters), MPI_BYTE, 0,
              comm_cart);
    break;
  default:
    break;
  }
}

void set_prefactor(double prefactor) {
  if (prefactor < 0.0) {
    throw std::domain_error("Coulomb prefactor has to be >= 0");
  }

  coulomb.prefactor = prefactor;
  mpi_bcast_coulomb_params();
}

void deactivate_method() {
  coulomb.prefactor = 0;

  Coulomb::deactivate();

  mpi_bcast_coulomb_params();
  coulomb.method = COULOMB_NONE;
  mpi_bcast_coulomb_params();
}
} // namespace Coulomb
#endif // ELECTROSTATICS
