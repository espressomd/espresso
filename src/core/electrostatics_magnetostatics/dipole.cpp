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

#include "electrostatics_magnetostatics/dipole.hpp"

// Real space cutoff of long range methods
double dipolar_cutoff;

#ifdef DIPOLES

#include "actor/DipolarBarnesHut.hpp"
#include "actor/DipolarDirectSum.hpp"
#include "electrostatics_magnetostatics/magnetic_non_p3m_methods.hpp"
#include "electrostatics_magnetostatics/mdlc_correction.hpp"
#include "electrostatics_magnetostatics/p3m-dipolar.hpp"
#include "electrostatics_magnetostatics/scafacos.hpp"

#include "communication.hpp"
#include "errorhandling.hpp"
#include "grid.hpp"
#include "integrate.hpp"
#include "npt.hpp"

#include <utils/constants.hpp>

#include <boost/mpi/collectives.hpp>
#include <cstdio>

Dipole_parameters dipole = {
    0.0,
    DIPOLAR_NONE,
};

namespace Dipole {
int pressure_n() { return 0; }

void calc_pressure_long_range() {
  switch (dipole.method) {
  case DIPOLAR_NONE:
    return;
  default:
    runtimeWarningMsg()
        << "WARNING: pressure calculated, but pressure not implemented.\n";
    return;
  }
}

void nonbonded_sanity_check(int &state) {
#ifdef DP3M
  switch (dipole.method) {
  case DIPOLAR_MDLC_P3M:
    if (mdlc_sanity_checks())
      state = 0; // fall through
  case DIPOLAR_P3M:
    if (dp3m_sanity_checks(node_grid))
      state = 0;
    break;
  case DIPOLAR_MDLC_DS:
    if (mdlc_sanity_checks())
      state = 0; // fall through
  case DIPOLAR_DS:
    if (magnetic_dipolar_direct_sum_sanity_checks())
      state = 0;
    break;
  default:
    break;
  }
#endif
}

double cutoff(const Utils::Vector3d &box_l) {
  switch (dipole.method) {
#ifdef DP3M
  case DIPOLAR_MDLC_P3M:
    // fall through
  case DIPOLAR_P3M: {
    /* do not use precalculated r_cut here, might not be set yet */
    return dp3m.params.r_cut_iL * box_l[0];
  }
#endif /*ifdef DP3M */
  default:
    return -1.;
  }
}

void on_observable_calc() {
  switch (dipole.method) {
#ifdef DP3M
  case DIPOLAR_MDLC_P3M:
    // fall through
  case DIPOLAR_P3M:
    dp3m_count_magnetic_particles();
    break;
#endif
  default:
    break;
  }
}

void on_coulomb_change() {
  switch (dipole.method) {
#ifdef DP3M
  case DIPOLAR_MDLC_P3M:
    // fall through
  case DIPOLAR_P3M:
    dp3m_init();
    break;
#endif
  default:
    break;
  }
}

void on_boxl_change() {
  switch (dipole.method) {
#ifdef DP3M
  case DIPOLAR_MDLC_P3M:
    // fall through
  case DIPOLAR_P3M:
    dp3m_scaleby_box_l();
    break;
#endif
#ifdef SCAFACOS
  case DIPOLAR_SCAFACOS:
    Scafacos::update_system_params();
    break;
#endif
  default:
    break;
  }
}

void init() {
  switch (dipole.method) {
#ifdef DP3M
  case DIPOLAR_MDLC_P3M:
    // fall through
  case DIPOLAR_P3M:
    dp3m_init();
    break;
#endif
  default:
    break;
  }
}

void calc_long_range_force(const ParticleRange &particles) {
  switch (dipole.method) {
#ifdef DP3M
  case DIPOLAR_MDLC_P3M:
    add_mdlc_force_corrections(particles);
    // fall through
  case DIPOLAR_P3M:
    dp3m_dipole_assign(particles);
#ifdef NPT
    if (integ_switch == INTEG_METHOD_NPT_ISO) {
      nptiso.p_vir[0] += dp3m_calc_kspace_forces(true, true, particles);
      fprintf(stderr, "dipolar_P3M at this moment is added to p_vir[0]\n");
    } else
#endif
      dp3m_calc_kspace_forces(true, false, particles);

    break;
#endif
  case DIPOLAR_ALL_WITH_ALL_AND_NO_REPLICA:
    dawaanr_calculations(true, false, particles);
    break;
#ifdef DP3M
  case DIPOLAR_MDLC_DS:
    add_mdlc_force_corrections(particles);
    // fall through
#endif
  case DIPOLAR_DS:
    magnetic_dipolar_direct_sum_calculations(true, false, particles);
    break;
  case DIPOLAR_DS_GPU: // NOLINT(bugprone-branch-clone)
    // do nothing: it's an actor
    break;
#ifdef DIPOLAR_BARNES_HUT
  case DIPOLAR_BH_GPU:
    // do nothing: it's an actor
    break;
#endif
#ifdef SCAFACOS_DIPOLES
  case DIPOLAR_SCAFACOS:
    assert(Scafacos::dipolar());
    Scafacos::add_long_range_force();
#endif
  case DIPOLAR_NONE:
    break;
  default:
    runtimeErrorMsg() << "unknown dipolar method";
    break;
  }
}

void calc_energy_long_range(Observable_stat &energy,
                            const ParticleRange &particles) {
  switch (dipole.method) {
#ifdef DP3M
  case DIPOLAR_P3M:
    dp3m_dipole_assign(particles);
    energy.dipolar[1] = dp3m_calc_kspace_forces(false, true, particles);
    break;
  case DIPOLAR_MDLC_P3M:
    dp3m_dipole_assign(particles);
    energy.dipolar[1] = dp3m_calc_kspace_forces(false, true, particles);
    energy.dipolar[2] = add_mdlc_energy_corrections(particles);
    break;
#endif
  case DIPOLAR_ALL_WITH_ALL_AND_NO_REPLICA:
    energy.dipolar[1] = dawaanr_calculations(false, true, particles);
    break;
#ifdef DP3M
  case DIPOLAR_MDLC_DS:
    energy.dipolar[1] =
        magnetic_dipolar_direct_sum_calculations(false, true, particles);
    energy.dipolar[2] = add_mdlc_energy_corrections(particles);
    break;
#endif
  case DIPOLAR_DS:
    energy.dipolar[1] =
        magnetic_dipolar_direct_sum_calculations(false, true, particles);
    break;
  case DIPOLAR_DS_GPU: // NOLINT(bugprone-branch-clone)
    // do nothing: it's an actor
    break;
#ifdef DIPOLAR_BARNES_HUT
  case DIPOLAR_BH_GPU:
    // do nothing: it's an actor.
    break;
#endif
#ifdef SCAFACOS_DIPOLES
  case DIPOLAR_SCAFACOS:
    assert(Scafacos::dipolar());
    energy.dipolar[1] = Scafacos::long_range_energy();
#endif
  case DIPOLAR_NONE:
    break;
  default:
    runtimeErrorMsg()
        << "energy calculation not implemented for dipolar method.";
    break;
  }
}

void energy_n(int &n_dipolar) {
  switch (dipole.method) {
  case DIPOLAR_NONE:
    n_dipolar = 1; // because there may be an external magnetic field
    break;
  case DIPOLAR_P3M:
  case DIPOLAR_ALL_WITH_ALL_AND_NO_REPLICA:
  case DIPOLAR_DS:
  case DIPOLAR_DS_GPU:
#ifdef DIPOLAR_BARNES_HUT
  case DIPOLAR_BH_GPU:
#endif
  case DIPOLAR_SCAFACOS:
    n_dipolar = 2;
    break;
  case DIPOLAR_MDLC_P3M:
  case DIPOLAR_MDLC_DS:
    n_dipolar = 3;
    break;
  default:
    break;
  }
}

int set_mesh() {
  switch (dipole.method) {
#ifdef DP3M
  case DIPOLAR_MDLC_P3M:
  case DIPOLAR_P3M:
    set_method_local(DIPOLAR_MDLC_P3M);
    return 0;
#endif
  case DIPOLAR_MDLC_DS:
  case DIPOLAR_DS:
    set_method_local(DIPOLAR_MDLC_DS);
    return 0;
  default:
    return 1;
  }
}

void bcast_params(const boost::mpi::communicator &comm) {
  namespace mpi = boost::mpi;

  switch (dipole.method) {
#ifdef DP3M
  case DIPOLAR_MDLC_P3M:
    mpi::broadcast(comm, dlc_params, 0);
    // fall through
  case DIPOLAR_P3M:
    mpi::broadcast(comm, dp3m.params, 0);
    break;
#endif
  default:
    break;
  }
}

int set_Dprefactor(double prefactor) {
  if (prefactor < 0.0) {
    runtimeErrorMsg() << "Dipolar prefactor has to be >= 0";
    return ES_ERROR;
  }

  dipole.prefactor = prefactor;

  mpi_bcast_coulomb_params();
  return ES_OK;
}

void set_method_local(DipolarInteraction method) {
#ifdef DIPOLAR_DIRECT_SUM
  if ((dipole.method == DIPOLAR_DS_GPU) && (method != DIPOLAR_DS_GPU)) {
    deactivate_dipolar_direct_sum_gpu();
  }
#endif
#ifdef DIPOLAR_BARNES_HUT
  if ((dipole.method == DIPOLAR_BH_GPU) && (method != DIPOLAR_BH_GPU)) {
    deactivate_dipolar_barnes_hut();
  }
#endif // BARNES_HUT
  dipole.method = method;
}

} // namespace Dipole
#endif // DIPOLES
