
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

#include <boost/mpi/collectives.hpp>

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
    return 0;
  }
}

void integrate_sanity_check() {
  switch (dipole.method) {
  case DIPOLAR_NONE:
    break;
#ifdef DP3M
  case DIPOLAR_P3M:
    break;
#endif /* DP3M */
  default: {
    runtimeErrorMsg()
        << "NpT does not work with your dipolar method, please use P3M.";
  }
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

void calc_long_range_force() {
  switch (dipole.method) {
#ifdef DP3M
  case DIPOLAR_MDLC_P3M:
    add_mdlc_force_corrections();
    // fall through
  case DIPOLAR_P3M:
    dp3m_dipole_assign();
#ifdef NPT
    if (integ_switch == INTEG_METHOD_NPT_ISO) {
      nptiso.p_vir[0] += dp3m_calc_kspace_forces(1, 1);
      fprintf(stderr, "dipolar_P3M at this moment is added to p_vir[0]\n");
    } else
#endif
      dp3m_calc_kspace_forces(1, 0);

    break;
#endif
  case DIPOLAR_ALL_WITH_ALL_AND_NO_REPLICA:
    dawaanr_calculations(1, 0);
    break;
#ifdef DP3M
  case DIPOLAR_MDLC_DS:
    add_mdlc_force_corrections();
    // fall through
#endif
  case DIPOLAR_DS:
    magnetic_dipolar_direct_sum_calculations(1, 0);
    break;
  case DIPOLAR_DS_GPU:
    // Do nothing. It's an actor
    break;
#ifdef DIPOLAR_BARNES_HUT
  case DIPOLAR_BH_GPU:
    // Do nothing, it's an actor.
    break;
#endif // BARNES_HUT
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

void calc_energy_long_range(Observable_stat &energy) {
  switch (dipole.method) {
#ifdef DP3M
  case DIPOLAR_P3M:
    dp3m_dipole_assign();
    energy.dipolar[1] = dp3m_calc_kspace_forces(0, 1);
    break;
  case DIPOLAR_MDLC_P3M:
    dp3m_dipole_assign();
    energy.dipolar[1] = dp3m_calc_kspace_forces(0, 1);
    energy.dipolar[2] = add_mdlc_energy_corrections();
    break;
#endif
  case DIPOLAR_ALL_WITH_ALL_AND_NO_REPLICA:
    energy.dipolar[1] = dawaanr_calculations(0, 1);
    break;
#ifdef DP3M
  case DIPOLAR_MDLC_DS:
    energy.dipolar[1] = magnetic_dipolar_direct_sum_calculations(0, 1);
    energy.dipolar[2] = add_mdlc_energy_corrections();
    break;
#endif
  case DIPOLAR_DS:
    energy.dipolar[1] = magnetic_dipolar_direct_sum_calculations(0, 1);
    break;
  case DIPOLAR_DS_GPU:
    break;
#ifdef DIPOLAR_BARNES_HUT
  case DIPOLAR_BH_GPU:
    // Do nothing, it's an actor.
    break;
#endif // DIPOLAR_BARNES_HUT
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
  case DIPOLAR_MDLC_P3M:
    n_dipolar = 3;
    break;
  case DIPOLAR_P3M:
    n_dipolar = 2;
    break;
  case DIPOLAR_ALL_WITH_ALL_AND_NO_REPLICA:
    n_dipolar = 2;
    break;
  case DIPOLAR_MDLC_DS:
    n_dipolar = 3;
    break;
  case DIPOLAR_DS:
    n_dipolar = 2;
    break;
  case DIPOLAR_DS_GPU:
    n_dipolar = 2;
    break;
#ifdef DIPOLAR_BARNES_HUT
  case DIPOLAR_BH_GPU:
    n_dipolar = 2;
    break;
#endif
  case DIPOLAR_SCAFACOS:
    n_dipolar = 2;
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
  case DIPOLAR_NONE:
    break;
#ifdef DP3M
  case DIPOLAR_MDLC_P3M:
    mpi::broadcast(comm, dlc_params, 0);
    // fall through
  case DIPOLAR_P3M:
    mpi::broadcast(comm, dp3m.params, 0);
    break;
#endif
  case DIPOLAR_ALL_WITH_ALL_AND_NO_REPLICA:
    break;
  case DIPOLAR_MDLC_DS:
    // fall trough
  case DIPOLAR_DS:
    break;
  case DIPOLAR_DS_GPU:
    break;
#ifdef DIPOLAR_BARNES_HUT
  case DIPOLAR_BH_GPU:
    break;
#endif
  case DIPOLAR_SCAFACOS:
    break;
  default:
    break;
  }
}

int set_Dprefactor(double prefactor) {
  if (prefactor < 0.0) {
    runtimeErrorMsg() << "Dipolar prefactor has to be >=0";
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