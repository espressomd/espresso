#ifndef ESPRESSO_DIPOLE_SWITCH_HPP
#define ESPRESSO_DIPOLE_SWITCH_HPP

#include "integrate.hpp"                                     // integ_switch
#include "electrostatics_magnetostatics/p3m-dipolar.hpp"     // dp3m
#include "statistics.hpp" // Observable_stat
#include "npt.hpp"                                           // nptiso

#ifdef ELECTROSTATICS
#ifdef DIPOLES

namespace Dipole {
    // pressure
    void pressure_n_dipolar(int &n_dipolar);
    void calc_long_range_dipole_force(Observable_stat &virials,
                                      Observable_stat &p_tensor);

    // nonbonded_interaction_data
    void nonbonded_sanity_check(int &state);
    void nonbonded_calc_cutoff(double &ret);

    // integrate
    void integrate_dipole_sanity_check();

    // initialize
    void on_observable_calc();
    void on_coulomb_change();
    void on_boxl_change();
    void init_dipole();

    // forces
    void calc_long_range_dipole();

    // energy
    void calc_long_range_dipole_energy(Observable_stat &energy);
    void energy_n_dipolar(int &n_dipolar);

    // mdlc_correction
    int set_dipole_mesh();

    // communication
    void bcast_dipole_params();

    // forces_inline
    inline void calc_pair_dipole_force(Particle *p1, Particle *p2, double *d,
                                       double dist, double dist2,
                                       Vector3d &force) {
      switch (coulomb.Dmethod) {
#ifdef DP3M
        case DIPOLAR_MDLC_P3M:
          // fall trough
        case DIPOLAR_P3M: {
#ifdef NPT
          double eng = dp3m_add_pair_force(p1, p2, d, dist2, dist, force.data());
          if (integ_switch == INTEG_METHOD_NPT_ISO)
            nptiso.p_vir[0] += eng;
#else
          dp3m_add_pair_force(p1, p2, d, dist2, dist, force.data());
#endif
          break;
        }
#endif /*ifdef DP3M */
        default:
          break;
      }
    }

    // energy_inline
    inline void add_pair_dipole_energy(Particle *p1, Particle *p2, double *d,
                                       double dist, double dist2,
                                       Observable_stat &energy) {
      double ret = 0;
      if (coulomb.Dmethod != DIPOLAR_NONE) {
        // ret=0;
        switch (coulomb.Dmethod) {
#ifdef DP3M
          case DIPOLAR_MDLC_P3M:
            // fall trough
          case DIPOLAR_P3M:
            ret = dp3m_pair_energy(p1, p2, d, dist2, dist);
            break;
#endif
          default:
            ret = 0;
        }
        energy.dipolar[0] += ret;
      }
    }

} // namespace Dipole




#endif // DIPOLES
#endif // ELECTROSTATICS
#endif // ESPRESSO_DIPOLE_SWITCH_HPP
