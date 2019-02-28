#ifndef ESPRESSO_COULOMB_SWITCH_HPP
#define ESPRESSO_COULOMB_SWITCH_HPP

#include "electrostatics_magnetostatics/debye_hueckel.hpp" // Debye_hueckel_params
#include "electrostatics_magnetostatics/elc.hpp"           // elc_params
#include "electrostatics_magnetostatics/mmm1d.hpp"         // MMM1D_sanity_check
#include "electrostatics_magnetostatics/mmm2d.hpp"         // MMM2D_sanity_check
#include "electrostatics_magnetostatics/p3m.hpp"           // p3m_params
#include "electrostatics_magnetostatics/reaction_field.hpp" // Reaction_field_params
#include "electrostatics_magnetostatics/scafacos.hpp"       // scafacos
#include "integrate.hpp"                                    // skin
#include "npt.hpp"                                          // nptiso
#include "statistics.hpp"                                   // Observable_stat

#ifdef ELECTROSTATICS

namespace Coulomb {

// pressure.cpp
void pressure_n(int &n_coulomb);
void calc_pressure_long_range(Observable_stat &virials,
                              Observable_stat &p_tensor);

// nonbonded_interaction_data
void sanity_checks(int &state);
void cutoff(double &ret);
void deactivate();

// integrate
void integrate_sanity_check();

// initialize
void on_observable_calc();
void on_coulomb_change();
void on_resort_particles();
void on_boxl_change();
void init();

// forces_calc
void calc_long_range_force();

// energy
void calc_energy_long_range(Observable_stat &energy);
void energy_n(int &n_coulomb);

// icc
int iccp3m_sanity_check();

// elc
int elc_sanity_check();

// communication
void bcast_coulomb_params();

// forces_inline
    inline void calc_pair_force(Particle *p1, Particle *p2, double *d, double dist,
                                double dist2, Vector3d &force) {
        auto const q1q2 = p1->p.q * p2->p.q;

        if (q1q2 != 0) {
            switch (coulomb.method) {
#ifdef P3M
                case COULOMB_ELC_P3M: {
                    p3m_add_pair_force(q1q2, d, dist2, dist, force.data());

                    // forces from the virtual charges
                    // they go directly onto the particles, since they are not pairwise forces
                    if (elc_params.dielectric_contrast_on)
                        ELC_P3M_dielectric_layers_force_contribution(p1, p2, p1->f.f.data(),
                                                                     p2->f.f.data());
                    break;
                }
                case COULOMB_P3M_GPU:
                case COULOMB_P3M: {
#ifdef NPT
                    double eng = p3m_add_pair_force(q1q2, d, dist2, dist, force.data());
                    if (integ_switch == INTEG_METHOD_NPT_ISO)
                        nptiso.p_vir[0] += eng;
#else
                    p3m_add_pair_force(q1q2, d, dist2, dist, force.data());
#endif
                    break;
                }
#endif
                case COULOMB_MMM1D:
                    add_mmm1d_coulomb_pair_force(q1q2, d, dist2, dist, force.data());
                    break;
                case COULOMB_MMM2D:
                    add_mmm2d_coulomb_pair_force(q1q2, d, dist2, dist, force.data());
                    break;
                case COULOMB_DH:
                    add_dh_coulomb_pair_force(p1, p2, d, dist, force.data());
                    break;
                case COULOMB_RF:
                    add_rf_coulomb_pair_force(p1, p2, d, dist, force.data());
                    break;
#ifdef SCAFACOS
                case COULOMB_SCAFACOS:
                    Scafacos::add_pair_force(p1, p2, d, dist, force.data());
                    break;
#endif
                default:
                    break;
            }
        }
    }

// pressure_inline.hpp
inline void add_pair_pressure(Particle *p1, Particle *p2, double *d,
                              double dist, double dist2,
                              Observable_stat &virials,
                              Observable_stat &p_tensor) {
  switch (coulomb.method) {
#ifdef P3M
  case COULOMB_P3M_GPU:
  case COULOMB_P3M:
#endif
  case COULOMB_MMM1D:
  case COULOMB_DH:
  case COULOMB_RF: {
    Vector3d force{};
    calc_pair_force(p1, p2, d, dist, dist2, force);

    /* Calculate the virial pressure */
    for (int k = 0; k < 3; k++) {
        for (int l = 0; l < 3; l++) {
            p_tensor.coulomb[k * 3 + l] += force[k] * d[l];
        }
        virials.coulomb[0] += d[k] * force[k];
    }
    break;
  }
  default:
    fprintf(stderr, "calculating pressure for electrostatics method that "
                    "doesn't have it implemented\n");
    break;
  }
}

// energy_inline
inline void add_pair_energy(Particle *p1, Particle *p2, double *d, double dist,
                            double dist2, Observable_stat &energy) {
  double ret = 0;
  if (coulomb.method != COULOMB_NONE) {
    /* real space Coulomb */
    switch (coulomb.method) {
#ifdef P3M
    case COULOMB_P3M_GPU:
    case COULOMB_P3M:
      ret = p3m_pair_energy(p1->p.q * p2->p.q, dist);
      break;
    case COULOMB_ELC_P3M:
      ret = p3m_pair_energy(p1->p.q * p2->p.q, dist);
      if (elc_params.dielectric_contrast_on)
        ret += 0.5 * ELC_P3M_dielectric_layers_energy_contribution(p1, p2);
      break;
#endif
#ifdef SCAFACOS
    case COULOMB_SCAFACOS:
      ret += Scafacos::pair_energy(p1, p2, dist);
      break;
#endif
    case COULOMB_DH:
      ret = dh_coulomb_pair_energy(p1, p2, dist);
      break;
    case COULOMB_RF:
      ret = rf_coulomb_pair_energy(p1, p2, dist);
      break;
    case COULOMB_MMM1D:
      ret = mmm1d_coulomb_pair_energy(p1, p2, d, dist2, dist);
      break;
    case COULOMB_MMM2D:
      ret = mmm2d_coulomb_pair_energy(p1->p.q * p2->p.q, d, dist2, dist);
      break;
    default:
      ret = 0.;
    }
    energy.coulomb[0] += ret;
  }
}

} // namespace Coulomb
#endif // ELECTROSTATICS
#endif // ESPRESSO_COULOMB_SWITCH_HPP