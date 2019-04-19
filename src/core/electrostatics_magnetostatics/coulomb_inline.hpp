#ifndef ESPRESSO_COULOMB_INLINE_HPP
#define ESPRESSO_COULOMB_INLINE_HPP

#include "electrostatics_magnetostatics/coulomb.hpp"

#ifdef ELECTROSTATICS

#include "electrostatics_magnetostatics/debye_hueckel.hpp"
#include "electrostatics_magnetostatics/elc.hpp"
#include "electrostatics_magnetostatics/mmm1d.hpp"
#include "electrostatics_magnetostatics/mmm2d.hpp"
#include "electrostatics_magnetostatics/p3m.hpp"
#include "electrostatics_magnetostatics/reaction_field.hpp"
#include "electrostatics_magnetostatics/scafacos.hpp"

namespace Coulomb {
// forces_inline
inline void calc_pair_force(Particle *p1, Particle *p2, double const q1q2,
                            double *d, double dist, double const dist2,
                            Utils::Vector3d &force) {

  if (q1q2 != 0) {
    Utils::Vector3d f{};

    switch (coulomb.method) {
#ifdef P3M
    case COULOMB_ELC_P3M: {
      p3m_add_pair_force(q1q2, d, dist2, dist, f.data());

      // forces from the virtual charges
      // they go directly onto the particles, since they are not pairwise forces
      if (elc_params.dielectric_contrast_on) {
        Utils::Vector3d f1{};
        Utils::Vector3d f2{};

        ELC_P3M_dielectric_layers_force_contribution(p1, p2, f1.data(),
                                                     f2.data());

        p1->f.f += coulomb.prefactor * f1;
        p2->f.f += coulomb.prefactor * f2;
      }
      break;
    }
    case COULOMB_P3M_GPU:
    case COULOMB_P3M: {
      p3m_add_pair_force(q1q2, d, dist2, dist, f.data());
      break;
    }
#endif
    case COULOMB_MMM1D:
      add_mmm1d_coulomb_pair_force(q1q2, d, dist2, dist, f.data());
      break;
    case COULOMB_MMM2D:
      add_mmm2d_coulomb_pair_force(q1q2, d, dist2, dist, f.data());
      break;
    case COULOMB_DH:
      add_dh_coulomb_pair_force(q1q2, d, dist, f.data());
      break;
    case COULOMB_RF:
      add_rf_coulomb_pair_force(q1q2, d, dist, f.data());
      break;
#ifdef SCAFACOS
    case COULOMB_SCAFACOS:
      Scafacos::add_pair_force(p1, p2, d, dist, f.data());
      break;
#endif
    default:
      break;
    }

    force += coulomb.prefactor * f;
  }
}

// pressure_inline.hpp
inline void add_pair_pressure(Particle *p1, Particle *p2, double q1q2,
                              double *d, double dist, double dist2,
                              Observable_stat &virials,
                              Observable_stat &p_tensor) {
  switch (coulomb.method) {
  case COULOMB_NONE:
    break;
#ifdef P3M
  case COULOMB_P3M_GPU:
  case COULOMB_P3M:
#endif
  case COULOMB_MMM1D:
  case COULOMB_DH:
  case COULOMB_RF: {
    Utils::Vector3d force{};
    calc_pair_force(p1, p2, q1q2, d, dist, dist2, force);

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
inline double add_pair_energy(Particle *p1, Particle *p2, double const q1q2,
                              double *d, double dist, double dist2) {
  /* real space Coulomb */
  auto E = [&]() {
    switch (coulomb.method) {
#ifdef P3M
    case COULOMB_P3M_GPU:
    case COULOMB_P3M:
      // TODO some energy functions include the prefactor, some don't
      return p3m_pair_energy(q1q2, dist);
    case COULOMB_ELC_P3M:
      if (elc_params.dielectric_contrast_on) {
        return 0.5 * ELC_P3M_dielectric_layers_energy_contribution(p1, p2) +
               p3m_pair_energy(q1q2, dist);
      } else {
        return p3m_pair_energy(q1q2, dist);
      }
#endif
#ifdef SCAFACOS
    case COULOMB_SCAFACOS:
      return Scafacos::pair_energy(p1, p2, dist);
#endif
    case COULOMB_DH:
      return dh_coulomb_pair_energy(q1q2, dist);
    case COULOMB_RF:
      return rf_coulomb_pair_energy(q1q2, dist);
    case COULOMB_MMM1D:
      return mmm1d_coulomb_pair_energy(q1q2, d, dist2, dist);
    case COULOMB_MMM2D:
      return mmm2d_coulomb_pair_energy(q1q2, d, dist2, dist);
    default:
      return 0.;
    }
  }();
  return coulomb.prefactor * E;
}
} // namespace Coulomb

#endif
#endif // ESPRESSO_COULOMB_INLINE_HPP
