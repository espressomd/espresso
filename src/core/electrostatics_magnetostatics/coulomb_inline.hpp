#ifndef ESPRESSO_COULOMB_INLINE_HPP
#define ESPRESSO_COULOMB_INLINE_HPP

#include "electrostatics_magnetostatics/coulomb.hpp"

#ifdef ELECTROSTATICS
#include <utils/Vector.hpp>
#include <utils/math/tensor_product.hpp>

#include "electrostatics_magnetostatics/debye_hueckel.hpp"
#include "electrostatics_magnetostatics/elc.hpp"
#include "electrostatics_magnetostatics/mmm1d.hpp"
#include "electrostatics_magnetostatics/mmm2d.hpp"
#include "electrostatics_magnetostatics/p3m.hpp"
#include "electrostatics_magnetostatics/reaction_field.hpp"
#include "electrostatics_magnetostatics/scafacos.hpp"

namespace Coulomb {
inline Utils::Vector3d central_force(double const q1q2, const double *d,
                                     double dist) {
  Utils::Vector3d f{};

  switch (coulomb.method) {
#ifdef P3M
  case COULOMB_P3M_GPU:
  case COULOMB_P3M:
  case COULOMB_ELC_P3M:
    p3m_add_pair_force(q1q2, d, dist, f.data());
    break;
#endif
  case COULOMB_MMM1D:
    add_mmm1d_coulomb_pair_force(q1q2, d, dist, f.data());
    break;
  case COULOMB_MMM2D:
    add_mmm2d_coulomb_pair_force(q1q2, d, dist, f.data());
    break;
  case COULOMB_DH:
    add_dh_coulomb_pair_force(q1q2, d, dist, f.data());
    break;
  case COULOMB_RF:
    add_rf_coulomb_pair_force(q1q2, d, dist, f.data());
    break;
#ifdef SCAFACOS
  case COULOMB_SCAFACOS:
    Scafacos::add_pair_force(q1q2, d, dist, f.data());
    break;
#endif
  default:
    break;
  }

  return coulomb.prefactor * f;
}

inline void calc_pair_force(Particle *p1, Particle *p2, const double *d,
                            double dist, Utils::Vector3d &force) {
  auto const q1q2 = p1->p.q * p2->p.q;

  if (q1q2 == 0)
    return;

  force += central_force(q1q2, d, dist);

#ifdef P3M
  if ((coulomb.method == COULOMB_ELC_P3M) &&
      (elc_params.dielectric_contrast_on)) {
    // forces from the virtual charges
    // they go directly onto the particles, since they are not pairwise forces
    Utils::Vector3d f1{};
    Utils::Vector3d f2{};

    ELC_P3M_dielectric_layers_force_contribution(p1, p2, f1.data(), f2.data());

    p1->f.f += coulomb.prefactor * f1;
    p2->f.f += coulomb.prefactor * f2;
  }
#endif
}

/**
 * @brief Pair contribution to the pressure tensor.
 *
 * If supported by the method, this returns the virial
 * contribution to the pressure tensor for this pair.
 *
 * @param p1 Particle
 * @param p2 Particle
 * @param d  Distance
 * @param dist |d|
 * @return Contribution to the pressure tensor.
 */
inline Utils::Vector<Utils::Vector3d, 3> pair_pressure(const Particle *p1,
                                                       const Particle *p2,
                                                       const Utils::Vector3d &d,
                                                       double dist) {
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
    auto const force = central_force(p1->p.q * p2->p.q, d.data(), dist);

    return Utils::tensor_product(force, d);
  }
  default:
    fprintf(stderr, "calculating pressure for electrostatics method that "
                    "doesn't have it implemented\n");
    break;
  }

  return {};
}

// energy_inline
inline double pair_energy(const Particle *p1, const Particle *p2,
                          double const q1q2, const double *d, double dist,
                          double dist2) {
  /* real space Coulomb */
  auto E = [&]() {
    switch (coulomb.method) {
#ifdef P3M
    case COULOMB_P3M_GPU:
    case COULOMB_P3M:
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
      return Scafacos::pair_energy(q1q2, dist);
#endif
    case COULOMB_DH:
      return dh_coulomb_pair_energy(q1q2, dist);
    case COULOMB_RF:
      return rf_coulomb_pair_energy(q1q2, dist);
    case COULOMB_MMM1D:
      return mmm1d_coulomb_pair_energy(q1q2, d, dist2, dist);
    case COULOMB_MMM2D:
      return mmm2d_coulomb_pair_energy(q1q2, d, dist);
    default:
      return 0.;
    }
  }();
  return coulomb.prefactor * E;
}
} // namespace Coulomb

#endif
#endif // ESPRESSO_COULOMB_INLINE_HPP
