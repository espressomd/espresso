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
#ifndef ESPRESSO_COULOMB_INLINE_HPP
#define ESPRESSO_COULOMB_INLINE_HPP

#include "electrostatics_magnetostatics/coulomb.hpp"

#ifdef ELECTROSTATICS
#include <utils/Vector.hpp>
#include <utils/math/tensor_product.hpp>

#include "electrostatics_magnetostatics/debye_hueckel.hpp"
#include "electrostatics_magnetostatics/elc.hpp"
#include "electrostatics_magnetostatics/mmm1d.hpp"
#include "electrostatics_magnetostatics/p3m.hpp"
#include "electrostatics_magnetostatics/reaction_field.hpp"
#include "electrostatics_magnetostatics/scafacos.hpp"

namespace Coulomb {
inline Utils::Vector3d central_force(double const q1q2,
                                     Utils::Vector3d const &d, double dist) {
  Utils::Vector3d f{};

  switch (coulomb.method) {
#ifdef P3M
  case COULOMB_P3M_GPU:
  case COULOMB_P3M:
  case COULOMB_ELC_P3M:
    p3m_add_pair_force(q1q2, d, dist, f);
    break;
#endif
  case COULOMB_MMM1D:
    add_mmm1d_coulomb_pair_force(q1q2, d, dist, f);
    break;
  case COULOMB_DH:
    add_dh_coulomb_pair_force(q1q2, d, dist, f);
    break;
  case COULOMB_RF:
    add_rf_coulomb_pair_force(q1q2, d, dist, f);
    break;
#ifdef SCAFACOS
  case COULOMB_SCAFACOS:
    Scafacos::add_pair_force(q1q2, d.data(), dist, f.data());
    break;
#endif
  default:
    break;
  }

  return coulomb.prefactor * f;
}

inline std::tuple<Utils::Vector3d, Utils::Vector3d, Utils::Vector3d>
pair_force(Particle const &p1, Particle const &p2, Utils::Vector3d const &d,
           double dist) {
  auto const q1q2 = p1.p.q * p2.p.q;

  if (q1q2 == 0) {
    return std::make_tuple(Utils::Vector3d{}, Utils::Vector3d{},
                           Utils::Vector3d{});
  }

  auto const force = central_force(q1q2, d, dist);
  Utils::Vector3d f1{};
  Utils::Vector3d f2{};

#ifdef P3M
  if ((coulomb.method == COULOMB_ELC_P3M) &&
      (elc_params.dielectric_contrast_on)) {
    // forces from the virtual charges
    // they go directly onto the particles, since they are not pairwise forces

    ELC_P3M_dielectric_layers_force_contribution(p1, p2, f1, f2);

    f1 *= coulomb.prefactor;
    f2 *= coulomb.prefactor;
  }
#endif

  return std::make_tuple(force, f1, f2);
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
inline Utils::Vector<Utils::Vector3d, 3> pair_pressure(Particle const &p1,
                                                       Particle const &p2,
                                                       Utils::Vector3d const &d,
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
    auto const force = central_force(p1.p.q * p2.p.q, d, dist);

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
inline double pair_energy(Particle const &p1, Particle const &p2,
                          double const q1q2, Utils::Vector3d const &d,
                          double dist, double dist2) {
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
    default:
      return 0.;
    }
  }();
  return coulomb.prefactor * E;
}
} // namespace Coulomb

#endif
#endif // ESPRESSO_COULOMB_INLINE_HPP
