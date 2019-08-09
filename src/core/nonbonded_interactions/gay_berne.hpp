/*
  Copyright (C) 2010-2018 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
    Max-Planck-Institute for Polymer Research, Theory Group

  This file is part of ESPResSo.

  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef _GB_HPP
#define _GB_HPP

/** \file
 *  Routines to calculate the Gay-Berne potential between particle pairs.
 *
 *  Implementation in \ref gay_berne.cpp.
 */

#include "nonbonded_interaction_data.hpp"
#include "particle_data.hpp"
#include <utils/math/int_pow.hpp>
#include <utils/math/sqr.hpp>

#ifdef GAY_BERNE

int gay_berne_set_params(int part_type_a, int part_type_b, double eps,
                         double sig, double cut, double k1, double k2,
                         double mu, double nu);

inline void add_gb_pair_force(Particle const *const p1,
                              Particle const *const p2,
                              IA_parameters const *const ia_params,
                              Utils::Vector3d const &d, double dist,
                              Utils::Vector3d &force,
                              Utils::Vector3d *const torque1,
                              Utils::Vector3d *const torque2) {
  using Utils::int_pow;
  using Utils::sqr;

  if (dist >= ia_params->gay_berne.cut) {
    return;
  }

  auto const u1 = p1->r.calc_director();
  auto const u2 = p2->r.calc_director();
  auto const a = d * u1;
  auto const b = d * u2;
  auto const c = u1 * u2;
  auto const E1 = 1 / sqrt(1 - ia_params->gay_berne.chi1 *
                                   ia_params->gay_berne.chi1 * c * c);
  auto const Plus1 = (a + b) / (1 + ia_params->gay_berne.chi1 * c);
  auto const Plus2 = (a + b) / (1 + ia_params->gay_berne.chi2 * c);
  auto const Minus1 = (a - b) / (1 - ia_params->gay_berne.chi1 * c);
  auto const Minus2 = (a - b) / (1 - ia_params->gay_berne.chi2 * c);
  auto const Brhi2 = (ia_params->gay_berne.chi2 / dist / dist) *
                     (Plus2 * (a + b) + Minus2 * (a - b));
  auto const E2 = 1 - 0.5 * Brhi2;
  auto const E = 4 * ia_params->gay_berne.eps *
                 pow(E1, ia_params->gay_berne.nu) *
                 pow(E2, ia_params->gay_berne.mu);
  auto const Brhi1 = (ia_params->gay_berne.chi1 / dist / dist) *
                     (Plus1 * (a + b) + Minus1 * (a - b));
  auto const Sigma = ia_params->gay_berne.sig / sqrt(1 - 0.5 * Brhi1);
  auto Koef1 = ia_params->gay_berne.mu / E2;
  auto Koef2 = int_pow<3>(Sigma) * 0.5;

  auto const X =
      ia_params->gay_berne.sig / (dist - Sigma + ia_params->gay_berne.sig);
  auto const Xcut =
      ia_params->gay_berne.sig /
      (ia_params->gay_berne.cut - Sigma + ia_params->gay_berne.sig);

  if (X < 1.25) { /* 1.25 corresponds to the interparticle penetration of 0.2
                    units of length.
                    If they are not that close, the GB forces and torques are
                    calculated */

    auto const X6 = int_pow<6>(X);
    auto const Xcut6 = int_pow<6>(Xcut);

    auto const Bra12 = 6 * X6 * X * (2 * X6 - 1);
    auto const Bra12Cut = 6 * Xcut6 * Xcut * (2 * Xcut6 - 1);
    auto const Brack = X6 * (X6 - 1);
    auto const BrackCut = Xcut6 * (Xcut6 - 1);

    /*-------- Here we calculate derivatives -----------------------------*/

    auto const dU_dr = E *
                       (Koef1 * Brhi2 * (Brack - BrackCut) -
                        Koef2 * Brhi1 * (Bra12 - Bra12Cut) - Bra12 * dist) /
                       sqr(dist);
    Koef1 *= ia_params->gay_berne.chi2 / sqr(dist);
    Koef2 *= ia_params->gay_berne.chi1 / sqr(dist);
    auto const dU_da = E * (Koef1 * (Minus2 + Plus2) * (BrackCut - Brack) +
                            Koef2 * (Plus1 + Minus1) * (Bra12 - Bra12Cut));
    auto const dU_db = E * (Koef1 * (Minus2 - Plus2) * (Brack - BrackCut) +
                            Koef2 * (Plus1 - Minus1) * (Bra12 - Bra12Cut));
    auto const dU_dc =
        E * ((Brack - BrackCut) * (ia_params->gay_berne.nu *
                                       sqr(E1 * ia_params->gay_berne.chi1) * c +
                                   0.5 * Koef1 * ia_params->gay_berne.chi2 *
                                       (sqr(Plus2) - sqr(Minus2))) -
             (Bra12 - Bra12Cut) * 0.5 * Koef2 * ia_params->gay_berne.chi1 *
                 (sqr(Plus1) - sqr(Minus1)));

    /*--------------------------------------------------------------------*/

    force -= dU_dr * d + dU_da * u1 + dU_db * u2;

    if (torque1 != nullptr) {
      /* calculate torque:  torque = u_1 x G   */
      auto const G2 = -dU_da * d - dU_dc * u2;
      *torque1 += vector_product(u1, G2);

      if (torque2 != nullptr) {
        /* calculate torque:  torque = u_2 x G     */
        auto const G1 = -dU_db * d - dU_dc * u1;
        *torque2 += vector_product(u2, G1);
      }
    }
  } else { /* the particles are too close to each other */
    Koef1 = 100;
    force += Koef1 * d;
  }
}

inline double gb_pair_energy(Particle const *const p1, Particle const *const p2,
                             IA_parameters const *const ia_params,
                             Utils::Vector3d const &d, double dist) {
  using Utils::int_pow;
  using Utils::sqr;

  if (dist >= ia_params->gay_berne.cut) {
    return 0.0;
  }

  auto const e0 = ia_params->gay_berne.eps;
  auto const s0 = ia_params->gay_berne.sig;
  auto const chi1 = ia_params->gay_berne.chi1;
  auto const chi2 = ia_params->gay_berne.chi2;
  auto const mu = ia_params->gay_berne.mu;
  auto const nu = ia_params->gay_berne.nu;
  auto const r = Utils::Vector3d({d[0], d[1], d[2]}).normalize();

  auto const ui = p1->r.calc_director();
  auto const uj = p2->r.calc_director();
  auto const uij = ui * uj;
  auto const rui = r * ui;
  auto const ruj = r * uj;

  auto const e1 = std::pow(1. - sqr(chi1 * uij), -0.5 * nu);

  auto const t1 = sqr(rui + ruj) / (1 + chi2 * uij);
  auto const t2 = sqr(rui - ruj) / (1 - chi2 * uij);
  auto const e2 = std::pow(1. - 0.5 * chi2 * (t1 + t2), mu);

  auto const e = e0 * e1 * e2;

  auto const s = s0 / std::sqrt(1. - 0.5 * chi1 *
                                         (sqr(rui + ruj) / (1 + chi1 * uij) +
                                          sqr(rui - ruj) / (1 - chi1 * uij)));

  auto r_eff = [=](double r) { return (r - s + s0) / s0; };
  auto E = [=](double r) {
    return 4. * e * (int_pow<12>(1. / r) - int_pow<6>(1. / r));
  };

  return E(r_eff(dist)) - E(r_eff(ia_params->gay_berne.cut));
}

#endif
#endif
