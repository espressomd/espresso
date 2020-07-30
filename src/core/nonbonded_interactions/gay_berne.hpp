/*
 * Copyright (C) 2010-2019 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
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
#ifndef _GB_HPP
#define _GB_HPP

/** \file
 *  Routines to calculate the Gay-Berne potential between particle pairs.
 *
 *  Please note that contrary to the original Gay and Berne 1981 paper, here
 *  the GB prefactor @f$ \varepsilon_0 @f$ is not raised to the power of
 *  @f$ \nu @f$, in agreement with the convention used in the GB literature.
 *
 *  Implementation in \ref gay_berne.cpp.
 */

#include "nonbonded_interaction_data.hpp"
#include <utils/math/int_pow.hpp>
#include <utils/math/sqr.hpp>

#ifdef GAY_BERNE

int gay_berne_set_params(int part_type_a, int part_type_b, double eps,
                         double sig, double cut, double k1, double k2,
                         double mu, double nu);

/** Calculate Gay-Berne force and torques */
inline std::tuple<Utils::Vector3d, Utils::Vector3d, Utils::Vector3d>
gb_pair_force(Utils::Vector3d const &ui, Utils::Vector3d const &uj,
              IA_parameters const &ia_params, Utils::Vector3d const &d,
              double dist, bool calc_torque1, bool calc_torque2) {
  using Utils::int_pow;
  using Utils::sqr;

  if (dist >= ia_params.gay_berne.cut) {
    return std::make_tuple(Utils::Vector3d{}, Utils::Vector3d{},
                           Utils::Vector3d{});
  }

  auto const e0 = ia_params.gay_berne.eps;
  auto const s0 = ia_params.gay_berne.sig;
  auto const chi1 = ia_params.gay_berne.chi1;
  auto const chi2 = ia_params.gay_berne.chi2;
  auto const mu = ia_params.gay_berne.mu;
  auto const nu = ia_params.gay_berne.nu;
  auto const r = Utils::Vector3d(d).normalize();
  auto const dui = d * ui;
  auto const duj = d * uj;
  auto const rui = r * ui;
  auto const ruj = r * uj;
  auto const uij = ui * uj;
  auto const oo1 = (dui + duj) / (1 + chi1 * uij);
  auto const oo2 = (dui - duj) / (1 - chi1 * uij);
  auto const tt1 = (dui + duj) / (1 + chi2 * uij);
  auto const tt2 = (dui - duj) / (1 - chi2 * uij);
  auto const o1 = sqr(rui + ruj) / (1 + chi1 * uij);
  auto const o2 = sqr(rui - ruj) / (1 - chi1 * uij);
  auto const t1 = sqr(rui + ruj) / (1 + chi2 * uij);
  auto const t2 = sqr(rui - ruj) / (1 - chi2 * uij);
  auto const Brhi1 = chi1 * (o1 + o2);
  auto const Brhi2 = chi2 * (t1 + t2);

  auto const e1 = 1. / (1. - sqr(chi1 * uij));
  auto const e2 = 1. - 0.5 * Brhi2;
  auto const e = 4 * e0 * pow(e1, 0.5 * nu) * pow(e2, mu);

  auto const s1 = 1. / sqrt(1. - 0.5 * Brhi1);
  auto const s = s0 * s1;
  auto Koef1 = mu / e2;
  auto Koef2 = int_pow<3>(s1) * 0.5;

  Utils::Vector3d force, torque1, torque2;

  auto const X = s0 / (dist - s + s0);
  auto const Xcut = s0 / (ia_params.gay_berne.cut - s + s0);

  auto const X6 = int_pow<6>(X);
  auto const Xcut6 = int_pow<6>(Xcut);

  auto const Bra12 = 6 * X6 * X * (2 * X6 - 1);
  auto const Bra12Cut = 6 * Xcut6 * Xcut * (2 * Xcut6 - 1);
  auto const Brack = X6 * (X6 - 1);
  auto const BrackCut = Xcut6 * (Xcut6 - 1);

  /*-------- Here we calculate derivatives -----------------------------*/

  auto const dU_dr = e *
                     (Koef1 * Brhi2 * (Brack - BrackCut) -
                      Koef2 * Brhi1 * (Bra12 - Bra12Cut) - Bra12 * dist / s0) /
                     sqr(dist);

  Koef1 *= chi2 / sqr(dist);
  Koef2 *= chi1 / sqr(dist);

  auto const dU_da = e * (Koef1 * (tt1 + tt2) * (BrackCut - Brack) +
                          Koef2 * (oo1 + oo2) * (Bra12 - Bra12Cut));
  auto const dU_db = e * (Koef1 * (tt2 - tt1) * (Brack - BrackCut) +
                          Koef2 * (oo1 - oo2) * (Bra12 - Bra12Cut));
  auto const dU_dc =
      e * ((Brack - BrackCut) * (nu * e1 * sqr(chi1) * uij +
                                 0.5 * Koef1 * chi2 * (sqr(tt1) - sqr(tt2))) -
           (Bra12 - Bra12Cut) * 0.5 * Koef2 * chi1 * (sqr(oo1) - sqr(oo2)));

  /*--------------------------------------------------------------------*/

  force = -dU_dr * d - dU_da * ui - dU_db * uj;

  if (calc_torque1) {
    /* calculate torque:  torque = u_i x G   */
    auto const G2 = -dU_da * d - dU_dc * uj;
    torque1 = vector_product(ui, G2);

    if (calc_torque2) {
      /* calculate torque:  torque = u_j x G     */
      auto const G1 = -dU_db * d - dU_dc * ui;
      torque2 = vector_product(uj, G1);
    }
  }
  return std::make_tuple(force, torque1, torque2);
}

/** Calculate Gay-Berne energy */
inline double gb_pair_energy(Utils::Vector3d const &ui,
                             Utils::Vector3d const &uj,
                             IA_parameters const &ia_params,
                             Utils::Vector3d const &d, double dist) {
  using Utils::int_pow;
  using Utils::sqr;

  if (dist >= ia_params.gay_berne.cut) {
    return 0.0;
  }

  auto const e0 = ia_params.gay_berne.eps;
  auto const s0 = ia_params.gay_berne.sig;
  auto const chi1 = ia_params.gay_berne.chi1;
  auto const chi2 = ia_params.gay_berne.chi2;
  auto const mu = ia_params.gay_berne.mu;
  auto const nu = ia_params.gay_berne.nu;
  auto const r = Utils::Vector3d(d).normalize();

  auto const uij = ui * uj;
  auto const rui = r * ui;
  auto const ruj = r * uj;

  auto const o1 = sqr(rui + ruj) / (1 + chi1 * uij);
  auto const o2 = sqr(rui - ruj) / (1 - chi1 * uij);
  auto const t1 = sqr(rui + ruj) / (1 + chi2 * uij);
  auto const t2 = sqr(rui - ruj) / (1 - chi2 * uij);

  auto const e1 = std::pow(1. - sqr(chi1 * uij), -0.5 * nu);
  auto const e2 = std::pow(1. - 0.5 * chi2 * (t1 + t2), mu);
  auto const e = e0 * e1 * e2;

  auto const s1 = 1. / std::sqrt(1. - 0.5 * chi1 * (o1 + o2));
  auto const s = s0 * s1;

  auto r_eff = [=](double r) { return (r - s + s0) / s0; };
  auto E = [=](double r) {
    return 4. * e * (int_pow<12>(1. / r) - int_pow<6>(1. / r));
  };

  return E(r_eff(dist)) - E(r_eff(ia_params.gay_berne.cut));
}

#endif
#endif
