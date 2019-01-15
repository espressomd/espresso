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
 *  Routines to calculate the Gay-Berne energy and force
 *  for a pair of particles.
 *  \ref forces.cpp
 */

#include "nonbonded_interaction_data.hpp"
#include "particle_data.hpp"
#include "utils.hpp"
#include "utils/math/int_pow.hpp"
#include "utils/math/sqr.hpp"

#ifdef GAY_BERNE

///
int gay_berne_set_params(int part_type_a, int part_type_b, double eps,
                         double sig, double cut, double k1, double k2,
                         double mu, double nu);

inline void
add_gb_pair_force(const Particle *const p1, const Particle *const p2,
                  IA_parameters *ia_params, double d[3], double dist,
                  double force[3], double torque1[3], double torque2[3])

{
  if (!(dist < ia_params->GB_cut))
    return;

  double a, b, c, X, Xcut, Brack, BrackCut, Bra12, Bra12Cut, u1x, u1y, u1z, u2x,
      u2y, u2z, E, E1, E2, Sigma, Brhi1, Brhi2, Plus1, Minus1, Plus2, Minus2,
      Koef1, Koef2,               /*  mu/E2  and  Sigma^3/2  */
      dU_dr, dU_da, dU_db, dU_dc, /*  all derivatives        */
      FikX, FikY, FikZ,           /*  help for forces        */
      Gx, Gy, Gz;                 /*  help for torques       */

  u1x = p1->r.calc_director()[0];
  u1y = p1->r.calc_director()[1];
  u1z = p1->r.calc_director()[2];
  u2x = p2->r.calc_director()[0];
  u2y = p2->r.calc_director()[1];
  u2z = p2->r.calc_director()[2];

  a = d[0] * u1x + d[1] * u1y + d[2] * u1z;
  b = d[0] * u2x + d[1] * u2y + d[2] * u2z;
  c = u1x * u2x + u1y * u2y + u1z * u2z;
  E1 = 1 / sqrt(1 - ia_params->GB_chi1 * ia_params->GB_chi1 * c * c);
  Plus1 = (a + b) / (1 + ia_params->GB_chi1 * c);
  Plus2 = (a + b) / (1 + ia_params->GB_chi2 * c);
  Minus1 = (a - b) / (1 - ia_params->GB_chi1 * c);
  Minus2 = (a - b) / (1 - ia_params->GB_chi2 * c);
  Brhi2 =
      (ia_params->GB_chi2 / dist / dist) * (Plus2 * (a + b) + Minus2 * (a - b));
  E2 = 1 - 0.5 * Brhi2;
  E = 4 * ia_params->GB_eps * pow(E1, ia_params->GB_nu) *
      pow(E2, ia_params->GB_mu);
  Brhi1 =
      (ia_params->GB_chi1 / dist / dist) * (Plus1 * (a + b) + Minus1 * (a - b));
  Sigma = ia_params->GB_sig / sqrt(1 - 0.5 * Brhi1);
  Koef1 = ia_params->GB_mu / E2;
  Koef2 = Sigma * Sigma * Sigma * 0.5;

  X = ia_params->GB_sig / (dist - Sigma + ia_params->GB_sig);
  Xcut = ia_params->GB_sig / (ia_params->GB_cut - Sigma + ia_params->GB_sig);

  if (X < 1.25) { /* 1.25 corresponds to the interparticle penetration of 0.2
                    units of length.
                    If they are not that close, the GB forces and torques are
                    calculated */

    Brack = X * X * X;
    BrackCut = Xcut * Xcut * Xcut;
    Brack = Brack * Brack;
    BrackCut = BrackCut * BrackCut;

    Bra12 = 6 * Brack * X * (2 * Brack - 1);
    Bra12Cut = 6 * BrackCut * Xcut * (2 * BrackCut - 1);
    Brack = Brack * (Brack - 1);
    BrackCut = BrackCut * (BrackCut - 1);

    /*-------- Here we calculate derivatives -----------------------------*/

    dU_dr = E *
            (Koef1 * Brhi2 * (Brack - BrackCut) -
             Koef2 * Brhi1 * (Bra12 - Bra12Cut) - Bra12 * dist) /
            dist / dist;
    Koef1 = Koef1 * ia_params->GB_chi2 / dist / dist;
    Koef2 = Koef2 * ia_params->GB_chi1 / dist / dist;
    dU_da = E * (Koef1 * (Minus2 + Plus2) * (BrackCut - Brack) +
                 Koef2 * (Plus1 + Minus1) * (Bra12 - Bra12Cut));
    dU_db = E * (Koef1 * (Minus2 - Plus2) * (Brack - BrackCut) +
                 Koef2 * (Plus1 - Minus1) * (Bra12 - Bra12Cut));
    dU_dc = E * ((Brack - BrackCut) *
                     (ia_params->GB_nu * E1 * E1 * ia_params->GB_chi1 *
                          ia_params->GB_chi1 * c +
                      0.5 * Koef1 * ia_params->GB_chi2 *
                          (Plus2 * Plus2 - Minus2 * Minus2)) -
                 (Bra12 - Bra12Cut) * 0.5 * Koef2 * ia_params->GB_chi1 *
                     (Plus1 * Plus1 - Minus1 * Minus1));

    /*--------------------------------------------------------------------*/

    FikX = -dU_dr * d[0] - dU_da * u1x - dU_db * u2x;
    FikY = -dU_dr * d[1] - dU_da * u1y - dU_db * u2y;
    FikZ = -dU_dr * d[2] - dU_da * u1z - dU_db * u2z;

    force[0] += FikX;
    force[1] += FikY;
    force[2] += FikZ;

    if (torque1 != nullptr) {
      /* calculate torque:  torque = u_1 x G   */

      Gx = -dU_da * d[0] - dU_dc * u2x;
      Gy = -dU_da * d[1] - dU_dc * u2y;
      Gz = -dU_da * d[2] - dU_dc * u2z;

      torque1[0] += u1y * Gz - u1z * Gy;
      torque1[1] += u1z * Gx - u1x * Gz;
      torque1[2] += u1x * Gy - u1y * Gx;

      if (torque2 != nullptr) {
        /* calculate torque:  torque = u_2 x G     */

        Gx = -dU_db * d[0] - dU_dc * u1x;
        Gy = -dU_db * d[1] - dU_dc * u1y;
        Gz = -dU_db * d[2] - dU_dc * u1z;

        torque2[0] += u2y * Gz - u2z * Gy;
        torque2[1] += u2z * Gx - u2x * Gz;
        torque2[2] += u2x * Gy - u2y * Gx;
      }
    }
  } else { /* the particles are too close to each other */
    Koef1 = 100;

    force[0] += Koef1 * d[0];
    force[1] += Koef1 * d[1];
    force[2] += Koef1 * d[2];
  }
}

inline double gb_pair_energy(const Particle *p1, const Particle *p2,
                             const IA_parameters *ia_params, const double d[3],
                             double dist) {
  using Utils::int_pow;
  using Utils::sqr;

  if (!(dist < ia_params->GB_cut))
    return 0.0;

  auto const e0 = ia_params->GB_eps;
  auto const s0 = ia_params->GB_sig;
  auto const chi1 = ia_params->GB_chi1;
  auto const chi2 = ia_params->GB_chi2;
  auto const mu = ia_params->GB_mu;
  auto const nu = ia_params->GB_nu;
  auto const r = Vector3d({d[0], d[1], d[2]}).normalize();

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

  return E(r_eff(dist)) - E(r_eff(ia_params->GB_cut));
}

#endif
#endif
