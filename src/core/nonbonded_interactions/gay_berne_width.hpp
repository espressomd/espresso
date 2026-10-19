/*
 * Copyright (C) 2010-2026 The ESPResSo project
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

#pragma once

/** \file
 * Routines to calculate a modified Gay--Berne potential with an
 * independently adjustable radial width.
 *
 * The original ESPResSo Gay--Berne interaction uses the radial
 * coordinate
 *
 * \f[
 * x = \frac{r-s+\sigma}{\sigma},
 * \f]
 *
 * where \f$s\f$ is the orientation-dependent contact distance.
 *
 * This interaction replaces it with
 *
 * \f[
 * x_w =
 * 2^{1/6} + \frac{r-r_{\min}}{w},
 * \f]
 *
 * where
 *
 * \f[
 * r_{\min}
 * =
 * s + \left(2^{1/6}-1\right)\sigma.
 * \f]
 *
 * Therefore, \f$\sigma\f$ retains control of the position of the
 * potential minimum, while \f$w\f$ controls the radial width and
 * local stiffness of the well.
 *
 * Setting \f$w=\sigma\f$ recovers the radial form of the original
 * ESPResSo Gay--Berne interaction.
 */

#include "config/config.hpp"

#ifdef ESPRESSO_GAY_BERNE_WIDTH

#include "Particle.hpp"
#include "nonbonded_interaction_data.hpp"

#include <utils/Vector.hpp>
#include <utils/math/int_pow.hpp>
#include <utils/math/sqr.hpp>
#include <utils/quaternion.hpp>

#include <cmath>

/**
 * Calculate the force and torque for the modified Gay--Berne interaction.
 *
 * This is the main implementation and operates directly on the particle
 * directors. The director-based form is also required by pressure and
 * virial calculations.
 *
 * @param ui Director of particle i.
 * @param uj Director of particle j.
 * @param ia_params Non-bonded interaction parameters.
 * @param d Distance vector between the particles.
 * @param dist Magnitude of the distance vector.
 *
 * @return Force acting on the particle and its torque.
 */
inline ParticleForce gb_width_pair_force(Utils::Vector3d const &ui,
                                         Utils::Vector3d const &uj,
                                         IA_parameters const &ia_params,
                                         Utils::Vector3d const &d,
                                         double dist) {
  using Utils::int_pow;
  using Utils::sqr;

  auto const &params = ia_params.gay_berne_width;

  if (dist >= params.cut) {
    return {};
  }

  auto const e0 = params.eps;

  /*
   * sig controls the reference minimum position.
   */
  auto const s0 = params.sig;

  /*
   * wid controls the radial width and stiffness of the well.
   */
  auto const w0 = params.wid;

  auto const chi1 = params.chi1;
  auto const chi2 = params.chi2;
  auto const mu = params.mu;
  auto const nu = params.nu;

  /*
   * Unit vector along the particle-particle separation.
   */
  auto const r = Utils::Vector3d(d).normalize();

  auto const dui = d * ui;
  auto const duj = d * uj;

  auto const rui = r * ui;
  auto const ruj = r * uj;
  auto const uij = ui * uj;

  auto const oo1 = (dui + duj) / (1. + chi1 * uij);
  auto const oo2 = (dui - duj) / (1. - chi1 * uij);

  auto const tt1 = (dui + duj) / (1. + chi2 * uij);
  auto const tt2 = (dui - duj) / (1. - chi2 * uij);

  auto const o1 = sqr(rui + ruj) / (1. + chi1 * uij);
  auto const o2 = sqr(rui - ruj) / (1. - chi1 * uij);

  auto const t1 = sqr(rui + ruj) / (1. + chi2 * uij);
  auto const t2 = sqr(rui - ruj) / (1. - chi2 * uij);

  auto const Brhi1 = chi1 * (o1 + o2);
  auto const Brhi2 = chi2 * (t1 + t2);

  auto const e1 = 1. / (1. - sqr(chi1 * uij));
  auto const e2 = 1. - 0.5 * Brhi2;

  /*
   * As in ESPResSo's original Gay--Berne force implementation,
   * this energy prefactor already contains the Lennard-Jones factor 4.
   */
  auto const e = 4. * e0 * std::pow(e1, 0.5 * nu) * std::pow(e2, mu);

  /*
   * Orientation-dependent contact distance:
   *
   *     s = sig / sqrt(1 - Brhi1 / 2).
   */
  auto const s1 = 1. / std::sqrt(1. - 0.5 * Brhi1);
  auto const s = s0 * s1;

  /*
   * Dimensionless position of the Lennard-Jones minimum.
   */
  auto const alpha = std::pow(2., 1. / 6.);

  /*
   * Preserve the original Gay--Berne minimum position.
   */
  auto const r_min = s + (alpha - 1.) * s0;

  /*
   * Modified inverse radial coordinate:
   *
   *     X = 1 / x_w
   *
   * where
   *
   *     x_w = alpha + (r - r_min) / wid.
   *
   * When wid == sig:
   *
   *     x_w = (r - s + sig) / sig,
   *
   * which is the original ESPResSo Gay--Berne radial coordinate.
   */
  auto const X = 1. / (alpha + (dist - r_min) / w0);

  /*
   * Value of the modified inverse radial coordinate at the cutoff.
   * This is used to shift the energy continuously to zero at the cutoff.
   */
  auto const Xcut = 1. / (alpha + (params.cut - r_min) / w0);

  auto const X6 = int_pow<6>(X);
  auto const Xcut6 = int_pow<6>(Xcut);

  /*
   * Lennard-Jones-like radial terms:
   *
   *     Brack = X^12 - X^6
   *
   * and
   *
   *     Bra12 = -d(Brack)/dx_w
   *           = 6 X^7 (2 X^6 - 1).
   */
  auto const Bra12 = 6. * X6 * X * (2. * X6 - 1.);
  auto const Bra12Cut = 6. * Xcut6 * Xcut * (2. * Xcut6 - 1.);

  auto const Brack = X6 * (X6 - 1.);
  auto const BrackCut = Xcut6 * (Xcut6 - 1.);

  /*
   * Angular energy-prefactor derivatives are unchanged from the original
   * Gay--Berne interaction.
   */
  auto Koef1 = mu / e2;
  auto Koef2 = int_pow<3>(s1) * 0.5;

  /*
   * The contact-distance derivative terms acquire sig / wid because
   *
   *     dx_w / ds = -1 / wid
   *
   * instead of
   *
   *     dx / ds = -1 / sig.
   */
  auto const width_scale = s0 / w0;

  /* Radial derivative. */
  auto const dU_dr =
      e *
      (Koef1 * Brhi2 * (Brack - BrackCut) -
       width_scale * Koef2 * Brhi1 * (Bra12 - Bra12Cut) - Bra12 * dist / w0) /
      sqr(dist);

  Koef1 *= chi2 / sqr(dist);
  Koef2 *= chi1 / sqr(dist);

  /*
   * Derivatives with respect to the orientational scalar products.
   */
  auto const dU_da =
      e * (Koef1 * (tt1 + tt2) * (BrackCut - Brack) +
           width_scale * Koef2 * (oo1 + oo2) * (Bra12 - Bra12Cut));

  auto const dU_db =
      e * (Koef1 * (tt2 - tt1) * (Brack - BrackCut) +
           width_scale * Koef2 * (oo1 - oo2) * (Bra12 - Bra12Cut));

  auto const dU_dc =
      e * ((Brack - BrackCut) * (nu * e1 * sqr(chi1) * uij +
                                 0.5 * Koef1 * chi2 * (sqr(tt1) - sqr(tt2))) -
           width_scale * (Bra12 - Bra12Cut) * 0.5 * Koef2 * chi1 *
               (sqr(oo1) - sqr(oo2)));

  /*
   * Force and torque structure follows the original ESPResSo
   * Gay--Berne implementation.
   *
   * The torque is
   *
   *     tau_i = ui x G2.
   */
  auto const G2 = -dU_da * d - dU_dc * uj;

  ParticleForce pf{
      -dU_dr * d - dU_da * ui - dU_db * uj,
      vector_product(ui, G2),
  };

  return pf;
}

/**
 * Calculate the energy of the modified Gay--Berne interaction.
 *
 * This is the main energy implementation and operates directly on
 * particle directors.
 *
 * @param ui Director of particle i.
 * @param uj Director of particle j.
 * @param ia_params Non-bonded interaction parameters.
 * @param d Distance vector between the particles.
 * @param dist Magnitude of the distance vector.
 *
 * @return Interaction energy.
 */
inline double gb_width_pair_energy(Utils::Vector3d const &ui,
                                   Utils::Vector3d const &uj,
                                   IA_parameters const &ia_params,
                                   Utils::Vector3d const &d, double dist) {
  using Utils::int_pow;
  using Utils::sqr;

  auto const &params = ia_params.gay_berne_width;

  if (dist >= params.cut) {
    return {};
  }

  auto const e0 = params.eps;
  auto const s0 = params.sig;
  auto const w0 = params.wid;

  auto const chi1 = params.chi1;
  auto const chi2 = params.chi2;
  auto const mu = params.mu;
  auto const nu = params.nu;

  auto const r = Utils::Vector3d(d).normalize();

  auto const uij = ui * uj;
  auto const rui = r * ui;
  auto const ruj = r * uj;

  auto const o1 = sqr(rui + ruj) / (1. + chi1 * uij);
  auto const o2 = sqr(rui - ruj) / (1. - chi1 * uij);

  auto const t1 = sqr(rui + ruj) / (1. + chi2 * uij);
  auto const t2 = sqr(rui - ruj) / (1. - chi2 * uij);

  /*
   * Orientation-dependent energy prefactor.
   */
  auto const e1 = std::pow(1. - sqr(chi1 * uij), -0.5 * nu);

  auto const e2 = std::pow(1. - 0.5 * chi2 * (t1 + t2), mu);

  auto const e = e0 * e1 * e2;

  /*
   * Orientation-dependent contact distance.
   */
  auto const s1 = 1. / std::sqrt(1. - 0.5 * chi1 * (o1 + o2));

  auto const s = s0 * s1;

  auto const alpha = std::pow(2., 1. / 6.);

  /*
   * Original Gay--Berne minimum position.
   */
  auto const r_min = s + (alpha - 1.) * s0;

  auto r_eff = [=](double r) { return alpha + (r - r_min) / w0; };

  auto E = [=](double x) {
    return 4. * e * (int_pow<12>(1. / x) - int_pow<6>(1. / x));
  };

  return E(r_eff(dist)) - E(r_eff(params.cut));
}

/**
 * Quaternion wrapper for the modified Gay--Berne force and torque.
 */
inline ParticleForce gb_width_pair_force(Utils::Quaternion<double> const &qi,
                                         Utils::Quaternion<double> const &qj,
                                         IA_parameters const &ia_params,
                                         Utils::Vector3d const &d,
                                         double dist) {
  auto const ui = Utils::convert_quaternion_to_director(qi);
  auto const uj = Utils::convert_quaternion_to_director(qj);

  return gb_width_pair_force(ui, uj, ia_params, d, dist);
}

/**
 * Quaternion wrapper for the modified Gay--Berne energy.
 */
inline double gb_width_pair_energy(Utils::Quaternion<double> const &qi,
                                   Utils::Quaternion<double> const &qj,
                                   IA_parameters const &ia_params,
                                   Utils::Vector3d const &d, double dist) {
  auto const ui = Utils::convert_quaternion_to_director(qi);
  auto const uj = Utils::convert_quaternion_to_director(qj);

  return gb_width_pair_energy(ui, uj, ia_params, d, dist);
}

#endif // ESPRESSO_GAY_BERNE_WIDTH
