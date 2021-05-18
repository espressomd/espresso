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
/** \file
 *  MMM1D algorithm for long range %Coulomb interaction.
 *
 *  For more information about MMM1D, see \ref mmm1d.hpp "mmm1d.hpp".
 */

#include "config.hpp"

#ifdef ELECTROSTATICS

#include "electrostatics_magnetostatics/mmm1d.hpp"

#include "electrostatics_magnetostatics/common.hpp"
#include "electrostatics_magnetostatics/coulomb.hpp"
#include "electrostatics_magnetostatics/mmm-common.hpp"
#include "electrostatics_magnetostatics/mmm-modpsi.hpp"

#include "cells.hpp"
#include "errorhandling.hpp"
#include "grid.hpp"
#include "specfunc.hpp"
#include "tuning.hpp"

#include <utils/Vector.hpp>
#include <utils/constants.hpp>
#include <utils/math/sqr.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <tuple>
#include <vector>

/** How many trial calculations in @ref mmm1d_tune */
#define TEST_INTEGRATIONS 1000

/** Largest numerically stable cutoff for Bessel function. Don't
 *  change without improving the formulas.
 */
#define MAXIMAL_B_CUT 30

/** Minimal radius for the far formula in multiples of box_l[2] */
#define MIN_RAD 0.01

/* if you define this, the Bessel functions are calculated up
 * to machine precision, otherwise 10^-14, which should be
 * definitely enough for daily life. */
#undef BESSEL_MACHINE_PREC

#ifndef BESSEL_MACHINE_PREC
#define K0 LPK0
#define K1 LPK1
#endif

/** @name inverse box dimensions and other constants */
/**@{*/
static double uz2, prefuz2, prefL3_i;
/**@}*/

MMM1D_struct mmm1d_params = {0.05, 1e-5, 0};
/** From which distance a certain Bessel cutoff is valid. Can't be part of the
    params since these get broadcasted. */
static std::vector<double> bessel_radii;

static double far_error(int P, double minrad) {
  // this uses an upper bound to all force components and the potential
  auto const rhores = 2 * Utils::pi() * box_geo.length_inv()[2] * minrad;
  auto const pref = 4 * box_geo.length_inv()[2] *
                    std::max(1.0, 2 * Utils::pi() * box_geo.length_inv()[2]);

  return pref * K1(rhores * P) * exp(rhores) / rhores * (P - 1 + 1 / rhores);
}

static double determine_minrad(double maxPWerror, int P) {
  // bisection to search for where the error is maxPWerror
  double const rgranularity = MIN_RAD * box_geo.length()[2];
  double rmin = rgranularity;
  double rmax = std::min(box_geo.length()[0], box_geo.length()[1]);
  double const errmin = far_error(P, rmin);
  double const errmax = far_error(P, rmax);
  if (errmin < maxPWerror) {
    // we can do almost all radii with this P
    return rmin;
  }
  if (errmax > maxPWerror) {
    // make sure that this switching radius cannot be reached
    return 2 * std::max(box_geo.length()[0], box_geo.length()[1]);
  }

  while (rmax - rmin > rgranularity) {
    auto const c = 0.5 * (rmin + rmax);
    auto const errc = far_error(P, c);
    if (errc > maxPWerror) {
      rmin = c;
    } else {
      rmax = c;
    }
  }
  return 0.5 * (rmin + rmax);
}

static void determine_bessel_radii(double maxPWerror, int maxP) {
  bessel_radii.resize(maxP);
  for (int P = 1; P <= maxP; ++P) {
    bessel_radii[P - 1] = determine_minrad(maxPWerror, P);
  }
}

static void prepare_polygamma_series(double maxPWerror, double maxrad2) {
  /* polygamma, determine order */
  double err;
  auto const rhomax2 = uz2 * maxrad2;
  /* rhomax2 < 1, so rhomax2m2 falls monotonously */
  int n = 1;
  auto rhomax2nm2 = 1.0;
  do {
    create_mod_psi_up_to(n + 1);

    /* |uz*z| <= 0.5 */
    err = 2 * n * fabs(mod_psi_even(n, 0.5)) * rhomax2nm2;
    rhomax2nm2 *= rhomax2;
    n++;
  } while (err > 0.1 * maxPWerror);
}

void MMM1D_set_params(double switch_rad, double maxPWerror) {
  mmm1d_params.far_switch_radius_2 =
      (switch_rad > 0) ? Utils::sqr(switch_rad) : -1;
  mmm1d_params.maxPWerror = maxPWerror;
  coulomb.method = COULOMB_MMM1D;

  mpi_bcast_coulomb_params();
}

int MMM1D_sanity_checks() {
  if (box_geo.periodic(0) || box_geo.periodic(1) || !box_geo.periodic(2)) {
    runtimeErrorMsg() << "MMM1D requires periodicity (0, 0, 1)";
    return ES_ERROR;
  }
  if (cell_structure.decomposition_type() != CELL_STRUCTURE_NSQUARE) {
    runtimeErrorMsg() << "MMM1D requires the N-square cellsystem";
    return ES_ERROR;
  }
  return ES_OK;
}

int MMM1D_init() {
  if (MMM1D_sanity_checks())
    return ES_ERROR;

  if (mmm1d_params.far_switch_radius_2 >= Utils::sqr(box_geo.length()[2]))
    mmm1d_params.far_switch_radius_2 = 0.8 * Utils::sqr(box_geo.length()[2]);

  uz2 = Utils::sqr(box_geo.length_inv()[2]);
  prefuz2 = coulomb.prefactor * uz2;
  prefL3_i = prefuz2 * box_geo.length_inv()[2];

  determine_bessel_radii(mmm1d_params.maxPWerror, MAXIMAL_B_CUT);
  prepare_polygamma_series(mmm1d_params.maxPWerror,
                           mmm1d_params.far_switch_radius_2);
  return ES_OK;
}

void add_mmm1d_coulomb_pair_force(double chpref, Utils::Vector3d const &d,
                                  double r, Utils::Vector3d &force) {
  constexpr double c_2pi = 2 * Utils::pi();
  auto const n_modPsi = static_cast<int>(modPsi.size() >> 1);
  auto const rxy2 = d[0] * d[0] + d[1] * d[1];
  auto const rxy2_d = rxy2 * uz2;
  auto const z_d = d[2] * box_geo.length_inv()[2];
  Utils::Vector3d F;

  if (rxy2 <= mmm1d_params.far_switch_radius_2) {
    /* polygamma summation */
    double sr = 0;
    double sz = mod_psi_odd(0, z_d);
    double r2nm1 = 1.0;
    for (int n = 1; n < n_modPsi; n++) {
      auto const deriv = static_cast<double>(2 * n);
      auto const mpe = mod_psi_even(n, z_d);
      auto const mpo = mod_psi_odd(n, z_d);
      auto const r2n = r2nm1 * rxy2_d;

      sz += r2n * mpo;
      sr += deriv * r2nm1 * mpe;

      if (fabs(deriv * r2nm1 * mpe) < mmm1d_params.maxPWerror)
        break;

      r2nm1 = r2n;
    }

    double Fx = prefL3_i * sr * d[0];
    double Fy = prefL3_i * sr * d[1];
    double Fz = prefuz2 * sz;

    /* real space parts */

    double pref, rt, rt2, shift_z;

    pref = 1. / (r * r * r);
    Fx += pref * d[0];
    Fy += pref * d[1];
    Fz += pref * d[2];

    shift_z = d[2] + box_geo.length()[2];
    rt2 = rxy2 + shift_z * shift_z;
    rt = sqrt(rt2);
    pref = 1. / (rt2 * rt);
    Fx += pref * d[0];
    Fy += pref * d[1];
    Fz += pref * shift_z;

    shift_z = d[2] - box_geo.length()[2];
    rt2 = rxy2 + shift_z * shift_z;
    rt = sqrt(rt2);
    pref = 1. / (rt2 * rt);
    Fx += pref * d[0];
    Fy += pref * d[1];
    Fz += pref * shift_z;

    F = {Fx, Fy, Fz};
  } else {
    /* far range formula */
    auto const rxy = sqrt(rxy2);
    auto const rxy_d = rxy * box_geo.length_inv()[2];
    double sr = 0, sz = 0;

    for (int bp = 1; bp < MAXIMAL_B_CUT; bp++) {
      if (bessel_radii[bp - 1] < rxy)
        break;

      auto const fq = c_2pi * bp;
      double k0, k1;
#ifdef BESSEL_MACHINE_PREC
      k0 = K0(fq * rxy_d);
      k1 = K1(fq * rxy_d);
#else
      std::tie(k0, k1) = LPK01(fq * rxy_d);
#endif
      sr += bp * k1 * cos(fq * z_d);
      sz += bp * k0 * sin(fq * z_d);
    }
    sr *= uz2 * 4 * c_2pi;
    sz *= uz2 * 4 * c_2pi;

    auto const pref = sr / rxy + 2 * box_geo.length_inv()[2] / rxy2;

    F = {pref * d[0], pref * d[1], sz};
  }

  force += chpref * F;
}

double mmm1d_coulomb_pair_energy(double const chpref, Utils::Vector3d const &d,
                                 double r2, double r) {
  if (chpref == 0)
    return 0;

  constexpr double c_2pi = 2 * Utils::pi();
  auto const n_modPsi = static_cast<int>(modPsi.size() >> 1);
  auto const rxy2 = d[0] * d[0] + d[1] * d[1];
  auto const rxy2_d = rxy2 * uz2;
  auto const z_d = d[2] * box_geo.length_inv()[2];
  double E;

  if (rxy2 <= mmm1d_params.far_switch_radius_2) {
    /* near range formula */
    E = -2 * Utils::gamma();

    /* polygamma summation */
    double r2n = 1.0;
    for (int n = 0; n < n_modPsi; n++) {
      auto const add = mod_psi_even(n, z_d) * r2n;
      E -= add;

      if (fabs(add) < mmm1d_params.maxPWerror)
        break;

      r2n *= rxy2_d;
    }
    E *= box_geo.length_inv()[2];

    /* real space parts */

    double rt, shift_z;

    E += 1 / r;

    shift_z = d[2] + box_geo.length()[2];
    rt = sqrt(rxy2 + shift_z * shift_z);
    E += 1 / rt;

    shift_z = d[2] - box_geo.length()[2];
    rt = sqrt(rxy2 + shift_z * shift_z);
    E += 1 / rt;
  } else {
    /* far range formula */
    auto const rxy = sqrt(rxy2);
    auto const rxy_d = rxy * box_geo.length_inv()[2];
    /* The first Bessel term will compensate a little bit the
       log term, so add them close together */
    E = -0.25 * log(rxy2_d) + 0.5 * (Utils::ln_2() - Utils::gamma());
    for (int bp = 1; bp < MAXIMAL_B_CUT; bp++) {
      if (bessel_radii[bp - 1] < rxy)
        break;

      auto const fq = c_2pi * bp;
      E += K0(fq * rxy_d) * cos(fq * z_d);
    }
    E *= 4 * box_geo.length_inv()[2];
  }

  return chpref * E;
}

int mmm1d_tune(bool verbose) {
  if (MMM1D_sanity_checks())
    return ES_ERROR;
  double min_time = std::numeric_limits<double>::infinity();
  double min_rad = -1;
  auto const maxrad = box_geo.length()[2];
  double switch_radius;

  if (mmm1d_params.far_switch_radius_2 < 0) {
    /* determine besselcutoff and optimal switching radius. Should be around
     * 0.33 */
    for (switch_radius = 0.2 * maxrad; switch_radius < 0.4 * maxrad;
         switch_radius += 0.025 * maxrad) {
      if (switch_radius <= bessel_radii[MAXIMAL_B_CUT - 1]) {
        // this switching radius is too small for our Bessel series
        continue;
      }

      mmm1d_params.far_switch_radius_2 = Utils::sqr(switch_radius);

      coulomb.method = COULOMB_MMM1D;

      /* initialize mmm1d temporary structures */
      mpi_bcast_coulomb_params();

      /* perform force calculation test */
      double int_time = time_force_calc(TEST_INTEGRATIONS);

      /* exit on errors */
      if (int_time < 0)
        return ES_ERROR;

      if (verbose) {
        std::printf("r= %f t= %f ms\n", switch_radius, int_time);
      }

      if (int_time < min_time) {
        min_time = int_time;
        min_rad = switch_radius;
      }
      /* stop if all hope is vain... */
      else if (int_time > 2 * min_time)
        break;
    }
    switch_radius = min_rad;
    mmm1d_params.far_switch_radius_2 = Utils::sqr(switch_radius);
  } else if (mmm1d_params.far_switch_radius_2 <=
             Utils::sqr(bessel_radii[MAXIMAL_B_CUT - 1])) {
    // this switching radius is too small for our Bessel series
    runtimeErrorMsg() << "could not find reasonable bessel cutoff";
    return ES_ERROR;
  }

  coulomb.method = COULOMB_MMM1D;

  mpi_bcast_coulomb_params();

  return ES_OK;
}

#endif
