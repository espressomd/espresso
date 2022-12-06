/*
 * Copyright (C) 2010-2022 The ESPResSo project
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

#include "config/config.hpp"

#ifdef ELECTROSTATICS

#include "electrostatics/mmm1d.hpp"

#include "electrostatics/coulomb.hpp"
#include "electrostatics/mmm-common.hpp"
#include "electrostatics/mmm-modpsi.hpp"

#include "Particle.hpp"
#include "cell_system/CellStructureType.hpp"
#include "errorhandling.hpp"
#include "event.hpp"
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
#include <vector>

/* if you define this feature, the Bessel functions are calculated up
 * to machine precision, otherwise 10^-14, which should be
 * definitely enough for daily life. */
#ifndef MMM1D_MACHINE_PREC
#define K0 LPK0
#define K1 LPK1
#endif

static double far_error(int P, double minrad) {
  auto const wavenumber = 2. * Utils::pi() * box_geo.length_inv()[2];
  // this uses an upper bound to all force components and the potential
  auto const rhores = wavenumber * minrad;
  auto const pref = 4. * box_geo.length_inv()[2] * std::max(1., wavenumber);

  return pref * K1(rhores * P) * exp(rhores) / rhores * (P - 1. + 1. / rhores);
}

static double determine_minrad(double maxPWerror, int P) {
  // bisection to search for where the error is maxPWerror
  auto constexpr min_rad = 0.01;
  auto const rgranularity = min_rad * box_geo.length()[2];
  auto rmin = rgranularity;
  auto rmax = std::min(box_geo.length()[0], box_geo.length()[1]);
  auto const errmin = far_error(P, rmin);
  auto const errmax = far_error(P, rmax);
  if (errmin < maxPWerror) {
    // we can do almost all radii with this P
    return rmin;
  }
  if (errmax > maxPWerror) {
    // make sure that this switching radius cannot be reached
    return 2. * std::max(box_geo.length()[0], box_geo.length()[1]);
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

void CoulombMMM1D::determine_bessel_radii() {
  for (int i = 0; i < MAXIMAL_B_CUT; ++i) {
    bessel_radii[i] = determine_minrad(maxPWerror, i + 1);
  }
}

void CoulombMMM1D::prepare_polygamma_series() {
  /* polygamma, determine order */
  double err;
  auto const rhomax2 = uz2 * far_switch_radius_sq;
  /* rhomax2 < 1, so rhomax2m2 falls monotonously */
  int n = 1;
  auto rhomax2nm2 = 1.0;
  do {
    create_mod_psi_up_to(n + 1);

    /* |uz*z| <= 0.5 */
    err = 2. * static_cast<double>(n) * fabs(mod_psi_even(n, 0.5)) * rhomax2nm2;
    rhomax2nm2 *= rhomax2;
    n++;
  } while (err > 0.1 * maxPWerror);
}

CoulombMMM1D::CoulombMMM1D(double prefactor, double maxPWerror,
                           double switch_rad, int tune_timings,
                           bool tune_verbose)
    : maxPWerror{maxPWerror}, far_switch_radius{switch_rad},
      tune_timings{tune_timings}, tune_verbose{tune_verbose}, m_is_tuned{false},
      far_switch_radius_sq{-1.}, uz2{0.}, prefuz2{0.}, prefL3_i{0.} {
  set_prefactor(prefactor);
  if (maxPWerror <= 0.) {
    throw std::domain_error("Parameter 'maxPWerror' must be > 0");
  }
  if (far_switch_radius <= 0. and far_switch_radius != -1.) {
    throw std::domain_error("Parameter 'far_switch_radius' must be > 0");
  }
  if (far_switch_radius > 0.) {
    far_switch_radius_sq = Utils::sqr(far_switch_radius);
  }
  if (tune_timings <= 0) {
    throw std::domain_error("Parameter 'timings' must be > 0");
  }
}

void CoulombMMM1D::sanity_checks_periodicity() const {
  if (box_geo.periodic(0) || box_geo.periodic(1) || !box_geo.periodic(2)) {
    throw std::runtime_error("MMM1D requires periodicity (False, False, True)");
  }
}

void CoulombMMM1D::sanity_checks_cell_structure() const {
  if (local_geo.cell_structure_type() !=
      CellStructureType::CELL_STRUCTURE_NSQUARE) {
    throw std::runtime_error("MMM1D requires the N-square cellsystem");
  }
}

void CoulombMMM1D::recalc_boxl_parameters() {
  if (far_switch_radius_sq >= Utils::sqr(box_geo.length()[2]))
    far_switch_radius_sq = 0.8 * Utils::sqr(box_geo.length()[2]);

  uz2 = Utils::sqr(box_geo.length_inv()[2]);
  prefuz2 = prefactor * uz2;
  prefL3_i = prefuz2 * box_geo.length_inv()[2];

  determine_bessel_radii();
  prepare_polygamma_series();
}

Utils::Vector3d CoulombMMM1D::pair_force(double q1q2, Utils::Vector3d const &d,
                                         double dist) const {
  auto constexpr c_2pi = 2. * Utils::pi();
  auto const n_modPsi = static_cast<int>(modPsi.size()) >> 1;
  auto const rxy2 = d[0] * d[0] + d[1] * d[1];
  auto const rxy2_d = rxy2 * uz2;
  auto const z_d = d[2] * box_geo.length_inv()[2];
  Utils::Vector3d force;

  if (rxy2 <= far_switch_radius_sq) {
    /* polygamma summation */
    auto sr = 0.;
    auto sz = mod_psi_odd(0, z_d);
    auto r2nm1 = 1.;
    for (int n = 1; n < n_modPsi; n++) {
      auto const deriv = static_cast<double>(2 * n);
      auto const mpe = mod_psi_even(n, z_d);
      auto const mpo = mod_psi_odd(n, z_d);
      auto const r2n = r2nm1 * rxy2_d;

      sz += r2n * mpo;
      sr += deriv * r2nm1 * mpe;

      if (fabs(deriv * r2nm1 * mpe) < maxPWerror)
        break;

      r2nm1 = r2n;
    }

    double Fx = prefL3_i * sr * d[0];
    double Fy = prefL3_i * sr * d[1];
    double Fz = prefuz2 * sz;

    /* real space parts */

    double pref, rt, rt2, shift_z;

    pref = 1. / Utils::int_pow<3>(dist);
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

    force = {Fx, Fy, Fz};
  } else {
    /* far range formula */
    auto const rxy = sqrt(rxy2);
    auto const rxy_d = rxy * box_geo.length_inv()[2];
    auto sr = 0., sz = 0.;

    for (int bp = 1; bp < MAXIMAL_B_CUT; bp++) {
      if (bessel_radii[bp - 1] < rxy)
        break;

      auto const fq = c_2pi * bp;
#ifdef MMM1D_MACHINE_PREC
      auto const k0 = K0(fq * rxy_d);
      auto const k1 = K1(fq * rxy_d);
#else
      auto const [k0, k1] = LPK01(fq * rxy_d);
#endif
      sr += bp * k1 * cos(fq * z_d);
      sz += bp * k0 * sin(fq * z_d);
    }
    sr *= uz2 * 4. * c_2pi;
    sz *= uz2 * 4. * c_2pi;

    auto const pref = sr / rxy + 2. * box_geo.length_inv()[2] / rxy2;

    force = {pref * d[0], pref * d[1], sz};
  }

  return (prefactor * q1q2) * force;
}

double CoulombMMM1D::pair_energy(double const q1q2, Utils::Vector3d const &d,
                                 double const dist) const {
  if (q1q2 == 0.)
    return 0.;

  auto constexpr c_2pi = 2. * Utils::pi();
  auto const n_modPsi = static_cast<int>(modPsi.size()) >> 1;
  auto const rxy2 = d[0] * d[0] + d[1] * d[1];
  auto const rxy2_d = rxy2 * uz2;
  auto const z_d = d[2] * box_geo.length_inv()[2];
  double energy;

  if (rxy2 <= far_switch_radius_sq) {
    /* near range formula */
    energy = -2. * Utils::gamma();

    /* polygamma summation */
    double r2n = 1.0;
    for (int n = 0; n < n_modPsi; n++) {
      auto const add = mod_psi_even(n, z_d) * r2n;
      energy -= add;

      if (fabs(add) < maxPWerror)
        break;

      r2n *= rxy2_d;
    }
    energy *= box_geo.length_inv()[2];

    /* real space parts */

    double rt, shift_z;

    energy += 1. / dist;

    shift_z = d[2] + box_geo.length()[2];
    rt = sqrt(rxy2 + shift_z * shift_z);
    energy += 1. / rt;

    shift_z = d[2] - box_geo.length()[2];
    rt = sqrt(rxy2 + shift_z * shift_z);
    energy += 1. / rt;
  } else {
    /* far range formula */
    auto const rxy = sqrt(rxy2);
    auto const rxy_d = rxy * box_geo.length_inv()[2];
    /* The first Bessel term will compensate a little bit the
       log term, so add them close together */
    energy = -0.25 * log(rxy2_d) + 0.5 * (Utils::ln_2() - Utils::gamma());
    for (int bp = 1; bp < MAXIMAL_B_CUT; bp++) {
      if (bessel_radii[bp - 1] < rxy)
        break;

      auto const fq = c_2pi * bp;
      energy += K0(fq * rxy_d) * cos(fq * z_d);
    }
    energy *= 4. * box_geo.length_inv()[2];
  }

  return prefactor * q1q2 * energy;
}

void CoulombMMM1D::tune() {
  if (is_tuned()) {
    return;
  }
  recalc_boxl_parameters();

  if (far_switch_radius_sq < 0.) {
    auto const maxrad = box_geo.length()[2];
    auto min_time = std::numeric_limits<double>::infinity();
    auto min_rad = -1.;
    auto switch_radius = 0.2 * maxrad;
    /* determine optimal switching radius. Should be around 0.33 */
    while (switch_radius < 0.4 * maxrad) {
      if (switch_radius > bessel_radii.back()) {
        // this switching radius is large enough for our Bessel series
        far_switch_radius_sq = Utils::sqr(switch_radius);
        on_coulomb_change();

        /* perform force calculation test */
        auto const int_time = benchmark_integration_step(tune_timings);

        if (tune_verbose) {
          std::printf("r= %f t= %f ms\n", switch_radius, int_time);
        }

        if (int_time < min_time) {
          min_time = int_time;
          min_rad = switch_radius;
        } else if (int_time > 2. * min_time) {
          // simple heuristic to skip remaining radii when performance drops
          break;
        }
      }
      switch_radius += 0.025 * maxrad;
    }
    switch_radius = min_rad;
    far_switch_radius_sq = Utils::sqr(switch_radius);
  } else if (far_switch_radius_sq <= Utils::sqr(bessel_radii.back())) {
    // this switching radius is too small for our Bessel series
    throw std::runtime_error("MMM1D could not find a reasonable Bessel cutoff");
  }

  m_is_tuned = true;
  on_coulomb_change();
}

#endif // ELECTROSTATICS
