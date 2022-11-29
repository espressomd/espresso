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

#ifdef P3M

#include "electrostatics/elc.hpp"

#include "electrostatics/coulomb.hpp"
#include "electrostatics/mmm-common.hpp"
#include "electrostatics/p3m.hpp"
#include "electrostatics/p3m_gpu.hpp"

#include "Particle.hpp"
#include "ParticleRange.hpp"
#include "cells.hpp"
#include "communication.hpp"
#include "errorhandling.hpp"
#include "event.hpp"
#include "grid.hpp"

#include <utils/Vector.hpp>
#include <utils/constants.hpp>
#include <utils/math/sqr.hpp>

#include <boost/mpi/collectives/all_reduce.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <functional>
#include <vector>

/** \name Product decomposition data organization
 *  For the cell blocks it is assumed that the lower blocks part is in the
 *  lower half. This has to have positive sign, so that has to be first.
 */
/**@{*/
#define POQESP 0
#define POQECP 1
#define POQESM 2
#define POQECM 3

#define PQESSP 0
#define PQESCP 1
#define PQECSP 2
#define PQECCP 3
#define PQESSM 4
#define PQESCM 5
#define PQECSM 6
#define PQECCM 7
/**@}*/

/** ELC axes (x and y directions)*/
enum class PoQ : int { P, Q };
/** ELC charge sum/assign protocol: real charges, image charges, or both. */
enum class ChargeProtocol : int { REAL, IMAGE, BOTH };

/** temporary buffers for product decomposition */
static std::vector<double> partblk;
/** collected data from the other cells */
static double gblcblk[8];

/** structure for caching sin and cos values */
struct SCCache {
  double s, c;
};

/** Cached sin/cos values along the x-axis and y-axis */
/**@{*/
static std::vector<SCCache> scxcache;
static std::vector<SCCache> scycache;
/**@}*/

/**
 * @brief Calculate cached sin/cos values for one direction.
 *
 * @tparam dir Index of the dimension to consider (e.g. 0 for x ...).
 *
 * @param particles Particle to calculate values for
 * @param n_freq Number of frequencies to calculate per particle
 * @param u Inverse box length
 * @return Calculated values.
 */
template <std::size_t dir>
static std::vector<SCCache> calc_sc_cache(ParticleRange const &particles,
                                          std::size_t n_freq, double u) {
  auto constexpr c_2pi = 2. * Utils::pi();
  auto const n_part = particles.size();
  std::vector<SCCache> ret(n_freq * n_part);

  for (std::size_t freq = 1; freq <= n_freq; freq++) {
    auto const pref = c_2pi * u * static_cast<double>(freq);

    std::size_t o = (freq - 1) * n_part;
    for (auto const &p : particles) {
      auto const arg = pref * p.pos()[dir];
      ret[o++] = {sin(arg), cos(arg)};
    }
  }

  return ret;
}

static std::pair<std::size_t, std::size_t>
prepare_sc_cache(ParticleRange const &particles, double far_cut) {
  assert(far_cut >= 0.);
  auto const n_freq_x =
      static_cast<std::size_t>(std::ceil(far_cut * box_geo.length()[0]) + 1.);
  auto const n_freq_y =
      static_cast<std::size_t>(std::ceil(far_cut * box_geo.length()[1]) + 1.);
  auto const u_x = box_geo.length_inv()[0];
  auto const u_y = box_geo.length_inv()[1];
  scxcache = calc_sc_cache<0>(particles, n_freq_x, u_x);
  scycache = calc_sc_cache<1>(particles, n_freq_y, u_y);
  return {n_freq_x, n_freq_y};
}

/*****************************************************************/
/* data distribution */
/*****************************************************************/

static void clear_vec(double *pdc, std::size_t size) {
  for (std::size_t i = 0; i < size; i++)
    pdc[i] = 0.;
}

static void copy_vec(double *pdc_d, double const *pdc_s, std::size_t size) {
  for (std::size_t i = 0; i < size; i++)
    pdc_d[i] = pdc_s[i];
}

static void add_vec(double *pdc_d, double const *pdc_s1, double const *pdc_s2,
                    std::size_t size) {
  for (std::size_t i = 0; i < size; i++)
    pdc_d[i] = pdc_s1[i] + pdc_s2[i];
}

static void addscale_vec(double *pdc_d, double scale, double const *pdc_s1,
                         double const *pdc_s2, std::size_t size) {
  for (std::size_t i = 0; i < size; i++)
    pdc_d[i] = scale * pdc_s1[i] + pdc_s2[i];
}

static void scale_vec(double scale, double *pdc, std::size_t size) {
  for (std::size_t i = 0; i < size; i++)
    pdc[i] *= scale;
}

static double *block(double *p, std::size_t index, std::size_t size) {
  return &p[index * size];
}

static void distribute(std::size_t size) {
  assert(size <= 8);
  double send_buf[8];
  copy_vec(send_buf, gblcblk, size);
  boost::mpi::all_reduce(comm_cart, send_buf, static_cast<int>(size), gblcblk,
                         std::plus<>());
}

void ElectrostaticLayerCorrection::check_gap(Particle const &p) const {
  if (p.q() != 0.) {
    auto const z = p.pos()[2];
    if (z < 0. or z > elc.box_h) {
      runtimeErrorMsg() << "Particle " << p.id() << " entered ELC gap "
                        << "region by " << ((z < 0.) ? z : z - elc.box_h);
    }
  }
}

/*****************************************************************/
/* dipole terms */
/*****************************************************************/

/** Calculate the dipole force.
 *  See @cite yeh99a.
 */
void ElectrostaticLayerCorrection::add_dipole_force(
    ParticleRange const &particles) const {
  constexpr std::size_t size = 3;
  auto const pref = prefactor * 4. * Utils::pi() / box_geo.volume();

  /* for non-neutral systems, this shift gives the background contribution
   * (rsp. for this shift, the DM of the background is zero) */
  auto const shift = box_geo.length_half()[2];

  // collect moments

  gblcblk[0] = 0.; // sum q_i (z_i - L/2)
  gblcblk[1] = 0.; // sum q_i z_i
  gblcblk[2] = 0.; // sum q_i

  for (auto const &p : particles) {
    check_gap(p);
    auto const q = p.q();
    auto const z = p.pos()[2];

    gblcblk[0] += q * (z - shift);
    gblcblk[1] += q * z;
    gblcblk[2] += q;

    if (elc.dielectric_contrast_on) {
      if (z < elc.space_layer) {
        gblcblk[0] += elc.delta_mid_bot * q * (-z - shift);
        gblcblk[2] += elc.delta_mid_bot * q;
      }
      if (z > (elc.box_h - elc.space_layer)) {
        gblcblk[0] += elc.delta_mid_top * q * (2. * elc.box_h - z - shift);
        gblcblk[2] += elc.delta_mid_top * q;
      }
    }
  }

  gblcblk[0] *= pref;
  gblcblk[1] *= pref / elc.box_h * box_geo.length()[2];
  gblcblk[2] *= pref;

  distribute(size);

  // Yeh + Berkowitz dipole term @cite yeh99a
  auto field_tot = gblcblk[0];

  // Constand potential contribution
  if (elc.const_pot) {
    auto const field_induced = gblcblk[1];
    auto const field_applied = elc.pot_diff / elc.box_h;
    field_tot -= field_applied + field_induced;
  }

  for (auto &p : particles) {
    p.force()[2] -= field_tot * p.q();

    if (!elc.neutralize) {
      // SUBTRACT the forces of the P3M homogeneous neutralizing background
      p.force()[2] += gblcblk[2] * p.q() * (p.pos()[2] - shift);
    }
  }
}

/** Calculate the dipole energy.
 *  See @cite yeh99a.
 */
double ElectrostaticLayerCorrection::dipole_energy(
    ParticleRange const &particles) const {
  constexpr std::size_t size = 7;
  auto const pref = prefactor * 2. * Utils::pi() / box_geo.volume();
  auto const lz = box_geo.length()[2];
  /* for nonneutral systems, this shift gives the background contribution
     (rsp. for this shift, the DM of the background is zero) */
  auto const shift = box_geo.length_half()[2];

  // collect moments

  gblcblk[0] = 0.; // sum q_i               primary box
  gblcblk[1] = 0.; // sum q_i               boundary layers
  gblcblk[2] = 0.; // sum q_i (z_i - L/2)   primary box
  gblcblk[3] = 0.; // sum q_i (z_i - L/2)   boundary layers
  gblcblk[4] = 0.; // sum q_i (z_i - L/2)^2 primary box
  gblcblk[5] = 0.; // sum q_i (z_i - L/2)^2 boundary layers
  gblcblk[6] = 0.; // sum q_i z_i           primary box

  for (auto const &p : particles) {
    check_gap(p);
    auto const q = p.q();
    auto const z = p.pos()[2];

    gblcblk[0] += q;
    gblcblk[2] += q * (z - shift);
    gblcblk[4] += q * (Utils::sqr(z - shift));
    gblcblk[6] += q * z;

    if (elc.dielectric_contrast_on) {
      if (z < elc.space_layer) {
        gblcblk[1] += elc.delta_mid_bot * q;
        gblcblk[3] += elc.delta_mid_bot * q * (-z - shift);
        gblcblk[5] += elc.delta_mid_bot * q * (Utils::sqr(-z - shift));
      }
      if (z > (elc.box_h - elc.space_layer)) {
        gblcblk[1] += elc.delta_mid_top * q;
        gblcblk[3] += elc.delta_mid_top * q * (2. * elc.box_h - z - shift);
        gblcblk[5] +=
            elc.delta_mid_top * q * (Utils::sqr(2. * elc.box_h - z - shift));
      }
    }
  }

  distribute(size);

  // Yeh + Berkowitz term @cite yeh99a
  auto energy = 2. * pref * (Utils::sqr(gblcblk[2]) + gblcblk[2] * gblcblk[3]);

  if (!elc.neutralize) {
    // SUBTRACT the energy of the P3M homogeneous neutralizing background
    energy += 2. * pref *
              (-gblcblk[0] * gblcblk[4] -
               (.25 - .5 / 3.) * Utils::sqr(gblcblk[0] * lz));
  }

  if (elc.dielectric_contrast_on) {
    if (elc.const_pot) {
      // zero potential difference contribution
      energy += pref / elc.box_h * lz * Utils::sqr(gblcblk[6]);
      // external potential shift contribution
      energy -= 2. * elc.pot_diff / elc.box_h * gblcblk[6];
    }

    /* counter the P3M homogeneous background contribution to the
       boundaries. We never need that, since a homogeneous background
       spanning the artificial boundary layers is aphysical. */
    energy +=
        pref * (-(gblcblk[1] * gblcblk[4] + gblcblk[0] * gblcblk[5]) -
                (1. - 2. / 3.) * gblcblk[0] * gblcblk[1] * Utils::sqr(lz));
  }

  return this_node == 0 ? energy : 0.;
}

/*****************************************************************/

static auto image_sum_b(double q, double z, double d) {
  auto const shift = box_geo.length_half()[2];
  auto const lz = box_geo.length()[2];
  return q / (1. - d) * (z - 2. * d * lz / (1. - d)) - q * shift / (1. - d);
}

static auto image_sum_t(double q, double z, double d) {
  auto const shift = box_geo.length_half()[2];
  auto const lz = box_geo.length()[2];
  return q / (1. - d) * (z + 2. * d * lz / (1. - d)) - q * shift / (1. - d);
}

double
ElectrostaticLayerCorrection::z_energy(ParticleRange const &particles) const {
  constexpr std::size_t size = 4;
  auto const xy_area_inv = box_geo.length_inv()[0] * box_geo.length_inv()[1];
  auto const pref = prefactor * 2. * Utils::pi() * xy_area_inv;
  auto const delta = elc.delta_mid_top * elc.delta_mid_bot;
  auto const fac_delta_mid_bot = elc.delta_mid_bot / (1. - delta);
  auto const fac_delta_mid_top = elc.delta_mid_top / (1. - delta);
  auto const fac_delta = delta / (1. - delta);

  /* for non-neutral systems, this shift gives the background contribution
   * (rsp. for this shift, the DM of the background is zero) */
  double const shift = box_geo.length_half()[2];

  if (elc.dielectric_contrast_on) {
    if (elc.const_pot) {
      clear_vec(gblcblk, size);
      for (auto const &p : particles) {
        auto const z = p.pos()[2];
        auto const q = p.q();
        gblcblk[0] += q;
        gblcblk[1] += q * (z - shift);
        if (z < elc.space_layer) {
          gblcblk[2] -= elc.delta_mid_bot * q;
          gblcblk[3] -= elc.delta_mid_bot * q * (-z - shift);
        }
        if (z > (elc.box_h - elc.space_layer)) {
          gblcblk[2] += elc.delta_mid_top * q;
          gblcblk[3] += elc.delta_mid_top * q * (2. * elc.box_h - z - shift);
        }
      }
    } else {
      // metallic boundaries
      clear_vec(gblcblk, size);
      auto const h = elc.box_h;
      for (auto const &p : particles) {
        auto const z = p.pos()[2];
        auto const q = p.q();
        gblcblk[0] += q;
        gblcblk[1] += q * (z - shift);
        if (elc.dielectric_contrast_on) {
          if (z < elc.space_layer) {
            gblcblk[2] += fac_delta * (elc.delta_mid_bot + 1.) * q;
            gblcblk[3] += q * (image_sum_b(elc.delta_mid_bot * delta,
                                           -(2. * h + z), delta) +
                               image_sum_b(delta, -(2. * h - z), delta));
          } else {
            gblcblk[2] += fac_delta_mid_bot * (1. + elc.delta_mid_top) * q;
            gblcblk[3] += q * (image_sum_b(elc.delta_mid_bot, -z, delta) +
                               image_sum_b(delta, -(2. * h - z), delta));
          }
          if (z > (h - elc.space_layer)) {
            // note the minus sign here which is required due to |z_i-z_j|
            gblcblk[2] -= fac_delta * (elc.delta_mid_top + 1.) * q;
            gblcblk[3] -=
                q * (image_sum_t(elc.delta_mid_top * delta, 4. * h - z, delta) +
                     image_sum_t(delta, 2. * h + z, delta));
          } else {
            // note the minus sign here which is required due to |z_i-z_j|
            gblcblk[2] -= fac_delta_mid_top * (1. + elc.delta_mid_bot) * q;
            gblcblk[3] -=
                q * (image_sum_t(elc.delta_mid_top, 2. * h - z, delta) +
                     image_sum_t(delta, 2. * h + z, delta));
          }
        }
      }
    }
  }
  distribute(size);

  auto const energy = gblcblk[1] * gblcblk[2] - gblcblk[0] * gblcblk[3];
  return (this_node == 0) ? -pref * energy : 0.;
}

void ElectrostaticLayerCorrection::add_z_force(
    ParticleRange const &particles) const {
  constexpr std::size_t size = 1;
  auto const xy_area_inv = box_geo.length_inv()[0] * box_geo.length_inv()[1];
  auto const pref = prefactor * 2. * Utils::pi() * xy_area_inv;
  auto const delta = elc.delta_mid_top * elc.delta_mid_bot;
  auto const fac_delta_mid_bot = elc.delta_mid_bot / (1. - delta);
  auto const fac_delta_mid_top = elc.delta_mid_top / (1. - delta);
  auto const fac_delta = delta / (1. - delta);

  if (elc.dielectric_contrast_on) {
    if (elc.const_pot) {
      clear_vec(gblcblk, size);
      /* just counter the 2 pi |z| contribution stemming from P3M */
      for (auto const &p : particles) {
        auto const z = p.pos()[2];
        auto const q = p.q();
        if (z < elc.space_layer)
          gblcblk[0] -= elc.delta_mid_bot * q;
        if (z > (elc.box_h - elc.space_layer))
          gblcblk[0] += elc.delta_mid_top * q;
      }
    } else {
      clear_vec(gblcblk, size);
      for (auto const &p : particles) {
        auto const z = p.pos()[2];
        auto const q = p.q();
        if (z < elc.space_layer) {
          gblcblk[0] += fac_delta * (elc.delta_mid_bot + 1.) * q;
        } else {
          gblcblk[0] += fac_delta_mid_bot * (elc.delta_mid_top + 1.) * q;
        }
        if (z > (elc.box_h - elc.space_layer)) {
          // note the minus sign here which is required due to |z_i-z_j|
          gblcblk[0] -= fac_delta * (elc.delta_mid_top + 1.) * q;
        } else {
          // note the minus sign here which is required due to |z_i-z_j|
          gblcblk[0] -= fac_delta_mid_top * (elc.delta_mid_bot + 1.) * q;
        }
      }
    }

    gblcblk[0] *= pref;

    distribute(size);

    for (auto &p : particles) {
      p.force()[2] += gblcblk[0] * p.q();
    }
  }
}

/*****************************************************************/
/* PoQ exp sum */
/*****************************************************************/

/** \name q=0 or p=0 per frequency code */
/**@{*/
template <PoQ axis>
void setup_PoQ(elc_data const &elc, double prefactor, std::size_t index,
               double omega, ParticleRange const &particles) {
  assert(index >= 1);
  constexpr std::size_t size = 4;
  auto const xy_area_inv = box_geo.length_inv()[0] * box_geo.length_inv()[1];
  auto const pref_di = prefactor * 4. * Utils::pi() * xy_area_inv;
  auto const pref = -pref_di / expm1(omega * box_geo.length()[2]);
  double lclimgebot[4], lclimgetop[4], lclimge[4];
  double fac_delta_mid_bot = 1., fac_delta_mid_top = 1., fac_delta = 1.;

  if (elc.dielectric_contrast_on) {
    auto const delta = elc.delta_mid_top * elc.delta_mid_bot;
    auto const fac_elc = 1. / (1. - delta * exp(-omega * 2. * elc.box_h));
    fac_delta_mid_bot = elc.delta_mid_bot * fac_elc;
    fac_delta_mid_top = elc.delta_mid_top * fac_elc;
    fac_delta = fac_delta_mid_bot * elc.delta_mid_top;
  }

  clear_vec(lclimge, size);
  clear_vec(gblcblk, size);
  auto const &sc_cache = (axis == PoQ::P) ? scxcache : scycache;

  std::size_t ic = 0;
  auto const o = (index - 1) * particles.size();
  for (auto const &p : particles) {
    auto const z = p.pos()[2];
    auto const q = p.q();
    auto e = exp(omega * z);

    partblk[size * ic + POQESM] = q * sc_cache[o + ic].s / e;
    partblk[size * ic + POQESP] = q * sc_cache[o + ic].s * e;
    partblk[size * ic + POQECM] = q * sc_cache[o + ic].c / e;
    partblk[size * ic + POQECP] = q * sc_cache[o + ic].c * e;

    add_vec(gblcblk, gblcblk, block(partblk.data(), ic, size), size);

    if (elc.dielectric_contrast_on) {
      if (z < elc.space_layer) { // handle the lower case first
        // negative sign is okay here as the image is located at -z

        e = exp(-omega * z);

        auto const scale = q * elc.delta_mid_bot;

        lclimgebot[POQESM] = sc_cache[o + ic].s / e;
        lclimgebot[POQESP] = sc_cache[o + ic].s * e;
        lclimgebot[POQECM] = sc_cache[o + ic].c / e;
        lclimgebot[POQECP] = sc_cache[o + ic].c * e;

        addscale_vec(gblcblk, scale, lclimgebot, gblcblk, size);

        e = (exp(omega * (-z - 2. * elc.box_h)) * elc.delta_mid_bot +
             exp(omega * (+z - 2. * elc.box_h))) *
            fac_delta;
      } else {
        e = (exp(-omega * z) +
             exp(omega * (z - 2. * elc.box_h)) * elc.delta_mid_top) *
            fac_delta_mid_bot;
      }

      lclimge[POQESP] += q * sc_cache[o + ic].s * e;
      lclimge[POQECP] += q * sc_cache[o + ic].c * e;

      if (z > (elc.box_h - elc.space_layer)) { // handle the upper case now
        e = exp(omega * (2. * elc.box_h - z));

        auto const scale = q * elc.delta_mid_top;

        lclimgetop[POQESM] = sc_cache[o + ic].s / e;
        lclimgetop[POQESP] = sc_cache[o + ic].s * e;
        lclimgetop[POQECM] = sc_cache[o + ic].c / e;
        lclimgetop[POQECP] = sc_cache[o + ic].c * e;

        addscale_vec(gblcblk, scale, lclimgetop, gblcblk, size);

        e = (exp(omega * (+z - 4. * elc.box_h)) * elc.delta_mid_top +
             exp(omega * (-z - 2. * elc.box_h))) *
            fac_delta;
      } else {
        e = (exp(omega * (+z - 2. * elc.box_h)) +
             exp(omega * (-z - 2. * elc.box_h)) * elc.delta_mid_bot) *
            fac_delta_mid_top;
      }

      lclimge[POQESM] += q * sc_cache[o + ic].s * e;
      lclimge[POQECM] += q * sc_cache[o + ic].c * e;
    }

    ++ic;
  }

  scale_vec(pref, gblcblk, size);

  if (elc.dielectric_contrast_on) {
    scale_vec(pref_di, lclimge, size);
    add_vec(gblcblk, gblcblk, lclimge, size);
  }
}

template <PoQ axis> void add_PoQ_force(ParticleRange const &particles) {
  constexpr auto i = static_cast<int>(axis);
  constexpr std::size_t size = 4;

  std::size_t ic = 0;
  for (auto &p : particles) {
    auto &force = p.force();
    force[i] += partblk[size * ic + POQESM] * gblcblk[POQECP] -
                partblk[size * ic + POQECM] * gblcblk[POQESP] +
                partblk[size * ic + POQESP] * gblcblk[POQECM] -
                partblk[size * ic + POQECP] * gblcblk[POQESM];
    force[2] += partblk[size * ic + POQECM] * gblcblk[POQECP] +
                partblk[size * ic + POQESM] * gblcblk[POQESP] -
                partblk[size * ic + POQECP] * gblcblk[POQECM] -
                partblk[size * ic + POQESP] * gblcblk[POQESM];
    ++ic;
  }
}

static double PoQ_energy(double omega, std::size_t n_part) {
  constexpr std::size_t size = 4;

  auto energy = 0.;
  for (std::size_t ic = 0; ic < n_part; ic++) {
    energy += partblk[size * ic + POQECM] * gblcblk[POQECP] +
              partblk[size * ic + POQESM] * gblcblk[POQESP] +
              partblk[size * ic + POQECP] * gblcblk[POQECM] +
              partblk[size * ic + POQESP] * gblcblk[POQESM];
  }

  return energy / omega;
}
/**@}*/

/*****************************************************************/
/* PQ particle blocks */
/*****************************************************************/

/** \name p,q <> 0 per frequency code */
/**@{*/
static void setup_PQ(elc_data const &elc, double prefactor, std::size_t index_p,
                     std::size_t index_q, double omega,
                     ParticleRange const &particles) {
  assert(index_p >= 1);
  assert(index_q >= 1);
  constexpr std::size_t size = 8;
  auto const xy_area_inv = box_geo.length_inv()[0] * box_geo.length_inv()[1];
  auto const pref_di = prefactor * 8 * Utils::pi() * xy_area_inv;
  auto const pref = -pref_di / expm1(omega * box_geo.length()[2]);
  double lclimgebot[8], lclimgetop[8], lclimge[8];
  double fac_delta_mid_bot = 1, fac_delta_mid_top = 1, fac_delta = 1;
  if (elc.dielectric_contrast_on) {
    auto const delta = elc.delta_mid_top * elc.delta_mid_bot;
    auto const fac_elc = 1. / (1. - delta * exp(-omega * 2. * elc.box_h));
    fac_delta_mid_bot = elc.delta_mid_bot * fac_elc;
    fac_delta_mid_top = elc.delta_mid_top * fac_elc;
    fac_delta = fac_delta_mid_bot * elc.delta_mid_top;
  }

  clear_vec(lclimge, size);
  clear_vec(gblcblk, size);

  std::size_t ic = 0;
  auto const ox = (index_p - 1) * particles.size();
  auto const oy = (index_q - 1) * particles.size();
  for (auto const &p : particles) {
    auto const z = p.pos()[2];
    auto const q = p.q();
    auto e = exp(omega * z);

    partblk[size * ic + PQESSM] =
        scxcache[ox + ic].s * scycache[oy + ic].s * q / e;
    partblk[size * ic + PQESCM] =
        scxcache[ox + ic].s * scycache[oy + ic].c * q / e;
    partblk[size * ic + PQECSM] =
        scxcache[ox + ic].c * scycache[oy + ic].s * q / e;
    partblk[size * ic + PQECCM] =
        scxcache[ox + ic].c * scycache[oy + ic].c * q / e;

    partblk[size * ic + PQESSP] =
        scxcache[ox + ic].s * scycache[oy + ic].s * q * e;
    partblk[size * ic + PQESCP] =
        scxcache[ox + ic].s * scycache[oy + ic].c * q * e;
    partblk[size * ic + PQECSP] =
        scxcache[ox + ic].c * scycache[oy + ic].s * q * e;
    partblk[size * ic + PQECCP] =
        scxcache[ox + ic].c * scycache[oy + ic].c * q * e;

    add_vec(gblcblk, gblcblk, block(partblk.data(), ic, size), size);

    if (elc.dielectric_contrast_on) {
      if (z < elc.space_layer) { // handle the lower case first
        // change e to take into account the z position of the images

        e = exp(-omega * z);
        auto const scale = q * elc.delta_mid_bot;

        lclimgebot[PQESSM] = scxcache[ox + ic].s * scycache[oy + ic].s / e;
        lclimgebot[PQESCM] = scxcache[ox + ic].s * scycache[oy + ic].c / e;
        lclimgebot[PQECSM] = scxcache[ox + ic].c * scycache[oy + ic].s / e;
        lclimgebot[PQECCM] = scxcache[ox + ic].c * scycache[oy + ic].c / e;

        lclimgebot[PQESSP] = scxcache[ox + ic].s * scycache[oy + ic].s * e;
        lclimgebot[PQESCP] = scxcache[ox + ic].s * scycache[oy + ic].c * e;
        lclimgebot[PQECSP] = scxcache[ox + ic].c * scycache[oy + ic].s * e;
        lclimgebot[PQECCP] = scxcache[ox + ic].c * scycache[oy + ic].c * e;

        addscale_vec(gblcblk, scale, lclimgebot, gblcblk, size);

        e = (exp(omega * (-z - 2. * elc.box_h)) * elc.delta_mid_bot +
             exp(omega * (+z - 2. * elc.box_h))) *
            fac_delta * q;

      } else {

        e = (exp(-omega * z) +
             exp(omega * (z - 2. * elc.box_h)) * elc.delta_mid_top) *
            fac_delta_mid_bot * q;
      }

      lclimge[PQESSP] += scxcache[ox + ic].s * scycache[oy + ic].s * e;
      lclimge[PQESCP] += scxcache[ox + ic].s * scycache[oy + ic].c * e;
      lclimge[PQECSP] += scxcache[ox + ic].c * scycache[oy + ic].s * e;
      lclimge[PQECCP] += scxcache[ox + ic].c * scycache[oy + ic].c * e;

      if (z > (elc.box_h - elc.space_layer)) { // handle the upper case now

        e = exp(omega * (2. * elc.box_h - z));
        auto const scale = q * elc.delta_mid_top;

        lclimgetop[PQESSM] = scxcache[ox + ic].s * scycache[oy + ic].s / e;
        lclimgetop[PQESCM] = scxcache[ox + ic].s * scycache[oy + ic].c / e;
        lclimgetop[PQECSM] = scxcache[ox + ic].c * scycache[oy + ic].s / e;
        lclimgetop[PQECCM] = scxcache[ox + ic].c * scycache[oy + ic].c / e;

        lclimgetop[PQESSP] = scxcache[ox + ic].s * scycache[oy + ic].s * e;
        lclimgetop[PQESCP] = scxcache[ox + ic].s * scycache[oy + ic].c * e;
        lclimgetop[PQECSP] = scxcache[ox + ic].c * scycache[oy + ic].s * e;
        lclimgetop[PQECCP] = scxcache[ox + ic].c * scycache[oy + ic].c * e;

        addscale_vec(gblcblk, scale, lclimgetop, gblcblk, size);

        e = (exp(omega * (+z - 4. * elc.box_h)) * elc.delta_mid_top +
             exp(omega * (-z - 2. * elc.box_h))) *
            fac_delta * q;

      } else {

        e = (exp(omega * (+z - 2. * elc.box_h)) +
             exp(omega * (-z - 2. * elc.box_h)) * elc.delta_mid_bot) *
            fac_delta_mid_top * q;
      }

      lclimge[PQESSM] += scxcache[ox + ic].s * scycache[oy + ic].s * e;
      lclimge[PQESCM] += scxcache[ox + ic].s * scycache[oy + ic].c * e;
      lclimge[PQECSM] += scxcache[ox + ic].c * scycache[oy + ic].s * e;
      lclimge[PQECCM] += scxcache[ox + ic].c * scycache[oy + ic].c * e;
    }

    ic++;
  }

  scale_vec(pref, gblcblk, size);
  if (elc.dielectric_contrast_on) {
    scale_vec(pref_di, lclimge, size);
    add_vec(gblcblk, gblcblk, lclimge, size);
  }
}

static void add_PQ_force(std::size_t index_p, std::size_t index_q, double omega,
                         const ParticleRange &particles) {
  auto constexpr c_2pi = 2. * Utils::pi();
  auto const pref_x =
      c_2pi * box_geo.length_inv()[0] * static_cast<double>(index_p) / omega;
  auto const pref_y =
      c_2pi * box_geo.length_inv()[1] * static_cast<double>(index_q) / omega;
  constexpr std::size_t size = 8;

  std::size_t ic = 0;
  for (auto &p : particles) {
    auto &force = p.force();
    force[0] += pref_x * (partblk[size * ic + PQESCM] * gblcblk[PQECCP] +
                          partblk[size * ic + PQESSM] * gblcblk[PQECSP] -
                          partblk[size * ic + PQECCM] * gblcblk[PQESCP] -
                          partblk[size * ic + PQECSM] * gblcblk[PQESSP] +
                          partblk[size * ic + PQESCP] * gblcblk[PQECCM] +
                          partblk[size * ic + PQESSP] * gblcblk[PQECSM] -
                          partblk[size * ic + PQECCP] * gblcblk[PQESCM] -
                          partblk[size * ic + PQECSP] * gblcblk[PQESSM]);
    force[1] += pref_y * (partblk[size * ic + PQECSM] * gblcblk[PQECCP] +
                          partblk[size * ic + PQESSM] * gblcblk[PQESCP] -
                          partblk[size * ic + PQECCM] * gblcblk[PQECSP] -
                          partblk[size * ic + PQESCM] * gblcblk[PQESSP] +
                          partblk[size * ic + PQECSP] * gblcblk[PQECCM] +
                          partblk[size * ic + PQESSP] * gblcblk[PQESCM] -
                          partblk[size * ic + PQECCP] * gblcblk[PQECSM] -
                          partblk[size * ic + PQESCP] * gblcblk[PQESSM]);
    force[2] += (partblk[size * ic + PQECCM] * gblcblk[PQECCP] +
                 partblk[size * ic + PQECSM] * gblcblk[PQECSP] +
                 partblk[size * ic + PQESCM] * gblcblk[PQESCP] +
                 partblk[size * ic + PQESSM] * gblcblk[PQESSP] -
                 partblk[size * ic + PQECCP] * gblcblk[PQECCM] -
                 partblk[size * ic + PQECSP] * gblcblk[PQECSM] -
                 partblk[size * ic + PQESCP] * gblcblk[PQESCM] -
                 partblk[size * ic + PQESSP] * gblcblk[PQESSM]);
    ic++;
  }
}

static double PQ_energy(double omega, std::size_t n_part) {
  constexpr std::size_t size = 8;

  auto energy = 0.;
  for (std::size_t ic = 0; ic < n_part; ic++) {
    energy += partblk[size * ic + PQECCM] * gblcblk[PQECCP] +
              partblk[size * ic + PQECSM] * gblcblk[PQECSP] +
              partblk[size * ic + PQESCM] * gblcblk[PQESCP] +
              partblk[size * ic + PQESSM] * gblcblk[PQESSP] +
              partblk[size * ic + PQECCP] * gblcblk[PQECCM] +
              partblk[size * ic + PQECSP] * gblcblk[PQECSM] +
              partblk[size * ic + PQESCP] * gblcblk[PQESCM] +
              partblk[size * ic + PQESSP] * gblcblk[PQESSM];
  }
  return energy / omega;
}
/**@}*/

void ElectrostaticLayerCorrection::add_force(
    ParticleRange const &particles) const {
  auto constexpr c_2pi = 2. * Utils::pi();
  auto const n_freqs = prepare_sc_cache(particles, elc.far_cut);
  auto const n_scxcache = std::get<0>(n_freqs);
  auto const n_scycache = std::get<1>(n_freqs);
  partblk.resize(particles.size() * 8);

  add_dipole_force(particles);
  add_z_force(particles);

  /* the second condition is just for the case of numerical accident */
  for (std::size_t p = 1;
       box_geo.length_inv()[0] * static_cast<double>(p - 1) < elc.far_cut &&
       p <= n_scxcache;
       p++) {
    auto const omega = c_2pi * box_geo.length_inv()[0] * static_cast<double>(p);
    setup_PoQ<PoQ::P>(elc, prefactor, p, omega, particles);
    distribute(4);
    add_PoQ_force<PoQ::P>(particles);
  }

  for (std::size_t q = 1;
       box_geo.length_inv()[1] * static_cast<double>(q - 1) < elc.far_cut &&
       q <= n_scycache;
       q++) {
    auto const omega = c_2pi * box_geo.length_inv()[1] * static_cast<double>(q);
    setup_PoQ<PoQ::Q>(elc, prefactor, q, omega, particles);
    distribute(4);
    add_PoQ_force<PoQ::Q>(particles);
  }

  for (std::size_t p = 1;
       box_geo.length_inv()[0] * static_cast<double>(p - 1) < elc.far_cut &&
       p <= n_scxcache;
       p++) {
    for (std::size_t q = 1;
         Utils::sqr(box_geo.length_inv()[0] * static_cast<double>(p - 1)) +
                 Utils::sqr(box_geo.length_inv()[1] *
                            static_cast<double>(q - 1)) <
             elc.far_cut2 &&
         q <= n_scycache;
         q++) {
      auto const omega =
          c_2pi *
          sqrt(Utils::sqr(box_geo.length_inv()[0] * static_cast<double>(p)) +
               Utils::sqr(box_geo.length_inv()[1] * static_cast<double>(q)));
      setup_PQ(elc, prefactor, p, q, omega, particles);
      distribute(8);
      add_PQ_force(p, q, omega, particles);
    }
  }
}

double ElectrostaticLayerCorrection::calc_energy(
    ParticleRange const &particles) const {
  auto constexpr c_2pi = 2. * Utils::pi();
  auto energy = dipole_energy(particles) + z_energy(particles);
  auto const n_freqs = prepare_sc_cache(particles, elc.far_cut);
  auto const n_scxcache = std::get<0>(n_freqs);
  auto const n_scycache = std::get<1>(n_freqs);

  auto const n_localpart = particles.size();
  partblk.resize(n_localpart * 8);

  /* the second condition is just for the case of numerical accident */
  for (std::size_t p = 1;
       box_geo.length_inv()[0] * static_cast<double>(p - 1) < elc.far_cut &&
       p <= n_scxcache;
       p++) {
    auto const omega = c_2pi * box_geo.length_inv()[0] * static_cast<double>(p);
    setup_PoQ<PoQ::P>(elc, prefactor, p, omega, particles);
    distribute(4);
    energy += PoQ_energy(omega, n_localpart);
  }

  for (std::size_t q = 1;
       box_geo.length_inv()[1] * static_cast<double>(q - 1) < elc.far_cut &&
       q <= n_scycache;
       q++) {
    auto const omega = c_2pi * box_geo.length_inv()[1] * static_cast<double>(q);
    setup_PoQ<PoQ::Q>(elc, prefactor, q, omega, particles);
    distribute(4);
    energy += PoQ_energy(omega, n_localpart);
  }

  for (std::size_t p = 1;
       box_geo.length_inv()[0] * static_cast<double>(p - 1) < elc.far_cut &&
       p <= n_scxcache;
       p++) {
    for (std::size_t q = 1;
         Utils::sqr(box_geo.length_inv()[0] * static_cast<double>(p - 1)) +
                 Utils::sqr(box_geo.length_inv()[1] *
                            static_cast<double>(q - 1)) <
             elc.far_cut2 &&
         q <= n_scycache;
         q++) {
      auto const omega =
          c_2pi *
          sqrt(Utils::sqr(box_geo.length_inv()[0] * static_cast<double>(p)) +
               Utils::sqr(box_geo.length_inv()[1] * static_cast<double>(q)));
      setup_PQ(elc, prefactor, p, q, omega, particles);
      distribute(8);
      energy += PQ_energy(omega, n_localpart);
    }
  }
  /* we count both i<->j and j<->i, so return just half of it */
  return 0.5 * energy;
}

double ElectrostaticLayerCorrection::tune_far_cut() const {
  // Largest reasonable cutoff for far formula
  auto constexpr maximal_far_cut = 50.;
  auto const box_l_x_inv = box_geo.length_inv()[0];
  auto const box_l_y_inv = box_geo.length_inv()[1];
  auto const min_inv_boxl = std::min(box_l_x_inv, box_l_y_inv);
  auto const box_l_z = box_geo.length()[2];
  // adjust lz according to dielectric layer method
  auto const lz =
      (elc.dielectric_contrast_on) ? elc.box_h + elc.space_layer : box_l_z;

  auto tuned_far_cut = min_inv_boxl;
  double err;
  do {
    auto const pref = 2. * Utils::pi() * tuned_far_cut;
    auto const sum = pref + 2. * (box_l_x_inv + box_l_y_inv);
    auto const den = -expm1(-pref * lz);
    auto const num1 = exp(pref * (elc.box_h - lz));
    auto const num2 = exp(-pref * (elc.box_h + lz));

    err = 0.5 / den *
          (num1 * (sum + 1. / (lz - elc.box_h)) / (lz - elc.box_h) +
           num2 * (sum + 1. / (lz + elc.box_h)) / (lz + elc.box_h));

    tuned_far_cut += min_inv_boxl;
  } while (err > elc.maxPWerror and tuned_far_cut < maximal_far_cut);
  if (tuned_far_cut >= maximal_far_cut) {
    throw std::runtime_error("ELC tuning failed: maxPWerror too small");
  }
  return tuned_far_cut - min_inv_boxl;
}

static auto calc_total_charge() {
  auto local_q = 0.;
  for (auto const &p : cell_structure.local_particles()) {
    local_q += p.q();
  }
  return boost::mpi::all_reduce(comm_cart, local_q, std::plus<>());
}

void ElectrostaticLayerCorrection::sanity_checks_periodicity() const {
  if (!box_geo.periodic(0) || !box_geo.periodic(1) || !box_geo.periodic(2)) {
    throw std::runtime_error("ELC: requires periodicity (True, True, True)");
  }
}

void ElectrostaticLayerCorrection::sanity_checks_dielectric_contrasts() const {
  if (elc.dielectric_contrast_on) {
    auto const precision_threshold = std::sqrt(ROUND_ERROR_PREC);
    auto const total_charge = std::abs(calc_total_charge());
    if (total_charge >= precision_threshold) {
      if (elc.const_pot) {
        // Disable this line to make ELC work again with non-neutral systems
        // and metallic boundaries
        throw std::runtime_error("ELC does not currently support non-neutral "
                                 "systems with a dielectric contrast.");
      }
      // ELC with non-neutral systems and no fully metallic boundaries
      // does not work
      throw std::runtime_error("ELC does not work for non-neutral systems and "
                               "non-metallic dielectric contrast.");
    }
  }
}

void ElectrostaticLayerCorrection::adapt_solver() {
  boost::apply_visitor(
      [this](auto &solver) {
        set_prefactor(solver->prefactor);
        solver->p3m.params.epsilon = P3M_EPSILON_METALLIC;
      },
      base_solver);
}

void ElectrostaticLayerCorrection::recalc_box_h() {
  auto const new_box_h = box_geo.length()[2] - elc.gap_size;
  if (new_box_h < 0.) {
    throw std::runtime_error("ELC gap size (" + std::to_string(elc.gap_size) +
                             ") larger than box length in z-direction (" +
                             std::to_string(box_geo.length()[2]) + ")");
  }
  elc.box_h = new_box_h;
}

void ElectrostaticLayerCorrection::recalc_space_layer() {
  if (elc.dielectric_contrast_on) {
    auto const p3m_r_cut = boost::apply_visitor(
        [](auto &solver) { return solver->p3m.params.r_cut; }, base_solver);
    // recalculate the space layer size:
    // 1. set the space_layer to be 1/3 of the gap size, so that box = layer
    elc.space_layer = (1. / 3.) * elc.gap_size;
    // 2. but make sure we don't overlap with the near-field formula
    auto const free_space = elc.gap_size - p3m_r_cut;
    // 3. and make sure the space layer is not bigger than half the actual
    // simulation box, to avoid overlaps
    auto const half_box_h = elc.box_h / 2.;
    auto const max_space_layer = std::min(free_space, half_box_h);
    if (elc.space_layer > max_space_layer) {
      if (max_space_layer <= 0.) {
        throw std::runtime_error("P3M real-space cutoff too large for ELC w/ "
                                 "dielectric contrast");
      }
      elc.space_layer = max_space_layer;
    }
    elc.space_box = elc.gap_size - 2. * elc.space_layer;
  }
}

elc_data::elc_data(double maxPWerror, double gap_size, double far_cut,
                   bool neutralize, double delta_top, double delta_bot,
                   bool with_const_pot, double potential_diff)
    : maxPWerror{maxPWerror}, gap_size{gap_size},
      box_h{box_geo.length()[2] - gap_size}, far_cut{far_cut}, far_cut2{-1.},
      far_calculated{far_cut == -1.}, dielectric_contrast_on{delta_top != 0. or
                                                             delta_bot != 0.},
      const_pot{with_const_pot and dielectric_contrast_on},
      neutralize{neutralize and !dielectric_contrast_on},
      delta_mid_top{std::clamp(delta_top, -1., +1.)}, delta_mid_bot{std::clamp(
                                                          delta_bot, -1., +1.)},
      pot_diff{(with_const_pot) ? potential_diff : 0.},
      // initial setup of parameters, may change later when P3M is finally tuned
      // set the space_layer to be 1/3 of the gap size, so that box = layer
      space_layer{(dielectric_contrast_on) ? gap_size / 3. : 0.},
      space_box{gap_size - ((dielectric_contrast_on) ? 2. * space_layer : 0.)} {

  auto const delta_range = 1. + std::sqrt(ROUND_ERROR_PREC);
  if (far_cut <= 0. and not far_calculated) {
    throw std::domain_error("Parameter 'far_cut' must be > 0");
  }
  if (maxPWerror <= 0.) {
    throw std::domain_error("Parameter 'maxPWerror' must be > 0");
  }
  if (gap_size <= 0.) {
    throw std::domain_error("Parameter 'gap_size' must be > 0");
  }
  if (potential_diff != 0. and not with_const_pot) {
    throw std::invalid_argument(
        "Parameter 'const_pot' must be True when 'pot_diff' is non-zero");
  }
  if (delta_top < -delta_range or delta_top > delta_range) {
    throw std::domain_error(
        "Parameter 'delta_mid_top' must be >= -1 and <= +1");
  }
  if (delta_bot < -delta_range or delta_bot > delta_range) {
    throw std::domain_error(
        "Parameter 'delta_mid_bot' must be >= -1 and <= +1");
  }
  /* Dielectric contrasts: the deltas should be either both -1 or both +1 when
   * no constant potential difference is applied. The case of two non-metallic
   * parallel boundaries can only be treated with a constant potential. */
  if (dielectric_contrast_on and not const_pot and
      (std::fabs(1. - delta_mid_top * delta_mid_bot) < ROUND_ERROR_PREC)) {
    throw std::domain_error("ELC with two parallel metallic boundaries "
                            "requires the const_pot option");
  }
}

ElectrostaticLayerCorrection::ElectrostaticLayerCorrection(
    elc_data &&parameters, BaseSolver &&solver)
    : elc{parameters}, base_solver{solver} {
  adapt_solver();
}

Utils::Vector3d elc_data::get_mi_vector(Utils::Vector3d const &a,
                                        Utils::Vector3d const &b) const {
  return box_geo.get_mi_vector(a, b);
}

static void p3m_assign_image_charge(elc_data const &elc, CoulombP3M &p3m,
                                    double q, Utils::Vector3d const &pos) {
  if (pos[2] < elc.space_layer) {
    auto const q_eff = elc.delta_mid_bot * q;
    p3m.assign_charge(q_eff, {pos[0], pos[1], -pos[2]});
  }
  if (pos[2] > (elc.box_h - elc.space_layer)) {
    auto const q_eff = elc.delta_mid_top * q;
    p3m.assign_charge(q_eff, {pos[0], pos[1], 2. * elc.box_h - pos[2]});
  }
}

template <ChargeProtocol protocol>
void charge_assign(elc_data const &elc, CoulombP3M &solver,
                   ParticleRange const &particles) {
  if (protocol == ChargeProtocol::BOTH or protocol == ChargeProtocol::IMAGE) {
    solver.p3m.inter_weights.reset(solver.p3m.params.cao);
  }
  /* prepare local FFT mesh */
  for (int i = 0; i < solver.p3m.local_mesh.size; i++)
    solver.p3m.rs_mesh[i] = 0.;

  for (auto const &p : particles) {
    if (p.q() != 0.) {
      if (protocol == ChargeProtocol::BOTH or
          protocol == ChargeProtocol::REAL) {
        solver.assign_charge(p.q(), p.pos(), solver.p3m.inter_weights);
      }
      if (protocol == ChargeProtocol::BOTH or
          protocol == ChargeProtocol::IMAGE) {
        p3m_assign_image_charge(elc, solver, p.q(), p.pos());
      }
    }
  }
}

template <ChargeProtocol protocol>
void modify_p3m_sums(elc_data const &elc, CoulombP3M &solver,
                     ParticleRange const &particles) {

  Utils::Vector3d node_sums{};
  for (auto const &p : particles) {
    auto const q = p.q();
    if (q != 0.) {
      auto const z = p.pos()[2];

      if (protocol == ChargeProtocol::BOTH or
          protocol == ChargeProtocol::REAL) {
        node_sums[0] += 1.;
        node_sums[1] += Utils::sqr(q);
        node_sums[2] += q;
      }

      if (protocol == ChargeProtocol::BOTH or
          protocol == ChargeProtocol::IMAGE) {
        if (z < elc.space_layer) {
          node_sums[0] += 1.;
          node_sums[1] += Utils::sqr(elc.delta_mid_bot * q);
          node_sums[2] += elc.delta_mid_bot * q;
        }

        if (z > (elc.box_h - elc.space_layer)) {
          node_sums[0] += 1.;
          node_sums[1] += Utils::sqr(elc.delta_mid_top * q);
          node_sums[2] += elc.delta_mid_top * q;
        }
      }
    }
  }

  auto const tot_sums =
      boost::mpi::all_reduce(comm_cart, node_sums, std::plus<>());
  solver.p3m.sum_qpart = static_cast<int>(tot_sums[0] + 0.1);
  solver.p3m.sum_q2 = tot_sums[1];
  solver.p3m.square_sum_q = Utils::sqr(tot_sums[2]);
}

double ElectrostaticLayerCorrection::long_range_energy(
    ParticleRange const &particles) const {
  auto const energy = boost::apply_visitor(
      [this, &particles](auto const &solver_ptr) {
        auto &solver = *solver_ptr;

        // assign the original charges (they may not have been assigned yet)
        solver.charge_assign(particles);

        if (!elc.dielectric_contrast_on) {
          return solver.long_range_energy(particles);
        }

        auto energy = 0.;
        energy += 0.5 * solver.long_range_energy(particles);
        energy += 0.5 * elc.dielectric_layers_self_energy(solver, particles);

        // assign both original and image charges
        charge_assign<ChargeProtocol::BOTH>(elc, solver, particles);
        modify_p3m_sums<ChargeProtocol::BOTH>(elc, solver, particles);
        energy += 0.5 * solver.long_range_energy(particles);

        // assign only the image charges now
        charge_assign<ChargeProtocol::IMAGE>(elc, solver, particles);
        modify_p3m_sums<ChargeProtocol::IMAGE>(elc, solver, particles);
        energy -= 0.5 * solver.long_range_energy(particles);

        // restore modified sums
        modify_p3m_sums<ChargeProtocol::REAL>(elc, solver, particles);

        return energy;
      },
      base_solver);
  return energy + calc_energy(particles);
}

void ElectrostaticLayerCorrection::add_long_range_forces(
    ParticleRange const &particles) const {
  boost::apply_visitor(
      [this, &particles](auto const &solver_ptr) {
        auto &solver = *solver_ptr;
        if (elc.dielectric_contrast_on) {
          modify_p3m_sums<ChargeProtocol::BOTH>(elc, solver, particles);
          charge_assign<ChargeProtocol::BOTH>(elc, solver, particles);
          elc.dielectric_layers_self_forces(solver, particles);
        } else {
          solver.charge_assign(particles);
        }
        solver.add_long_range_forces(particles);
        if (elc.dielectric_contrast_on) {
          modify_p3m_sums<ChargeProtocol::REAL>(elc, solver, particles);
        }
      },
      base_solver);
  add_force(particles);
}

#endif
