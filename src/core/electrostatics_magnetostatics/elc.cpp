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
 *  Implementation of \ref elc.hpp.
 */

#include "config.hpp"

#ifdef P3M

#include "electrostatics_magnetostatics/elc.hpp"

#include "electrostatics_magnetostatics/common.hpp"
#include "electrostatics_magnetostatics/coulomb.hpp"
#include "electrostatics_magnetostatics/mmm-common.hpp"
#include "electrostatics_magnetostatics/p3m.hpp"

#include "Particle.hpp"
#include "ParticleRange.hpp"
#include "communication.hpp"
#include "errorhandling.hpp"
#include "grid.hpp"

#include <utils/Vector.hpp>
#include <utils/constants.hpp>
#include <utils/math/sqr.hpp>

#include <mpi.h>

#include <cassert>
#include <cmath>
#include <cstddef>
#include <vector>

ELC_struct elc_params = {1e100, 10,    1, 0, true, true, false, 1,
                         1,     false, 0, 0, 0,    0,    0.0};

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

/** temporary buffers for product decomposition */
static std::vector<double> partblk;
/** collected data from the other cells */
static double gblcblk[8];

/** structure for caching sin and cos values */
typedef struct {
  double s, c;
} SCCache;

/** Cached sin/cos values along the x-axis and y-axis */
/**@{*/
static std::vector<SCCache> scxcache;
static std::vector<SCCache> scycache;
/**@}*/

/****************************************
 * LOCAL FUNCTIONS
 ****************************************/

static void distribute(std::size_t size);
static void add_dipole_force(const ParticleRange &particles);
static double dipole_energy(const ParticleRange &particles);
static double z_energy(const ParticleRange &particles);
static void add_z_force(const ParticleRange &particles);

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
template <size_t dir>
static std::vector<SCCache> calc_sc_cache(const ParticleRange &particles,
                                          std::size_t n_freq, double u) {
  constexpr double c_2pi = 2 * Utils::pi();
  auto const n_part = particles.size();
  std::vector<SCCache> ret(n_freq * n_part);

  for (std::size_t freq = 1; freq <= n_freq; freq++) {
    auto const pref = c_2pi * u * static_cast<double>(freq);

    size_t o = (freq - 1) * n_part;
    for (auto const &part : particles) {
      auto const arg = pref * part.r.p[dir];
      ret[o++] = {sin(arg), cos(arg)};
    }
  }

  return ret;
}

static void prepare_sc_cache(const ParticleRange &particles,
                             std::size_t n_freq_x, double u_x,
                             std::size_t n_freq_y, double u_y) {
  scxcache = calc_sc_cache<0>(particles, n_freq_x, u_x);
  scycache = calc_sc_cache<1>(particles, n_freq_y, u_y);
}

/*****************************************************************/
/* data distribution */
/*****************************************************************/

inline void clear_vec(double *pdc, std::size_t size) {
  for (std::size_t i = 0; i < size; i++)
    pdc[i] = 0;
}

inline void copy_vec(double *pdc_d, double const *pdc_s, std::size_t size) {
  for (std::size_t i = 0; i < size; i++)
    pdc_d[i] = pdc_s[i];
}

inline void add_vec(double *pdc_d, double const *pdc_s1, double const *pdc_s2,
                    std::size_t size) {
  for (std::size_t i = 0; i < size; i++)
    pdc_d[i] = pdc_s1[i] + pdc_s2[i];
}

inline void addscale_vec(double *pdc_d, double scale, double const *pdc_s1,
                         double const *pdc_s2, std::size_t size) {
  for (std::size_t i = 0; i < size; i++)
    pdc_d[i] = scale * pdc_s1[i] + pdc_s2[i];
}

inline void scale_vec(double scale, double *pdc, std::size_t size) {
  for (std::size_t i = 0; i < size; i++)
    pdc[i] *= scale;
}

inline double *block(double *p, std::size_t index, std::size_t size) {
  return &p[index * size];
}

void distribute(std::size_t size) {
  assert(size <= 8);
  double send_buf[8];
  copy_vec(send_buf, gblcblk, size);
  MPI_Allreduce(send_buf, gblcblk, static_cast<int>(size), MPI_DOUBLE, MPI_SUM,
                comm_cart);
}

/** Checks if a charged particle is in the forbidden gap region
 */
inline void check_gap_elc(const Particle &p) {
  if (p.p.q != 0) {
    if (p.r.p[2] < 0)
      runtimeErrorMsg() << "Particle " << p.p.identity << " entered ELC gap "
                        << "region by " << (p.r.p[2]);
    else if (p.r.p[2] > elc_params.h) {
      runtimeErrorMsg() << "Particle " << p.p.identity << " entered ELC gap "
                        << "region by " << (p.r.p[2] - elc_params.h);
    }
  }
}

/*****************************************************************/
/* dipole terms */
/*****************************************************************/

/** Calculate the dipole force.
 *  See @cite yeh99a.
 */
static void add_dipole_force(const ParticleRange &particles) {
  double const pref = coulomb.prefactor * 4 * Utils::pi() *
                      box_geo.length_inv()[0] * box_geo.length_inv()[1] *
                      box_geo.length_inv()[2];
  constexpr std::size_t size = 3;

  auto local_particles = particles;

  /* for nonneutral systems, this shift gives the background contribution
     (rsp. for this shift, the DM of the background is zero) */
  double const shift = box_geo.length_half()[2];

  // collect moments

  gblcblk[0] = 0; // sum q_i (z_i - L/2)
  gblcblk[1] = 0; // sum q_i z_i
  gblcblk[2] = 0; // sum q_i

  for (auto const &p : local_particles) {
    check_gap_elc(p);

    gblcblk[0] += p.p.q * (p.r.p[2] - shift);
    gblcblk[1] += p.p.q * p.r.p[2];
    gblcblk[2] += p.p.q;

    if (elc_params.dielectric_contrast_on) {
      if (p.r.p[2] < elc_params.space_layer) {
        gblcblk[0] += elc_params.delta_mid_bot * p.p.q * (-p.r.p[2] - shift);
        gblcblk[2] += elc_params.delta_mid_bot * p.p.q;
      }
      if (p.r.p[2] > (elc_params.h - elc_params.space_layer)) {
        gblcblk[0] += elc_params.delta_mid_top * p.p.q *
                      (2 * elc_params.h - p.r.p[2] - shift);
        gblcblk[2] += elc_params.delta_mid_top * p.p.q;
      }
    }
  }

  gblcblk[0] *= pref;
  gblcblk[1] *= pref / elc_params.h * box_geo.length()[2];
  gblcblk[2] *= pref;

  distribute(size);

  // Yeh + Berkowitz dipole term @cite yeh99a
  double field_tot = gblcblk[0];

  // Const. potential contribution
  if (elc_params.const_pot) {
    coulomb.field_induced = gblcblk[1];
    coulomb.field_applied = elc_params.pot_diff / elc_params.h;
    field_tot -= coulomb.field_applied + coulomb.field_induced;
  }

  for (auto &p : local_particles) {
    p.f.f[2] -= field_tot * p.p.q;

    if (!elc_params.neutralize) {
      // SUBTRACT the forces of the P3M homogeneous neutralizing background
      p.f.f[2] += gblcblk[2] * p.p.q * (p.r.p[2] - shift);
    }
  }
}

/** Calculate the dipole energy.
 *  See @cite yeh99a.
 */
static double dipole_energy(const ParticleRange &particles) {
  double const pref = coulomb.prefactor * 2 * Utils::pi() *
                      box_geo.length_inv()[0] * box_geo.length_inv()[1] *
                      box_geo.length_inv()[2];
  constexpr std::size_t size = 7;
  /* for nonneutral systems, this shift gives the background contribution
     (rsp. for this shift, the DM of the background is zero) */
  double const shift = box_geo.length_half()[2];

  // collect moments

  gblcblk[0] = 0; // sum q_i               primary box
  gblcblk[1] = 0; // sum q_i               boundary layers
  gblcblk[2] = 0; // sum q_i (z_i - L/2)   primary box
  gblcblk[3] = 0; // sum q_i (z_i - L/2)   boundary layers
  gblcblk[4] = 0; // sum q_i (z_i - L/2)^2 primary box
  gblcblk[5] = 0; // sum q_i (z_i - L/2)^2 boundary layers
  gblcblk[6] = 0; // sum q_i z_i           primary box

  for (auto &p : particles) {
    check_gap_elc(p);

    gblcblk[0] += p.p.q;
    gblcblk[2] += p.p.q * (p.r.p[2] - shift);
    gblcblk[4] += p.p.q * (Utils::sqr(p.r.p[2] - shift));
    gblcblk[6] += p.p.q * p.r.p[2];

    if (elc_params.dielectric_contrast_on) {
      if (p.r.p[2] < elc_params.space_layer) {
        gblcblk[1] += elc_params.delta_mid_bot * p.p.q;
        gblcblk[3] += elc_params.delta_mid_bot * p.p.q * (-p.r.p[2] - shift);
        gblcblk[5] +=
            elc_params.delta_mid_bot * p.p.q * (Utils::sqr(-p.r.p[2] - shift));
      }
      if (p.r.p[2] > (elc_params.h - elc_params.space_layer)) {
        gblcblk[1] += elc_params.delta_mid_top * p.p.q;
        gblcblk[3] += elc_params.delta_mid_top * p.p.q *
                      (2 * elc_params.h - p.r.p[2] - shift);
        gblcblk[5] += elc_params.delta_mid_top * p.p.q *
                      (Utils::sqr(2 * elc_params.h - p.r.p[2] - shift));
      }
    }
  }

  distribute(size);

  // Yeh + Berkowitz term @cite yeh99a
  double energy = 2 * pref * (Utils::sqr(gblcblk[2]) + gblcblk[2] * gblcblk[3]);

  if (!elc_params.neutralize) {
    // SUBTRACT the energy of the P3M homogeneous neutralizing background
    energy += 2 * pref *
              (-gblcblk[0] * gblcblk[4] -
               (.25 - .5 / 3.) * Utils::sqr(gblcblk[0] * box_geo.length()[2]));
  }

  if (elc_params.dielectric_contrast_on) {
    if (elc_params.const_pot) {
      // zero potential difference contribution
      energy +=
          pref / elc_params.h * box_geo.length()[2] * Utils::sqr(gblcblk[6]);
      // external potential shift contribution
      energy -= 2 * elc_params.pot_diff / elc_params.h * gblcblk[6];
    }

    /* counter the P3M homogeneous background contribution to the
       boundaries. We never need that, since a homogeneous background
       spanning the artificial boundary layers is aphysical. */
    energy += pref * (-(gblcblk[1] * gblcblk[4] + gblcblk[0] * gblcblk[5]) -
                      (1. - 2. / 3.) * gblcblk[0] * gblcblk[1] *
                          Utils::sqr(box_geo.length()[2]));
  }

  return this_node == 0 ? energy : 0;
}

/*****************************************************************/

inline double image_sum_b(double q, double z) {
  double const shift = box_geo.length_half()[2];
  double const fac = elc_params.delta_mid_top * elc_params.delta_mid_bot;
  double const image_sum =
      (q / (1.0 - fac) * (z - 2.0 * fac * box_geo.length()[2] / (1.0 - fac))) -
      q * shift / (1 - fac);
  return image_sum;
}

inline double image_sum_t(double q, double z) {
  double const shift = box_geo.length_half()[2];
  double const fac = elc_params.delta_mid_top * elc_params.delta_mid_bot;
  double const image_sum =
      (q / (1.0 - fac) * (z + 2.0 * fac * box_geo.length()[2] / (1.0 - fac))) -
      q * shift / (1 - fac);
  return image_sum;
}

/*****************************************************************/
static double z_energy(const ParticleRange &particles) {
  double const pref = coulomb.prefactor * 2 * Utils::pi() *
                      box_geo.length_inv()[0] * box_geo.length_inv()[1];
  constexpr std::size_t size = 4;

  /* for nonneutral systems, this shift gives the background contribution
     (rsp. for this shift, the DM of the background is zero) */
  double const shift = box_geo.length_half()[2];

  if (elc_params.dielectric_contrast_on) {
    if (elc_params.const_pot) {
      clear_vec(gblcblk, size);
      for (auto &p : particles) {
        gblcblk[0] += p.p.q;
        gblcblk[1] += p.p.q * (p.r.p[2] - shift);
        if (p.r.p[2] < elc_params.space_layer) {
          gblcblk[2] -= elc_params.delta_mid_bot * p.p.q;
          gblcblk[3] -= elc_params.delta_mid_bot * p.p.q * (-p.r.p[2] - shift);
        }
        if (p.r.p[2] > (elc_params.h - elc_params.space_layer)) {
          gblcblk[2] += elc_params.delta_mid_top * p.p.q;
          gblcblk[3] += elc_params.delta_mid_top * p.p.q *
                        (2 * elc_params.h - p.r.p[2] - shift);
        }
      }
    } else {
      double const delta = elc_params.delta_mid_top * elc_params.delta_mid_bot;
      double const fac_delta_mid_bot = elc_params.delta_mid_bot / (1 - delta);
      double const fac_delta_mid_top = elc_params.delta_mid_top / (1 - delta);
      double const fac_delta = delta / (1 - delta);

      clear_vec(gblcblk, size);
      for (auto &p : particles) {
        gblcblk[0] += p.p.q;
        gblcblk[1] += p.p.q * (p.r.p[2] - shift);
        if (elc_params.dielectric_contrast_on) {
          if (p.r.p[2] < elc_params.space_layer) {
            gblcblk[2] += fac_delta * (elc_params.delta_mid_bot + 1) * p.p.q;
            gblcblk[3] +=
                p.p.q * (image_sum_b(elc_params.delta_mid_bot * delta,
                                     -(2 * elc_params.h + p.r.p[2])) +
                         image_sum_b(delta, -(2 * elc_params.h - p.r.p[2])));
          } else {
            gblcblk[2] +=
                fac_delta_mid_bot * (1 + elc_params.delta_mid_top) * p.p.q;
            gblcblk[3] +=
                p.p.q * (image_sum_b(elc_params.delta_mid_bot, -p.r.p[2]) +
                         image_sum_b(delta, -(2 * elc_params.h - p.r.p[2])));
          }
          if (p.r.p[2] > (elc_params.h - elc_params.space_layer)) {
            // note the minus sign here which is required due to |z_i-z_j|
            gblcblk[2] -= fac_delta * (elc_params.delta_mid_top + 1) * p.p.q;
            gblcblk[3] -=
                p.p.q * (image_sum_t(elc_params.delta_mid_top * delta,
                                     4 * elc_params.h - p.r.p[2]) +
                         image_sum_t(delta, 2 * elc_params.h + p.r.p[2]));
          } else {
            // note the minus sign here which is required due to |z_i-z_j|
            gblcblk[2] -=
                fac_delta_mid_top * (1 + elc_params.delta_mid_bot) * p.p.q;
            gblcblk[3] -=
                p.p.q * (image_sum_t(elc_params.delta_mid_top,
                                     2 * elc_params.h - p.r.p[2]) +
                         image_sum_t(delta, 2 * elc_params.h + p.r.p[2]));
          }
        }
      }
    }
  }
  distribute(size);

  double energy = 0;
  if (this_node == 0)
    energy -= gblcblk[1] * gblcblk[2] - gblcblk[0] * gblcblk[3];

  return pref * energy;
}

/*****************************************************************/
static void add_z_force(const ParticleRange &particles) {
  double const pref = coulomb.prefactor * 2 * Utils::pi() *
                      box_geo.length_inv()[0] * box_geo.length_inv()[1];
  constexpr std::size_t size = 1;

  if (elc_params.dielectric_contrast_on) {
    auto local_particles = particles;
    if (elc_params.const_pot) {
      clear_vec(gblcblk, size);
      /* just counter the 2 pi |z| contribution stemming from P3M */
      for (auto &p : local_particles) {
        if (p.r.p[2] < elc_params.space_layer)
          gblcblk[0] -= elc_params.delta_mid_bot * p.p.q;
        if (p.r.p[2] > (elc_params.h - elc_params.space_layer))
          gblcblk[0] += elc_params.delta_mid_top * p.p.q;
      }
    } else {
      double const delta = elc_params.delta_mid_top * elc_params.delta_mid_bot;
      double const fac_delta_mid_bot = elc_params.delta_mid_bot / (1 - delta);
      double const fac_delta_mid_top = elc_params.delta_mid_top / (1 - delta);
      double const fac_delta = delta / (1 - delta);

      clear_vec(gblcblk, size);
      for (auto &p : local_particles) {
        if (p.r.p[2] < elc_params.space_layer) {
          gblcblk[0] += fac_delta * (elc_params.delta_mid_bot + 1) * p.p.q;
        } else {
          gblcblk[0] +=
              fac_delta_mid_bot * (1 + elc_params.delta_mid_top) * p.p.q;
        }

        if (p.r.p[2] > (elc_params.h - elc_params.space_layer)) {
          // note the minus sign here which is required due to |z_i-z_j|
          gblcblk[0] -= fac_delta * (elc_params.delta_mid_top + 1) * p.p.q;
        } else {
          // note the minus sign here which is required due to |z_i-z_j|
          gblcblk[0] -=
              fac_delta_mid_top * (1 + elc_params.delta_mid_bot) * p.p.q;
        }
      }
    }

    gblcblk[0] *= pref;

    distribute(size);

    for (auto &p : local_particles) {
      p.f.f[2] += gblcblk[0] * p.p.q;
    }
  }
}

/*****************************************************************/
/* PoQ exp sum */
/*****************************************************************/

/** \name q=0 or p=0 per frequency code */
/**@{*/
template <PoQ axis>
void setup_PoQ(std::size_t index, double omega,
               const ParticleRange &particles) {
  assert(index >= 1);
  double const pref_di = coulomb.prefactor * 4 * Utils::pi() *
                         box_geo.length_inv()[0] * box_geo.length_inv()[1];
  double const pref = -pref_di / expm1(omega * box_geo.length()[2]);
  constexpr std::size_t size = 4;
  double lclimgebot[4], lclimgetop[4], lclimge[4];
  double fac_delta_mid_bot = 1, fac_delta_mid_top = 1, fac_delta = 1;

  if (elc_params.dielectric_contrast_on) {
    double const fac_elc =
        1.0 / (1 - elc_params.delta_mid_top * elc_params.delta_mid_bot *
                       exp(-omega * 2 * elc_params.h));
    fac_delta_mid_bot = elc_params.delta_mid_bot * fac_elc;
    fac_delta_mid_top = elc_params.delta_mid_top * fac_elc;
    fac_delta = fac_delta_mid_bot * elc_params.delta_mid_top;
  }

  clear_vec(lclimge, size);
  clear_vec(gblcblk, size);
  auto &sc_cache = (axis == PoQ::P) ? scxcache : scycache;

  std::size_t ic = 0;
  auto const o = (index - 1) * particles.size();
  for (auto &p : particles) {
    double e = exp(omega * p.r.p[2]);

    partblk[size * ic + POQESM] = p.p.q * sc_cache[o + ic].s / e;
    partblk[size * ic + POQESP] = p.p.q * sc_cache[o + ic].s * e;
    partblk[size * ic + POQECM] = p.p.q * sc_cache[o + ic].c / e;
    partblk[size * ic + POQECP] = p.p.q * sc_cache[o + ic].c * e;

    add_vec(gblcblk, gblcblk, block(partblk.data(), ic, size), size);

    if (elc_params.dielectric_contrast_on) {
      if (p.r.p[2] < elc_params.space_layer) { // handle the lower case first
        // negative sign is okay here as the image is located at -p.r.p[2]

        e = exp(-omega * p.r.p[2]);

        double const scale = p.p.q * elc_params.delta_mid_bot;

        lclimgebot[POQESM] = sc_cache[o + ic].s / e;
        lclimgebot[POQESP] = sc_cache[o + ic].s * e;
        lclimgebot[POQECM] = sc_cache[o + ic].c / e;
        lclimgebot[POQECP] = sc_cache[o + ic].c * e;

        addscale_vec(gblcblk, scale, lclimgebot, gblcblk, size);

        e = (exp(omega * (-p.r.p[2] - 2 * elc_params.h)) *
                 elc_params.delta_mid_bot +
             exp(omega * (p.r.p[2] - 2 * elc_params.h))) *
            fac_delta;

      } else {

        e = (exp(omega * (-p.r.p[2])) +
             exp(omega * (p.r.p[2] - 2 * elc_params.h)) *
                 elc_params.delta_mid_top) *
            fac_delta_mid_bot;
      }

      lclimge[POQESP] += p.p.q * sc_cache[o + ic].s * e;
      lclimge[POQECP] += p.p.q * sc_cache[o + ic].c * e;

      if (p.r.p[2] > (elc_params.h -
                      elc_params.space_layer)) { // handle the upper case now

        e = exp(omega * (2 * elc_params.h - p.r.p[2]));

        double const scale = p.p.q * elc_params.delta_mid_top;

        lclimgetop[POQESM] = sc_cache[o + ic].s / e;
        lclimgetop[POQESP] = sc_cache[o + ic].s * e;
        lclimgetop[POQECM] = sc_cache[o + ic].c / e;
        lclimgetop[POQECP] = sc_cache[o + ic].c * e;

        addscale_vec(gblcblk, scale, lclimgetop, gblcblk, size);

        e = (exp(omega * (+p.r.p[2] - 4 * elc_params.h)) *
                 elc_params.delta_mid_top +
             exp(omega * (-p.r.p[2] - 2 * elc_params.h))) *
            fac_delta;

      } else {

        e = (exp(omega * (+p.r.p[2] - 2 * elc_params.h)) +
             exp(omega * (-p.r.p[2] - 2 * elc_params.h)) *
                 elc_params.delta_mid_bot) *
            fac_delta_mid_top;
      }

      lclimge[POQESM] += p.p.q * sc_cache[o + ic].s * e;
      lclimge[POQECM] += p.p.q * sc_cache[o + ic].c * e;
    }

    ic++;
  }

  scale_vec(pref, gblcblk, size);

  if (elc_params.dielectric_contrast_on) {
    scale_vec(pref_di, lclimge, size);
    add_vec(gblcblk, gblcblk, lclimge, size);
  }
}

template <PoQ axis> void add_PoQ_force(const ParticleRange &particles) {
  constexpr auto i = static_cast<int>(axis);
  constexpr std::size_t size = 4;

  std::size_t ic = 0;
  for (auto &p : particles) {
    p.f.f[i] += partblk[size * ic + POQESM] * gblcblk[POQECP] -
                partblk[size * ic + POQECM] * gblcblk[POQESP] +
                partblk[size * ic + POQESP] * gblcblk[POQECM] -
                partblk[size * ic + POQECP] * gblcblk[POQESM];
    p.f.f[2] += partblk[size * ic + POQECM] * gblcblk[POQECP] +
                partblk[size * ic + POQESM] * gblcblk[POQESP] -
                partblk[size * ic + POQECP] * gblcblk[POQECM] -
                partblk[size * ic + POQESP] * gblcblk[POQESM];
    ic++;
  }
}

static double PoQ_energy(double omega, std::size_t n_part) {
  constexpr std::size_t size = 4;

  double energy = 0;
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
static void setup_PQ(std::size_t index_p, std::size_t index_q, double omega,
                     const ParticleRange &particles) {
  assert(index_p >= 1);
  assert(index_q >= 1);
  double const pref_di = coulomb.prefactor * 8 * Utils::pi() *
                         box_geo.length_inv()[0] * box_geo.length_inv()[1];
  double const pref = -pref_di / expm1(omega * box_geo.length()[2]);
  constexpr std::size_t size = 8;
  double lclimgebot[8], lclimgetop[8], lclimge[8];
  double fac_delta_mid_bot = 1, fac_delta_mid_top = 1, fac_delta = 1;
  if (elc_params.dielectric_contrast_on) {
    double fac_elc =
        1.0 / (1 - elc_params.delta_mid_top * elc_params.delta_mid_bot *
                       exp(-omega * 2 * elc_params.h));
    fac_delta_mid_bot = elc_params.delta_mid_bot * fac_elc;
    fac_delta_mid_top = elc_params.delta_mid_top * fac_elc;
    fac_delta = fac_delta_mid_bot * elc_params.delta_mid_top;
  }

  clear_vec(lclimge, size);
  clear_vec(gblcblk, size);

  std::size_t ic = 0;
  auto const ox = (index_p - 1) * particles.size();
  auto const oy = (index_q - 1) * particles.size();
  for (auto const &p : particles) {
    double e = exp(omega * p.r.p[2]);

    partblk[size * ic + PQESSM] =
        scxcache[ox + ic].s * scycache[oy + ic].s * p.p.q / e;
    partblk[size * ic + PQESCM] =
        scxcache[ox + ic].s * scycache[oy + ic].c * p.p.q / e;
    partblk[size * ic + PQECSM] =
        scxcache[ox + ic].c * scycache[oy + ic].s * p.p.q / e;
    partblk[size * ic + PQECCM] =
        scxcache[ox + ic].c * scycache[oy + ic].c * p.p.q / e;

    partblk[size * ic + PQESSP] =
        scxcache[ox + ic].s * scycache[oy + ic].s * p.p.q * e;
    partblk[size * ic + PQESCP] =
        scxcache[ox + ic].s * scycache[oy + ic].c * p.p.q * e;
    partblk[size * ic + PQECSP] =
        scxcache[ox + ic].c * scycache[oy + ic].s * p.p.q * e;
    partblk[size * ic + PQECCP] =
        scxcache[ox + ic].c * scycache[oy + ic].c * p.p.q * e;

    add_vec(gblcblk, gblcblk, block(partblk.data(), ic, size), size);

    if (elc_params.dielectric_contrast_on) {
      if (p.r.p[2] < elc_params.space_layer) { // handle the lower case first
        // change e to take into account the z position of the images

        e = exp(-omega * p.r.p[2]);
        auto const scale = p.p.q * elc_params.delta_mid_bot;

        lclimgebot[PQESSM] = scxcache[ox + ic].s * scycache[oy + ic].s / e;
        lclimgebot[PQESCM] = scxcache[ox + ic].s * scycache[oy + ic].c / e;
        lclimgebot[PQECSM] = scxcache[ox + ic].c * scycache[oy + ic].s / e;
        lclimgebot[PQECCM] = scxcache[ox + ic].c * scycache[oy + ic].c / e;

        lclimgebot[PQESSP] = scxcache[ox + ic].s * scycache[oy + ic].s * e;
        lclimgebot[PQESCP] = scxcache[ox + ic].s * scycache[oy + ic].c * e;
        lclimgebot[PQECSP] = scxcache[ox + ic].c * scycache[oy + ic].s * e;
        lclimgebot[PQECCP] = scxcache[ox + ic].c * scycache[oy + ic].c * e;

        addscale_vec(gblcblk, scale, lclimgebot, gblcblk, size);

        e = (exp(omega * (-p.r.p[2] - 2 * elc_params.h)) *
                 elc_params.delta_mid_bot +
             exp(omega * (p.r.p[2] - 2 * elc_params.h))) *
            fac_delta * p.p.q;

      } else {

        e = (exp(omega * (-p.r.p[2])) +
             exp(omega * (p.r.p[2] - 2 * elc_params.h)) *
                 elc_params.delta_mid_top) *
            fac_delta_mid_bot * p.p.q;
      }

      lclimge[PQESSP] += scxcache[ox + ic].s * scycache[oy + ic].s * e;
      lclimge[PQESCP] += scxcache[ox + ic].s * scycache[oy + ic].c * e;
      lclimge[PQECSP] += scxcache[ox + ic].c * scycache[oy + ic].s * e;
      lclimge[PQECCP] += scxcache[ox + ic].c * scycache[oy + ic].c * e;

      if (p.r.p[2] > (elc_params.h -
                      elc_params.space_layer)) { // handle the upper case now

        e = exp(omega * (2 * elc_params.h - p.r.p[2]));
        auto const scale = p.p.q * elc_params.delta_mid_top;

        lclimgetop[PQESSM] = scxcache[ox + ic].s * scycache[oy + ic].s / e;
        lclimgetop[PQESCM] = scxcache[ox + ic].s * scycache[oy + ic].c / e;
        lclimgetop[PQECSM] = scxcache[ox + ic].c * scycache[oy + ic].s / e;
        lclimgetop[PQECCM] = scxcache[ox + ic].c * scycache[oy + ic].c / e;

        lclimgetop[PQESSP] = scxcache[ox + ic].s * scycache[oy + ic].s * e;
        lclimgetop[PQESCP] = scxcache[ox + ic].s * scycache[oy + ic].c * e;
        lclimgetop[PQECSP] = scxcache[ox + ic].c * scycache[oy + ic].s * e;
        lclimgetop[PQECCP] = scxcache[ox + ic].c * scycache[oy + ic].c * e;

        addscale_vec(gblcblk, scale, lclimgetop, gblcblk, size);

        e = (exp(omega * (p.r.p[2] - 4 * elc_params.h)) *
                 elc_params.delta_mid_top +
             exp(omega * (-p.r.p[2] - 2 * elc_params.h))) *
            fac_delta * p.p.q;

      } else {

        e = (exp(omega * (p.r.p[2] - 2 * elc_params.h)) +
             exp(omega * (-p.r.p[2] - 2 * elc_params.h)) *
                 elc_params.delta_mid_bot) *
            fac_delta_mid_top * p.p.q;
      }

      lclimge[PQESSM] += scxcache[ox + ic].s * scycache[oy + ic].s * e;
      lclimge[PQESCM] += scxcache[ox + ic].s * scycache[oy + ic].c * e;
      lclimge[PQECSM] += scxcache[ox + ic].c * scycache[oy + ic].s * e;
      lclimge[PQECCM] += scxcache[ox + ic].c * scycache[oy + ic].c * e;
    }

    ic++;
  }

  scale_vec(pref, gblcblk, size);
  if (elc_params.dielectric_contrast_on) {
    scale_vec(pref_di, lclimge, size);
    add_vec(gblcblk, gblcblk, lclimge, size);
  }
}

static void add_PQ_force(std::size_t index_p, std::size_t index_q, double omega,
                         const ParticleRange &particles) {
  constexpr double c_2pi = 2 * Utils::pi();
  double const pref_x =
      c_2pi * box_geo.length_inv()[0] * static_cast<double>(index_p) / omega;
  double const pref_y =
      c_2pi * box_geo.length_inv()[1] * static_cast<double>(index_q) / omega;
  constexpr std::size_t size = 8;

  std::size_t ic = 0;
  for (auto &p : particles) {
    p.f.f[0] += pref_x * (partblk[size * ic + PQESCM] * gblcblk[PQECCP] +
                          partblk[size * ic + PQESSM] * gblcblk[PQECSP] -
                          partblk[size * ic + PQECCM] * gblcblk[PQESCP] -
                          partblk[size * ic + PQECSM] * gblcblk[PQESSP] +
                          partblk[size * ic + PQESCP] * gblcblk[PQECCM] +
                          partblk[size * ic + PQESSP] * gblcblk[PQECSM] -
                          partblk[size * ic + PQECCP] * gblcblk[PQESCM] -
                          partblk[size * ic + PQECSP] * gblcblk[PQESSM]);
    p.f.f[1] += pref_y * (partblk[size * ic + PQECSM] * gblcblk[PQECCP] +
                          partblk[size * ic + PQESSM] * gblcblk[PQESCP] -
                          partblk[size * ic + PQECCM] * gblcblk[PQECSP] -
                          partblk[size * ic + PQESCM] * gblcblk[PQESSP] +
                          partblk[size * ic + PQECSP] * gblcblk[PQECCM] +
                          partblk[size * ic + PQESSP] * gblcblk[PQESCM] -
                          partblk[size * ic + PQECCP] * gblcblk[PQECSM] -
                          partblk[size * ic + PQESCP] * gblcblk[PQESSM]);
    p.f.f[2] += (partblk[size * ic + PQECCM] * gblcblk[PQECCP] +
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

  double energy = 0;
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

/*****************************************************************/
/* main loops */
/*****************************************************************/

void ELC_add_force(const ParticleRange &particles) {
  constexpr double c_2pi = 2 * Utils::pi();
  auto const n_scxcache =
      std::size_t(ceil(elc_params.far_cut * box_geo.length()[0]) + 1);
  auto const n_scycache =
      std::size_t(ceil(elc_params.far_cut * box_geo.length()[1]) + 1);

  prepare_sc_cache(particles, n_scxcache, box_geo.length_inv()[0], n_scycache,
                   box_geo.length_inv()[1]);
  partblk.resize(particles.size() * 8);

  add_dipole_force(particles);
  add_z_force(particles);

  /* the second condition is just for the case of numerical accident */
  for (std::size_t p = 1; box_geo.length_inv()[0] * static_cast<double>(p - 1) <
                              elc_params.far_cut &&
                          p <= n_scxcache;
       p++) {
    auto const omega = c_2pi * box_geo.length_inv()[0] * static_cast<double>(p);
    setup_PoQ<PoQ::P>(p, omega, particles);
    distribute(4);
    add_PoQ_force<PoQ::P>(particles);
  }

  for (std::size_t q = 1; box_geo.length_inv()[1] * static_cast<double>(q - 1) <
                              elc_params.far_cut &&
                          q <= n_scycache;
       q++) {
    auto const omega = c_2pi * box_geo.length_inv()[1] * static_cast<double>(q);
    setup_PoQ<PoQ::Q>(q, omega, particles);
    distribute(4);
    add_PoQ_force<PoQ::Q>(particles);
  }

  for (std::size_t p = 1; box_geo.length_inv()[0] * static_cast<double>(p - 1) <
                              elc_params.far_cut &&
                          p <= n_scxcache;
       p++) {
    for (std::size_t q = 1;
         Utils::sqr(box_geo.length_inv()[0] * static_cast<double>(p - 1)) +
                 Utils::sqr(box_geo.length_inv()[1] *
                            static_cast<double>(q - 1)) <
             elc_params.far_cut2 &&
         q <= n_scycache;
         q++) {
      auto const omega =
          c_2pi *
          sqrt(Utils::sqr(box_geo.length_inv()[0] * static_cast<double>(p)) +
               Utils::sqr(box_geo.length_inv()[1] * static_cast<double>(q)));
      setup_PQ(p, q, omega, particles);
      distribute(8);
      add_PQ_force(p, q, omega, particles);
    }
  }
}

double ELC_energy(const ParticleRange &particles) {
  constexpr double c_2pi = 2 * Utils::pi();
  auto energy = dipole_energy(particles);
  energy += z_energy(particles);

  auto const n_scxcache =
      std::size_t(ceil(elc_params.far_cut * box_geo.length()[0]) + 1);
  auto const n_scycache =
      std::size_t(ceil(elc_params.far_cut * box_geo.length()[1]) + 1);
  prepare_sc_cache(particles, n_scxcache, box_geo.length_inv()[0], n_scycache,
                   box_geo.length_inv()[1]);

  auto const n_localpart = particles.size();
  partblk.resize(n_localpart * 8);

  /* the second condition is just for the case of numerical accident */
  for (std::size_t p = 1; box_geo.length_inv()[0] * static_cast<double>(p - 1) <
                              elc_params.far_cut &&
                          p <= n_scxcache;
       p++) {
    auto const omega = c_2pi * box_geo.length_inv()[0] * static_cast<double>(p);
    setup_PoQ<PoQ::P>(p, omega, particles);
    distribute(4);
    energy += PoQ_energy(omega, n_localpart);
  }

  for (std::size_t q = 1; box_geo.length_inv()[1] * static_cast<double>(q - 1) <
                              elc_params.far_cut &&
                          q <= n_scycache;
       q++) {
    auto const omega = c_2pi * box_geo.length_inv()[1] * static_cast<double>(q);
    setup_PoQ<PoQ::Q>(q, omega, particles);
    distribute(4);
    energy += PoQ_energy(omega, n_localpart);
  }

  for (std::size_t p = 1; box_geo.length_inv()[0] * static_cast<double>(p - 1) <
                              elc_params.far_cut &&
                          p <= n_scxcache;
       p++) {
    for (std::size_t q = 1;
         Utils::sqr(box_geo.length_inv()[0] * static_cast<double>(p - 1)) +
                 Utils::sqr(box_geo.length_inv()[1] *
                            static_cast<double>(q - 1)) <
             elc_params.far_cut2 &&
         q <= n_scycache;
         q++) {
      auto const omega =
          c_2pi *
          sqrt(Utils::sqr(box_geo.length_inv()[0] * static_cast<double>(p)) +
               Utils::sqr(box_geo.length_inv()[1] * static_cast<double>(q)));
      setup_PQ(p, q, omega, particles);
      distribute(8);
      energy += PQ_energy(omega, n_localpart);
    }
  }
  /* we count both i<->j and j<->i, so return just half of it */
  return 0.5 * energy;
}

double ELC_tune_far_cut(ELC_struct const &params) {
  // Largest reasonable cutoff for far formula
  constexpr auto maximal_far_cut = 50.;
  double const h = params.h;
  double lz = box_geo.length()[2];
  double const min_inv_boxl =
      std::min(box_geo.length_inv()[0], box_geo.length_inv()[1]);

  if (params.dielectric_contrast_on) {
    // adjust lz according to dielectric layer method
    lz = params.h + params.space_layer;
  }

  if (h < 0.) {
    throw std::runtime_error("gap size too large");
  }

  auto far_cut = min_inv_boxl;
  double err;
  do {
    const auto prefactor = 2 * Utils::pi() * far_cut;

    const auto sum =
        prefactor + 2 * (box_geo.length_inv()[0] + box_geo.length_inv()[1]);
    const auto den = -expm1(-prefactor * lz);
    const auto num1 = exp(prefactor * (h - lz));
    const auto num2 = exp(-prefactor * (h + lz));

    err = 0.5 / den *
          (num1 * (sum + 1 / (lz - h)) / (lz - h) +
           num2 * (sum + 1 / (lz + h)) / (lz + h));

    far_cut += min_inv_boxl;
  } while (err > params.maxPWerror && far_cut < maximal_far_cut);
  if (far_cut >= maximal_far_cut) {
    throw std::runtime_error("ELC tuning failed: maxPWerror too small");
  }
  return far_cut - min_inv_boxl;
}

/****************************************
 * COMMON PARTS
 ****************************************/

void ELC_sanity_checks(ELC_struct const &params) {
  if (!box_geo.periodic(0) || !box_geo.periodic(1) || !box_geo.periodic(2)) {
    throw std::runtime_error("ELC requires periodicity 1 1 1");
  }
  /* The product of the two dielectric contrasts should be < 1 for ELC to
     work. This is not the case for two parallel boundaries, which can only
     be treated by the constant potential code */
  if (params.dielectric_contrast_on &&
      (fabs(1.0 - params.delta_mid_top * params.delta_mid_bot) <
       ROUND_ERROR_PREC) &&
      !params.const_pot) {
    throw std::runtime_error("ELC with two parallel metallic boundaries "
                             "requires the const_pot option");
  }

  // ELC with non-neutral systems and no fully metallic boundaries does not work
  if (params.dielectric_contrast_on && !params.const_pot &&
      p3m.square_sum_q > ROUND_ERROR_PREC) {
    throw std::runtime_error("ELC does not work for non-neutral systems and "
                             "non-metallic dielectric contrast.");
  }

  // Disable this line to make ELC work again with non-neutral systems and
  // metallic boundaries
  if (params.dielectric_contrast_on && params.const_pot &&
      p3m.square_sum_q > ROUND_ERROR_PREC) {
    throw std::runtime_error("ELC does not currently support non-neutral "
                             "systems with a dielectric contrast.");
  }
}

void ELC_init() {
  elc_params.h = box_geo.length()[2] - elc_params.gap_size;

  if (elc_params.dielectric_contrast_on) {
    // recalculate the space layer size
    // set the space_layer to be 1/3 of the gap size, so that box = layer
    elc_params.space_layer = (1. / 3.) * elc_params.gap_size;
    // but make sure we leave enough space to not have to bother with
    // overlapping realspace P3M
    double maxsl = elc_params.gap_size - p3m.params.r_cut;
    // and make sure the space layer is not bigger than half the actual
    // simulation box, to avoid overlaps
    if (maxsl > .5 * elc_params.h)
      maxsl = .5 * elc_params.h;
    if (elc_params.space_layer > maxsl) {
      if (maxsl <= 0) {
        runtimeErrorMsg() << "P3M real space cutoff too large for ELC w/ "
                             "dielectric contrast";
      } else
        elc_params.space_layer = maxsl;
    }

    // set the space_box
    elc_params.space_box = elc_params.gap_size - 2 * elc_params.space_layer;
    // reset minimal_dist for tuning
    elc_params.minimal_dist =
        std::min(elc_params.space_box, elc_params.space_layer);
  }

  if (elc_params.far_calculated && elc_params.dielectric_contrast_on) {
    try {
      elc_params.far_cut = ELC_tune_far_cut(elc_params);
      elc_params.far_cut2 = Utils::sqr(elc_params.far_cut);
    } catch (std::runtime_error const &err) {
      runtimeErrorMsg() << err.what() << " (during auto-retuning)";
    }
  }
}

void ELC_set_params(double maxPWerror, double gap_size, double far_cut,
                    bool neutralize, double delta_top, double delta_bot,
                    bool const_pot, double pot_diff) {
  assert(coulomb.method == COULOMB_ELC_P3M or coulomb.method == COULOMB_P3M);
  auto const h = box_geo.length()[2] - gap_size;
  if (maxPWerror <= 0.) {
    throw std::domain_error("maxPWerror must be > 0");
  }
  if (gap_size <= 0.) {
    throw std::domain_error("gap_size must be > 0");
  }
  if (h < 0.) {
    throw std::domain_error("gap size too large");
  }

  ELC_struct new_elc_params;
  if (delta_top != 0.0 || delta_bot != 0.0) {
    // setup with dielectric contrast (neutralize is automatic)

    // initial setup of parameters, may change later when P3M is finally tuned
    // set the space_layer to be 1/3 of the gap size, so that box = layer
    auto const space_layer = gap_size / 3.;
    auto const space_box = gap_size - 2. * space_layer;

    new_elc_params = ELC_struct{maxPWerror,
                                far_cut,
                                0.,
                                gap_size,
                                far_cut == -1.,
                                false,
                                true,
                                delta_top,
                                delta_bot,
                                const_pot,
                                (const_pot) ? pot_diff : 0.,
                                std::min(space_box, space_layer),
                                space_layer,
                                space_box,
                                h};
  } else {
    // setup without dielectric contrast
    new_elc_params =
        ELC_struct{maxPWerror, far_cut,  0., gap_size, far_cut == -1.,
                   neutralize, false,    0., 0.,       false,
                   0.,         gap_size, 0., gap_size, h};
  }

  ELC_sanity_checks(new_elc_params);

  if (new_elc_params.far_calculated) {
    new_elc_params.far_cut = ELC_tune_far_cut(new_elc_params);
  }
  new_elc_params.far_cut2 = Utils::sqr(new_elc_params.far_cut);

  // set new parameters
  elc_params = new_elc_params;
  p3m.params.epsilon = P3M_EPSILON_METALLIC;
  coulomb.method = COULOMB_ELC_P3M;

  mpi_bcast_coulomb_params();
}

////////////////////////////////////////////////////////////////////////////////////

void ELC_P3M_self_forces(const ParticleRange &particles) {
  for (auto &p : particles) {
    p.f.f += coulomb.prefactor * ELC_P3M_dielectric_layers_force_contribution(
                                     p.r.p, p.r.p, p.p.q * p.p.q);
  }
}

////////////////////////////////////////////////////////////////////////////////////

namespace {
void assign_image_charge(const Particle &p) {
  if (p.r.p[2] < elc_params.space_layer) {
    auto const q_eff = elc_params.delta_mid_bot * p.p.q;
    auto const pos = Utils::Vector3d{p.r.p[0], p.r.p[1], -p.r.p[2]};

    p3m_assign_charge(q_eff, pos);
  }

  if (p.r.p[2] > (elc_params.h - elc_params.space_layer)) {
    auto const q_eff = elc_params.delta_mid_top * p.p.q;
    auto const pos =
        Utils::Vector3d{p.r.p[0], p.r.p[1], 2 * elc_params.h - p.r.p[2]};

    p3m_assign_charge(q_eff, pos);
  }
}
} // namespace

void ELC_p3m_charge_assign_both(const ParticleRange &particles) {
  p3m.inter_weights.reset(p3m.params.cao);

  /* prepare local FFT mesh */
  for (int i = 0; i < p3m.local_mesh.size; i++)
    p3m.rs_mesh[i] = 0.0;

  for (auto const &p : particles) {
    if (p.p.q != 0.0) {
      p3m_assign_charge(p.p.q, p.r.p, p3m.inter_weights);
      assign_image_charge(p);
    }
  }
}

void ELC_p3m_charge_assign_image(const ParticleRange &particles) {
  /* prepare local FFT mesh */
  for (int i = 0; i < p3m.local_mesh.size; i++)
    p3m.rs_mesh[i] = 0.0;

  for (auto const &p : particles) {
    if (p.p.q != 0.0) {
      assign_image_charge(p);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////////

Utils::Vector3d ELC_P3M_dielectric_layers_force_contribution(
    const Utils::Vector3d &pos1, const Utils::Vector3d &pos2, double q1q2) {
  Utils::Vector3d force{};

  if (pos1[2] < elc_params.space_layer) {
    auto const q = elc_params.delta_mid_bot * q1q2;
    auto const d = box_geo.get_mi_vector(pos2, {pos1[0], pos1[1], -pos1[2]});

    p3m_add_pair_force(q, d, d.norm(), force);
  }

  if (pos1[2] > (elc_params.h - elc_params.space_layer)) {
    auto const q = elc_params.delta_mid_top * q1q2;
    auto const d = box_geo.get_mi_vector(
        pos2, {pos1[0], pos1[1], 2 * elc_params.h - pos1[2]});

    p3m_add_pair_force(q, d, d.norm(), force);
  }

  return force;
}

/////////////////////////////////////////////////////////////////////////////////////

double ELC_P3M_dielectric_layers_energy_contribution(
    Utils::Vector3d const &pos1, Utils::Vector3d const &pos2, double q1q2) {
  double eng = 0.0;

  if (pos1[2] < elc_params.space_layer) {
    auto const q = elc_params.delta_mid_bot * q1q2;

    eng += p3m_pair_energy(
        q, box_geo.get_mi_vector(pos2, {pos1[0], pos1[1], -pos1[2]}).norm());
  }

  if (pos1[2] > (elc_params.h - elc_params.space_layer)) {
    auto const q = elc_params.delta_mid_top * q1q2;
    eng += p3m_pair_energy(
        q,
        box_geo
            .get_mi_vector(pos2, {pos1[0], pos1[1], 2 * elc_params.h - pos1[2]})
            .norm());
  }

  return eng;
}

double ELC_P3M_dielectric_layers_energy_contribution(Particle const &p1,
                                                     Particle const &p2) {
  double eng = 0.0;

  auto const pos1 = p1.r.p;
  auto const pos2 = p2.r.p;
  auto const q1q2 = p1.p.q * p2.p.q;

  eng += ELC_P3M_dielectric_layers_energy_contribution(pos1, pos2, q1q2);
  eng += ELC_P3M_dielectric_layers_energy_contribution(pos2, pos1, q1q2);

  return eng;
}

//////////////////////////////////////////////////////////////////////////////////

double ELC_P3M_dielectric_layers_energy_self(ParticleRange const &particles) {
  double eng = 0.0;

  // Loop cell neighbors
  for (auto const &p : particles) {
    eng += ELC_P3M_dielectric_layers_energy_contribution(p.r.p, p.r.p,
                                                         p.p.q * p.p.q);
  }

  return eng;
}

/////////////////////////////////////////////////////////////////////////////////

void ELC_P3M_modify_p3m_sums_both(ParticleRange const &particles) {
  double node_sums[3], tot_sums[3];

  for (int i = 0; i < 3; i++) {
    node_sums[i] = 0.0;
    tot_sums[i] = 0.0;
  }

  for (auto const &p : particles) {
    if (p.p.q != 0.0) {

      node_sums[0] += 1.0;
      node_sums[1] += Utils::sqr(p.p.q);
      node_sums[2] += p.p.q;

      if (p.r.p[2] < elc_params.space_layer) {

        node_sums[0] += 1.0;
        node_sums[1] += Utils::sqr(elc_params.delta_mid_bot * p.p.q);
        node_sums[2] += elc_params.delta_mid_bot * p.p.q;
      }
      if (p.r.p[2] > (elc_params.h - elc_params.space_layer)) {

        node_sums[0] += 1.0;
        node_sums[1] += Utils::sqr(elc_params.delta_mid_top * p.p.q);
        node_sums[2] += elc_params.delta_mid_top * p.p.q;
      }
    }
  }

  MPI_Allreduce(node_sums, tot_sums, 3, MPI_DOUBLE, MPI_SUM, comm_cart);
  p3m.sum_qpart = (int)(tot_sums[0] + 0.1);
  p3m.sum_q2 = tot_sums[1];
  p3m.square_sum_q = Utils::sqr(tot_sums[2]);
}

void ELC_P3M_modify_p3m_sums_image(ParticleRange const &particles) {
  double node_sums[3], tot_sums[3];

  for (int i = 0; i < 3; i++) {
    node_sums[i] = 0.0;
    tot_sums[i] = 0.0;
  }

  for (auto const &p : particles) {
    if (p.p.q != 0.0) {

      if (p.r.p[2] < elc_params.space_layer) {

        node_sums[0] += 1.0;
        node_sums[1] += Utils::sqr(elc_params.delta_mid_bot * p.p.q);
        node_sums[2] += elc_params.delta_mid_bot * p.p.q;
      }
      if (p.r.p[2] > (elc_params.h - elc_params.space_layer)) {

        node_sums[0] += 1.0;
        node_sums[1] += Utils::sqr(elc_params.delta_mid_top * p.p.q);
        node_sums[2] += elc_params.delta_mid_top * p.p.q;
      }
    }
  }

  MPI_Allreduce(node_sums, tot_sums, 3, MPI_DOUBLE, MPI_SUM, comm_cart);

  p3m.sum_qpart = (int)(tot_sums[0] + 0.1);
  p3m.sum_q2 = tot_sums[1];
  p3m.square_sum_q = Utils::sqr(tot_sums[2]);
}

// this function is required in force.cpp for energy evaluation
void ELC_P3M_restore_p3m_sums(ParticleRange const &particles) {
  double node_sums[3], tot_sums[3];

  for (int i = 0; i < 3; i++) {
    node_sums[i] = 0.0;
    tot_sums[i] = 0.0;
  }

  for (auto const &p : particles) {
    if (p.p.q != 0.0) {

      node_sums[0] += 1.0;
      node_sums[1] += Utils::sqr(p.p.q);
      node_sums[2] += p.p.q;
    }
  }

  MPI_Allreduce(node_sums, tot_sums, 3, MPI_DOUBLE, MPI_SUM, comm_cart);

  p3m.sum_qpart = (int)(tot_sums[0] + 0.1);
  p3m.sum_q2 = tot_sums[1];
  p3m.square_sum_q = Utils::sqr(tot_sums[2]);
}

#endif
