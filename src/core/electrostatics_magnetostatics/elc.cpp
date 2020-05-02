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
#include "Particle.hpp"
#include "communication.hpp"
#include "errorhandling.hpp"
#include "grid.hpp"
#include "mmm-common.hpp"
#include "pressure.hpp"
#include <cmath>
#include <mpi.h>

#include "electrostatics_magnetostatics/coulomb.hpp"
#include "electrostatics_magnetostatics/elc.hpp"
#include "electrostatics_magnetostatics/p3m.hpp"

#ifdef P3M

/****************************************
 * LOCAL DEFINES
 ****************************************/

/** Largest reasonable cutoff for far formula */
#define MAXIMAL_FAR_CUT 50

/****************************************
 * LOCAL VARIABLES
 ****************************************/

/** \name Inverse box dimensions and derived constants */
/*@{*/
static double ux, uy, uz, height_inverse;
/*@}*/

ELC_struct elc_params = {1e100, 10,    1, 0, true, true, false, 1,
                         1,     false, 0, 0, 0,    0,    0.0};

/****************************************
 * LOCAL ARRAYS
 ****************************************/

/** \name Product decomposition data organization
 *  For the cell blocks it is assumed that the lower blocks part is in the
 *  lower half. This has to have positive sign, so that has to be first.
 */
/*@{*/
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
/*@}*/

/** temporary buffers for product decomposition */
static std::vector<double> partblk;
/** collected data from the other cells */
static double gblcblk[8];

/** structure for caching sin and cos values */
typedef struct {
  double s, c;
} SCCache;

/** Cached sin/cos values along the x-axis and y-axis */
/*@{*/
static std::vector<SCCache> scxcache;
static std::vector<SCCache> scycache;
/*@}*/

/****************************************
 * LOCAL FUNCTIONS
 ****************************************/

static void distribute(int size);
/** \name p=0 per frequency code */
/*@{*/
static void setup_P(int p, double omega, const ParticleRange &particles);
static void add_P_force(const ParticleRange &particles);
static double P_energy(double omega, int n_part);
/*@}*/
/** \name q=0 per frequency code */
/*@{*/
static void setup_Q(int q, double omega, const ParticleRange &particles);
static void add_Q_force(const ParticleRange &particles);
static double Q_energy(double omega, int n_part);
/*@}*/
/** \name p,q <> 0 per frequency code */
/*@{*/
static void setup_PQ(int p, int q, double omega,
                     const ParticleRange &particles);
static void add_PQ_force(int p, int q, double omega,
                         const ParticleRange &particles);
static double PQ_energy(double omega, int n_part);
/*@}*/
static void add_dipole_force(const ParticleRange &particles);
static double dipole_energy(const ParticleRange &particles);
static double z_energy(const ParticleRange &particles);
static void add_z_force(const ParticleRange &particles);

void ELC_setup_constants() {
  ux = 1 / box_geo.length()[0];
  uy = 1 / box_geo.length()[1];
  uz = 1 / box_geo.length()[2];

  height_inverse = 1 / elc_params.h;
}

/**
 * @brief Calculated cached sin/cos values for one direction.
 *
 * @tparam dir Index of the dimension to consider (e.g. 0 for x ...).
 *
 * @param particles Particle to calculate values for
 * @param n_freq Number of frequencies to calculate per particle
 * @param u Inverse box length
 * @return Calculated values.
 */
template <size_t dir>
static std::vector<SCCache> sc_cache(const ParticleRange &particles, int n_freq,
                                     double u) {
  auto const n_part = particles.size();
  std::vector<SCCache> ret(n_freq * n_part);

  for (size_t freq = 1; freq <= n_freq; freq++) {
    double pref = C_2PI * u * freq;

    size_t o = (freq - 1) * n_part;
    for (auto const &part : particles) {
      auto const arg = pref * part.r.p[dir];
      ret[o++] = {sin(arg), cos(arg)};
    }
  }

  return ret;
}

static void prepare_sc_cache(const ParticleRange &particles, int n_freq_x,
                             double u_x, int n_freq_y, double u_y) {
  scxcache = sc_cache<0>(particles, n_freq_x, u_x);
  scycache = sc_cache<1>(particles, n_freq_y, u_y);
}

/*****************************************************************/
/* data distribution */
/*****************************************************************/

inline void clear_vec(double *pdc, int size) {
  for (int i = 0; i < size; i++)
    pdc[i] = 0;
}

inline void copy_vec(double *pdc_d, double const *pdc_s, int size) {
  for (int i = 0; i < size; i++)
    pdc_d[i] = pdc_s[i];
}

inline void add_vec(double *pdc_d, double const *pdc_s1, double const *pdc_s2,
                    int size) {
  for (int i = 0; i < size; i++)
    pdc_d[i] = pdc_s1[i] + pdc_s2[i];
}

inline void addscale_vec(double *pdc_d, double scale, double const *pdc_s1,
                         double const *pdc_s2, int size) {
  for (int i = 0; i < size; i++)
    pdc_d[i] = scale * pdc_s1[i] + pdc_s2[i];
}

inline void scale_vec(double scale, double *pdc, int size) {
  for (int i = 0; i < size; i++)
    pdc[i] *= scale;
}

inline double *block(double *p, int index, int size) {
  return &p[index * size];
}

void distribute(int size) {
  double send_buf[8];
  copy_vec(send_buf, gblcblk, size);
  MPI_Allreduce(send_buf, gblcblk, size, MPI_DOUBLE, MPI_SUM, comm_cart);
}

/*****************************************************************/
/* dipole terms */
/*****************************************************************/

/** Calculate the dipole force.
 *  See @cite yeh99a.
 */
static void add_dipole_force(const ParticleRange &particles) {
  double pref = coulomb.prefactor * 4 * M_PI * ux * uy * uz;
  int size = 3;

  auto local_particles = particles;

  /* for nonneutral systems, this shift gives the background contribution
     (rsp. for this shift, the DM of the background is zero) */
  double shift = 0.5 * box_geo.length()[2];
  double field_tot = 0;

  // collect moments

  gblcblk[0] = 0; // sum q_i (z_i - L/2)
  gblcblk[1] = 0; // sum q_i z_i
  gblcblk[2] = 0; // sum q_i

  for (auto const &p : local_particles) {
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
  gblcblk[1] *= pref * height_inverse / uz;
  gblcblk[2] *= pref;

  distribute(size);

  // Yeh + Berkowitz dipole term @cite yeh99a
  field_tot = gblcblk[0];

  // Const. potential contribution
  if (elc_params.const_pot) {
    coulomb.field_induced = gblcblk[1];
    coulomb.field_applied = elc_params.pot_diff * height_inverse;
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
  double pref = coulomb.prefactor * 2 * M_PI * ux * uy * uz;
  int size = 7;
  /* for nonneutral systems, this shift gives the background contribution
     (rsp. for this shift, the DM of the background is zero) */
  double shift = 0.5 * box_geo.length()[2];

  // collect moments

  gblcblk[0] = 0; // sum q_i               primary box
  gblcblk[1] = 0; // sum q_i               boundary layers
  gblcblk[2] = 0; // sum q_i (z_i - L/2)   primary box
  gblcblk[3] = 0; // sum q_i (z_i - L/2)   boundary layers
  gblcblk[4] = 0; // sum q_i (z_i - L/2)^2 primary box
  gblcblk[5] = 0; // sum q_i (z_i - L/2)^2 boundary layers
  gblcblk[6] = 0; // sum q_i z_i           primary box

  for (auto &p : particles) {
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
  double eng = 2 * pref * (Utils::sqr(gblcblk[2]) + gblcblk[2] * gblcblk[3]);

  if (!elc_params.neutralize) {
    // SUBTRACT the energy of the P3M homogeneous neutralizing background
    eng += 2 * pref *
           (-gblcblk[0] * gblcblk[4] -
            (.25 - .5 / 3.) * Utils::sqr(gblcblk[0] * box_geo.length()[2]));
  }

  if (elc_params.dielectric_contrast_on) {
    if (elc_params.const_pot) {
      // zero potential difference contribution
      eng += pref * height_inverse / uz * Utils::sqr(gblcblk[6]);
      // external potential shift contribution
      eng -= 2 * elc_params.pot_diff * height_inverse * gblcblk[6];
    }

    /* counter the P3M homogeneous background contribution to the
       boundaries. We never need that, since a homogeneous background
       spanning the artificial boundary layers is aphysical. */
    eng += pref * (-(gblcblk[1] * gblcblk[4] + gblcblk[0] * gblcblk[5]) -
                   (1. - 2. / 3.) * gblcblk[0] * gblcblk[1] *
                       Utils::sqr(box_geo.length()[2]));
  }

  return this_node == 0 ? eng : 0;
}

/*****************************************************************/

inline double image_sum_b(double q, double z) {
  double shift = 0.5 * box_geo.length()[2];
  double fac = elc_params.delta_mid_top * elc_params.delta_mid_bot;
  double image_sum =
      (q / (1.0 - fac) * (z - 2.0 * fac * box_geo.length()[2] / (1.0 - fac))) -
      q * shift / (1 - fac);
  return image_sum;
}

inline double image_sum_t(double q, double z) {
  double shift = 0.5 * box_geo.length()[2];
  double fac = elc_params.delta_mid_top * elc_params.delta_mid_bot;
  double image_sum =
      (q / (1.0 - fac) * (z + 2.0 * fac * box_geo.length()[2] / (1.0 - fac))) -
      q * shift / (1 - fac);
  return image_sum;
}

/*****************************************************************/
static double z_energy(const ParticleRange &particles) {
  double pref = coulomb.prefactor * 2 * M_PI * ux * uy;
  int size = 4;

  double eng = 0;
  /* for nonneutral systems, this shift gives the background contribution
     (rsp. for this shift, the DM of the background is zero) */
  double shift = 0.5 * box_geo.length()[2];

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
      double delta = elc_params.delta_mid_top * elc_params.delta_mid_bot;
      double fac_delta_mid_bot = elc_params.delta_mid_bot / (1 - delta);
      double fac_delta_mid_top = elc_params.delta_mid_top / (1 - delta);
      double fac_delta = delta / (1 - delta);

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

  if (this_node == 0)
    eng -= pref * (gblcblk[1] * gblcblk[2] - gblcblk[0] * gblcblk[3]);

  return eng;
}

/*****************************************************************/
static void add_z_force(const ParticleRange &particles) {
  double pref = coulomb.prefactor * 2 * M_PI * ux * uy;

  if (elc_params.dielectric_contrast_on) {
    auto local_particles = particles;
    int size = 1;
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
      double delta = elc_params.delta_mid_top * elc_params.delta_mid_bot;
      double fac_delta_mid_bot = elc_params.delta_mid_bot / (1 - delta);
      double fac_delta_mid_top = elc_params.delta_mid_top / (1 - delta);
      double fac_delta = delta / (1 - delta);

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

static void setup_P(int p, double omega, const ParticleRange &particles) {
  double const pref = -coulomb.prefactor * 4 * M_PI * ux * uy /
                      (expm1(omega * box_geo.length()[2]));
  double const pref_di = coulomb.prefactor * 4 * M_PI * ux * uy;
  int const size = 4;
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

  int ic = 0;
  auto const o = static_cast<int>((p - 1) * particles.size());
  for (auto &p : particles) {
    double e = exp(omega * p.r.p[2]);

    partblk[size * ic + POQESM] = p.p.q * scxcache[o + ic].s / e;
    partblk[size * ic + POQESP] = p.p.q * scxcache[o + ic].s * e;
    partblk[size * ic + POQECM] = p.p.q * scxcache[o + ic].c / e;
    partblk[size * ic + POQECP] = p.p.q * scxcache[o + ic].c * e;

    add_vec(gblcblk, gblcblk, block(partblk.data(), ic, size), size);

    if (elc_params.dielectric_contrast_on) {
      if (p.r.p[2] < elc_params.space_layer) { // handle the lower case first
        // negative sign is okay here as the image is located at -p.r.p[2]

        e = exp(-omega * p.r.p[2]);

        double const scale = p.p.q * elc_params.delta_mid_bot;

        lclimgebot[POQESM] = scxcache[o + ic].s / e;
        lclimgebot[POQESP] = scxcache[o + ic].s * e;
        lclimgebot[POQECM] = scxcache[o + ic].c / e;
        lclimgebot[POQECP] = scxcache[o + ic].c * e;

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

      lclimge[POQESP] += p.p.q * scxcache[o + ic].s * e;
      lclimge[POQECP] += p.p.q * scxcache[o + ic].c * e;

      if (p.r.p[2] > (elc_params.h -
                      elc_params.space_layer)) { // handle the upper case now

        e = exp(omega * (2 * elc_params.h - p.r.p[2]));

        double const scale = p.p.q * elc_params.delta_mid_top;

        lclimgetop[POQESM] = scxcache[o + ic].s / e;
        lclimgetop[POQESP] = scxcache[o + ic].s * e;
        lclimgetop[POQECM] = scxcache[o + ic].c / e;
        lclimgetop[POQECP] = scxcache[o + ic].c * e;

        addscale_vec(gblcblk, scale, lclimgetop, gblcblk, size);

        e = (exp(omega * (p.r.p[2] - 4 * elc_params.h)) *
                 elc_params.delta_mid_top +
             exp(omega * (-p.r.p[2] - 2 * elc_params.h))) *
            fac_delta;

      } else {

        e = (exp(omega * (+p.r.p[2] - 2 * elc_params.h)) +
             exp(omega * (-p.r.p[2] - 2 * elc_params.h)) *
                 elc_params.delta_mid_bot) *
            fac_delta_mid_top;
      }

      lclimge[POQESM] += p.p.q * scxcache[o + ic].s * e;
      lclimge[POQECM] += p.p.q * scxcache[o + ic].c * e;
    }

    ic++;
  }

  scale_vec(pref, gblcblk, size);

  if (elc_params.dielectric_contrast_on) {
    scale_vec(pref_di, lclimge, size);
    add_vec(gblcblk, gblcblk, lclimge, size);
  }
}

static void setup_Q(int q, double omega, const ParticleRange &particles) {
  double const pref = -coulomb.prefactor * 4 * M_PI * ux * uy /
                      (expm1(omega * box_geo.length()[2]));
  double const pref_di = coulomb.prefactor * 4 * M_PI * ux * uy;
  int const size = 4;
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

  int ic = 0;
  auto const o = static_cast<int>((q - 1) * particles.size());
  for (auto &p : particles) {
    double e = exp(omega * p.r.p[2]);

    partblk[size * ic + POQESM] = p.p.q * scycache[o + ic].s / e;
    partblk[size * ic + POQESP] = p.p.q * scycache[o + ic].s * e;
    partblk[size * ic + POQECM] = p.p.q * scycache[o + ic].c / e;
    partblk[size * ic + POQECP] = p.p.q * scycache[o + ic].c * e;

    add_vec(gblcblk, gblcblk, block(partblk.data(), ic, size), size);

    if (elc_params.dielectric_contrast_on) {
      if (p.r.p[2] < elc_params.space_layer) { // handle the lower case first
        // negative sign before omega is okay here as the image is located
        // at -p.r.p[2]

        e = exp(-omega * p.r.p[2]);

        double const scale = p.p.q * elc_params.delta_mid_bot;

        lclimgebot[POQESM] = scycache[o + ic].s / e;
        lclimgebot[POQESP] = scycache[o + ic].s * e;
        lclimgebot[POQECM] = scycache[o + ic].c / e;
        lclimgebot[POQECP] = scycache[o + ic].c * e;

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

      lclimge[POQESP] += p.p.q * scycache[o + ic].s * e;
      lclimge[POQECP] += p.p.q * scycache[o + ic].c * e;

      if (p.r.p[2] > (elc_params.h -
                      elc_params.space_layer)) { // handle the upper case now

        e = exp(omega * (2 * elc_params.h - p.r.p[2]));

        double const scale = p.p.q * elc_params.delta_mid_top;

        lclimgetop[POQESM] = scycache[o + ic].s / e;
        lclimgetop[POQESP] = scycache[o + ic].s * e;
        lclimgetop[POQECM] = scycache[o + ic].c / e;
        lclimgetop[POQECP] = scycache[o + ic].c * e;

        addscale_vec(gblcblk, scale, lclimgetop, gblcblk, size);

        e = (exp(omega * (p.r.p[2] - 4 * elc_params.h)) *
                 elc_params.delta_mid_top +
             exp(omega * (-p.r.p[2] - 2 * elc_params.h))) *
            fac_delta;

      } else {

        e = (exp(omega * (p.r.p[2] - 2 * elc_params.h)) +
             exp(omega * (-p.r.p[2] - 2 * elc_params.h)) *
                 elc_params.delta_mid_bot) *
            fac_delta_mid_top;
      }

      lclimge[POQESM] += p.p.q * scycache[o + ic].s * e;
      lclimge[POQECM] += p.p.q * scycache[o + ic].c * e;
    }

    ic++;
  }

  scale_vec(pref, gblcblk, size);

  if (elc_params.dielectric_contrast_on) {
    scale_vec(pref_di, lclimge, size);
    add_vec(gblcblk, gblcblk, lclimge, size);
  }
}

static void add_P_force(const ParticleRange &particles) {
  int ic;
  int size = 4;

  ic = 0;
  for (auto &p : particles) {
    p.f.f[0] += partblk[size * ic + POQESM] * gblcblk[POQECP] -
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

static double P_energy(double omega, int n_part) {
  int size = 4;
  double eng = 0;
  double pref = 1 / omega;

  for (unsigned ic = 0; ic < n_part; ic++) {
    eng += pref * (partblk[size * ic + POQECM] * gblcblk[POQECP] +
                   partblk[size * ic + POQESM] * gblcblk[POQESP] +
                   partblk[size * ic + POQECP] * gblcblk[POQECM] +
                   partblk[size * ic + POQESP] * gblcblk[POQESM]);
  }

  return eng;
}

static void add_Q_force(const ParticleRange &particles) {
  int ic;
  int size = 4;

  ic = 0;
  for (auto &p : particles) {
    p.f.f[1] += partblk[size * ic + POQESM] * gblcblk[POQECP] -
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

static double Q_energy(double omega, int n_part) {
  int size = 4;
  double eng = 0;
  double pref = 1 / omega;

  for (unsigned ic = 0; ic < n_part; ic++) {
    eng += pref * (partblk[size * ic + POQECM] * gblcblk[POQECP] +
                   partblk[size * ic + POQESM] * gblcblk[POQESP] +
                   partblk[size * ic + POQECP] * gblcblk[POQECM] +
                   partblk[size * ic + POQESP] * gblcblk[POQESM]);
  }
  return eng;
}

/*****************************************************************/
/* PQ particle blocks */
/*****************************************************************/

static void setup_PQ(int p, int q, double omega,
                     const ParticleRange &particles) {
  double const pref = -coulomb.prefactor * 8 * M_PI * ux * uy /
                      (expm1(omega * box_geo.length()[2]));
  double const pref_di = coulomb.prefactor * 8 * M_PI * ux * uy;
  int const size = 8;
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

  int ic = 0;
  auto const ox = static_cast<int>((p - 1) * particles.size());
  auto const oy = static_cast<int>((q - 1) * particles.size());
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

static void add_PQ_force(int p, int q, double omega,
                         const ParticleRange &particles) {
  int ic;

  double pref_x = C_2PI * ux * p / omega;
  double pref_y = C_2PI * uy * q / omega;
  int size = 8;

  ic = 0;
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

static double PQ_energy(double omega, int n_part) {
  int size = 8;
  double eng = 0;
  double pref = 1 / omega;

  for (unsigned ic = 0; ic < n_part; ic++) {
    eng += pref * (partblk[size * ic + PQECCM] * gblcblk[PQECCP] +
                   partblk[size * ic + PQECSM] * gblcblk[PQECSP] +
                   partblk[size * ic + PQESCM] * gblcblk[PQESCP] +
                   partblk[size * ic + PQESSM] * gblcblk[PQESSP] +
                   partblk[size * ic + PQECCP] * gblcblk[PQECCM] +
                   partblk[size * ic + PQECSP] * gblcblk[PQECSM] +
                   partblk[size * ic + PQESCP] * gblcblk[PQESCM] +
                   partblk[size * ic + PQESSP] * gblcblk[PQESSM]);
  }
  return eng;
}

/*****************************************************************/
/* main loops */
/*****************************************************************/

void ELC_add_force(const ParticleRange &particles) {
  int p, q;
  double omega;

  auto const n_scxcache = int(ceil(elc_params.far_cut / ux) + 1);
  auto const n_scycache = int(ceil(elc_params.far_cut / uy) + 1);

  prepare_sc_cache(particles, n_scxcache, ux, n_scycache, uy);
  partblk.resize(particles.size() * 8);

  add_dipole_force(particles);
  add_z_force(particles);

  /* the second condition is just for the case of numerical accident */
  for (p = 1; ux * (p - 1) < elc_params.far_cut && p <= n_scxcache; p++) {
    omega = C_2PI * ux * p;
    setup_P(p, omega, particles);
    distribute(4);
    add_P_force(particles);
  }

  for (q = 1; uy * (q - 1) < elc_params.far_cut && q <= n_scycache; q++) {
    omega = C_2PI * uy * q;
    setup_Q(q, omega, particles);
    distribute(4);
    add_Q_force(particles);
  }

  for (p = 1; ux * (p - 1) < elc_params.far_cut && p <= n_scxcache; p++) {
    for (q = 1; Utils::sqr(ux * (p - 1)) + Utils::sqr(uy * (q - 1)) <
                    elc_params.far_cut2 &&
                q <= n_scycache;
         q++) {
      omega = C_2PI * sqrt(Utils::sqr(ux * p) + Utils::sqr(uy * q));
      setup_PQ(p, q, omega, particles);
      distribute(8);
      add_PQ_force(p, q, omega, particles);
    }
  }
}

double ELC_energy(const ParticleRange &particles) {
  double eng;
  int p, q;
  double omega;

  eng = dipole_energy(particles);
  eng += z_energy(particles);

  auto const n_scxcache = int(ceil(elc_params.far_cut / ux) + 1);
  auto const n_scycache = int(ceil(elc_params.far_cut / uy) + 1);
  prepare_sc_cache(particles, n_scxcache, ux, n_scycache, uy);

  auto const n_localpart = particles.size();
  partblk.resize(n_localpart * 8);

  /* the second condition is just for the case of numerical accident */
  for (p = 1; ux * (p - 1) < elc_params.far_cut && p <= n_scxcache; p++) {
    omega = C_2PI * ux * p;
    setup_P(p, omega, particles);
    distribute(4);
    eng += P_energy(omega, n_localpart);
  }
  for (q = 1; uy * (q - 1) < elc_params.far_cut && q <= n_scycache; q++) {
    omega = C_2PI * uy * q;
    setup_Q(q, omega, particles);
    distribute(4);
    eng += Q_energy(omega, n_localpart);
  }
  for (p = 1; ux * (p - 1) < elc_params.far_cut && p <= n_scxcache; p++) {
    for (q = 1; Utils::sqr(ux * (p - 1)) + Utils::sqr(uy * (q - 1)) <
                    elc_params.far_cut2 &&
                q <= n_scycache;
         q++) {
      omega = C_2PI * sqrt(Utils::sqr(ux * p) + Utils::sqr(uy * q));
      setup_PQ(p, q, omega, particles);
      distribute(8);
      eng += PQ_energy(omega, n_localpart);
    }
  }
  /* we count both i<->j and j<->i, so return just half of it */
  return 0.5 * eng;
}

int ELC_tune(double error) {
  double err;
  double h = elc_params.h, lz = box_geo.length()[2];
  double min_inv_boxl = std::min(ux, uy);

  if (elc_params.dielectric_contrast_on) {
    // adjust lz according to dielectric layer method
    lz = elc_params.h + elc_params.space_layer;
  }

  if (h < 0)
    return ES_ERROR;

  elc_params.far_cut = min_inv_boxl;

  do {
    const auto prefactor = 2 * Utils::pi() * elc_params.far_cut;

    const auto sum = prefactor + 2 * (ux + uy);
    const auto den = -expm1(-prefactor * lz);
    const auto num1 = exp(prefactor * (h - lz));
    const auto num2 = exp(-prefactor * (h + lz));

    err = 0.5 / den *
          (num1 * (sum + 1 / (lz - h)) / (lz - h) +
           num2 * (sum + 1 / (lz + h)) / (lz + h));

    elc_params.far_cut += min_inv_boxl;
  } while (err > error && elc_params.far_cut < MAXIMAL_FAR_CUT);
  if (elc_params.far_cut >= MAXIMAL_FAR_CUT)
    return ES_ERROR;
  elc_params.far_cut -= min_inv_boxl;
  elc_params.far_cut2 = Utils::sqr(elc_params.far_cut);

  return ES_OK;
}

/****************************************
 * COMMON PARTS
 ****************************************/

int ELC_sanity_checks() {
  if (!box_geo.periodic(0) || !box_geo.periodic(1) || !box_geo.periodic(2)) {
    runtimeErrorMsg() << "ELC requires periodicity 1 1 1";
    return ES_ERROR;
  }
  /* The product of the two dielectric contrasts should be < 1 for ELC to
     work. This is not the case for two parallel boundaries, which can only
     be treated by the constant potential code */
  if (elc_params.dielectric_contrast_on &&
      (fabs(1.0 - elc_params.delta_mid_top * elc_params.delta_mid_bot) <
       ROUND_ERROR_PREC) &&
      !elc_params.const_pot) {
    runtimeErrorMsg() << "ELC with two parallel metallic boundaries requires "
                         "the const_pot option";
    return ES_ERROR;
  }

  // ELC with non-neutral systems and no fully metallic boundaries does not work
  if (elc_params.dielectric_contrast_on && !elc_params.const_pot &&
      p3m.square_sum_q > ROUND_ERROR_PREC) {
    runtimeErrorMsg() << "ELC does not work for non-neutral systems and "
                         "non-metallic dielectric contrast.";
    return ES_ERROR;
  }

  // Disable this line to make ELC work again with non-neutral systems and
  // metallic boundaries
  if (elc_params.dielectric_contrast_on && elc_params.const_pot &&
      p3m.square_sum_q > ROUND_ERROR_PREC) {
    runtimeErrorMsg() << "ELC does not currently support non-neutral "
                         "systems with a dielectric contrast.";
    return ES_ERROR;
  }

  return ES_OK;
}

void ELC_init() {

  ELC_setup_constants();

  if (elc_params.dielectric_contrast_on) {
    // recalculate the space layer size
    // set the space_layer to be 1/3 of the gap size, so that box = layer
    elc_params.space_layer = (1. / 3.) * elc_params.gap_size;
    // but make sure we leave enough space to not have to bother with
    // overlapping
    // realspace P3M
    double maxsl = elc_params.gap_size - p3m.params.r_cut;
    // and make sure the space layer is not bigger than half the actual
    // simulation box,
    // to avoid overlaps
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

  if (elc_params.far_calculated && (elc_params.dielectric_contrast_on)) {
    if (ELC_tune(elc_params.maxPWerror) == ES_ERROR) {
      runtimeErrorMsg() << "ELC auto-retuning failed, gap size too small";
    }
  }
  if (elc_params.dielectric_contrast_on) {
    p3m.params.additional_mesh[0] = 0;
    p3m.params.additional_mesh[1] = 0;
    p3m.params.additional_mesh[2] = elc_params.space_layer;
  } else {
    p3m.params.additional_mesh[0] = 0;
    p3m.params.additional_mesh[1] = 0;
    p3m.params.additional_mesh[2] = 0;
  }
}

int ELC_set_params(double maxPWerror, double gap_size, double far_cut,
                   bool neutralize, double delta_top, double delta_bot,
                   bool const_pot, double pot_diff) {
  elc_params.maxPWerror = maxPWerror;
  elc_params.gap_size = gap_size;
  elc_params.h = box_geo.length()[2] - gap_size;

  if (delta_top != 0.0 || delta_bot != 0.0) {
    elc_params.dielectric_contrast_on = true;

    elc_params.delta_mid_top = delta_top;
    elc_params.delta_mid_bot = delta_bot;

    // neutralize is automatic with dielectric contrast
    elc_params.neutralize = false;
    // initial setup of parameters, may change later when P3M is finally tuned
    // set the space_layer to be 1/3 of the gap size, so that box = layer
    elc_params.space_layer = (1. / 3.) * gap_size;
    // set the space_box
    elc_params.space_box = gap_size - 2 * elc_params.space_layer;
    // reset minimal_dist for tuning
    elc_params.minimal_dist =
        std::min(elc_params.space_box, elc_params.space_layer);

    // Constant potential parameter setup
    if (const_pot) {
      elc_params.const_pot = true;
      elc_params.pot_diff = pot_diff;
    }
  } else {
    // setup without dielectric contrast
    elc_params.dielectric_contrast_on = false;
    elc_params.const_pot = false;
    elc_params.delta_mid_top = 0;
    elc_params.delta_mid_bot = 0;
    elc_params.neutralize = neutralize;
    elc_params.space_layer = 0;
    elc_params.space_box = elc_params.minimal_dist = gap_size;
  }

  ELC_setup_constants();

  Coulomb::elc_sanity_check();

  elc_params.far_cut = far_cut;
  if (far_cut != -1) {
    elc_params.far_cut2 = Utils::sqr(far_cut);
    elc_params.far_calculated = false;
  } else {
    elc_params.far_calculated = true;
    if (ELC_tune(elc_params.maxPWerror) == ES_ERROR) {
      runtimeErrorMsg() << "ELC tuning failed, gap size too small";
    }
  }
  mpi_bcast_coulomb_params();

  return ES_OK;
}

////////////////////////////////////////////////////////////////////////////////////

void ELC_P3M_self_forces(const ParticleRange &particles) {
  Utils::Vector3d pos;
  double q;

  for (auto &p : particles) {
    if (p.r.p[2] < elc_params.space_layer) {
      q = elc_params.delta_mid_bot * p.p.q * p.p.q;

      pos[0] = p.r.p[0];
      pos[1] = p.r.p[1];
      pos[2] = -p.r.p[2];
      auto const d = get_mi_vector(p.r.p, pos, box_geo);

      p3m_add_pair_force(q, d, d.norm(), p.f.f);
    }
    if (p.r.p[2] > (elc_params.h - elc_params.space_layer)) {
      q = elc_params.delta_mid_top * p.p.q * p.p.q;
      pos[0] = p.r.p[0];
      pos[1] = p.r.p[1];
      pos[2] = 2 * elc_params.h - p.r.p[2];
      auto const d = get_mi_vector(p.r.p, pos, box_geo);

      p3m_add_pair_force(q, d, d.norm(), p.f.f);
    }
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

  for (auto &p : particles) {
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

  for (auto &p : particles) {
    if (p.p.q != 0.0) {
      assign_image_charge(p);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////////

void ELC_P3M_dielectric_layers_force_contribution(Particle const &p1,
                                                  Particle const &p2,
                                                  Utils::Vector3d &force1,
                                                  Utils::Vector3d &force2) {
  Utils::Vector3d pos;
  double q;

  if (p1.r.p[2] < elc_params.space_layer) {
    q = elc_params.delta_mid_bot * p1.p.q * p2.p.q;
    pos[0] = p1.r.p[0];
    pos[1] = p1.r.p[1];
    pos[2] = -p1.r.p[2];
    auto const d = get_mi_vector(p2.r.p, pos, box_geo);

    p3m_add_pair_force(q, d, d.norm(), force2);
  }

  if (p1.r.p[2] > (elc_params.h - elc_params.space_layer)) {
    q = elc_params.delta_mid_top * p1.p.q * p2.p.q;
    pos[0] = p1.r.p[0];
    pos[1] = p1.r.p[1];
    pos[2] = 2 * elc_params.h - p1.r.p[2];
    auto const d = get_mi_vector(p2.r.p, pos, box_geo);

    p3m_add_pair_force(q, d, d.norm(), force2);
  }

  if (p2.r.p[2] < elc_params.space_layer) {
    q = elc_params.delta_mid_bot * p1.p.q * p2.p.q;
    pos[0] = p2.r.p[0];
    pos[1] = p2.r.p[1];
    pos[2] = -p2.r.p[2];
    auto const d = get_mi_vector(p1.r.p, pos, box_geo);

    p3m_add_pair_force(q, d, d.norm(), force1);
  }

  if (p2.r.p[2] > (elc_params.h - elc_params.space_layer)) {
    q = elc_params.delta_mid_top * p1.p.q * p2.p.q;
    pos[0] = p2.r.p[0];
    pos[1] = p2.r.p[1];
    pos[2] = 2 * elc_params.h - p2.r.p[2];
    auto const d = get_mi_vector(p1.r.p, pos, box_geo);

    p3m_add_pair_force(q, d, d.norm(), force1);
  }
}

/////////////////////////////////////////////////////////////////////////////////////

double ELC_P3M_dielectric_layers_energy_contribution(Particle const &p1,
                                                     Particle const &p2) {
  Utils::Vector3d pos;
  double q;
  double tp2;
  double eng = 0.0;

  tp2 = p2.r.p[2];

  if (p1.r.p[2] < elc_params.space_layer) {
    q = elc_params.delta_mid_bot * p1.p.q * p2.p.q;
    pos[0] = p1.r.p[0];
    pos[1] = p1.r.p[1];
    pos[2] = -p1.r.p[2];

    eng += p3m_pair_energy(q, get_mi_vector(p2.r.p, pos, box_geo).norm());
  }

  if (p1.r.p[2] > (elc_params.h - elc_params.space_layer)) {
    q = elc_params.delta_mid_top * p1.p.q * p2.p.q;
    pos[0] = p1.r.p[0];
    pos[1] = p1.r.p[1];
    pos[2] = 2 * elc_params.h - p1.r.p[2];

    eng += p3m_pair_energy(q, get_mi_vector(p2.r.p, pos, box_geo).norm());
  }

  if (tp2 < elc_params.space_layer) {
    q = elc_params.delta_mid_bot * p1.p.q * p2.p.q;
    pos[0] = p2.r.p[0];
    pos[1] = p2.r.p[1];
    pos[2] = -tp2;

    eng += p3m_pair_energy(q, get_mi_vector(p1.r.p, pos, box_geo).norm());
  }

  if (tp2 > (elc_params.h - elc_params.space_layer)) {
    q = elc_params.delta_mid_top * p1.p.q * p2.p.q;
    pos[0] = p2.r.p[0];
    pos[1] = p2.r.p[1];
    pos[2] = 2 * elc_params.h - tp2;

    eng += p3m_pair_energy(q, get_mi_vector(p1.r.p, pos, box_geo).norm());
  }

  return (eng);
}

//////////////////////////////////////////////////////////////////////////////////

double ELC_P3M_dielectric_layers_energy_self(ParticleRange const &particles) {
  Utils::Vector3d pos;
  double q;
  double eng = 0.0;

  // Loop cell neighbors
  for (auto const &p : particles) {
    // Loop neighbor cell particles

    if (p.r.p[2] < elc_params.space_layer) {
      q = elc_params.delta_mid_bot * p.p.q * p.p.q;
      pos[0] = p.r.p[0];
      pos[1] = p.r.p[1];
      pos[2] = -p.r.p[2];

      eng += p3m_pair_energy(q, get_mi_vector(p.r.p, pos, box_geo).norm());
    }

    if (p.r.p[2] > (elc_params.h - elc_params.space_layer)) {
      q = elc_params.delta_mid_top * p.p.q * p.p.q;
      pos[0] = p.r.p[0];
      pos[1] = p.r.p[1];
      pos[2] = 2 * elc_params.h - p.r.p[2];

      eng += p3m_pair_energy(q, get_mi_vector(p.r.p, pos, box_geo).norm());
    }
  }
  return (eng);
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
