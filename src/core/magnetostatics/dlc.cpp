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

#ifdef DIPOLES

#include "magnetostatics/dlc.hpp"

#include "magnetostatics/dds.hpp"
#include "magnetostatics/dipoles.hpp"
#include "magnetostatics/dp3m.hpp"
#include "p3m/common.hpp"

#include "Particle.hpp"
#include "cells.hpp"
#include "communication.hpp"
#include "errorhandling.hpp"
#include "grid.hpp"
#include "particle_data.hpp"

#include <utils/constants.hpp>
#include <utils/math/sqr.hpp>

#include <boost/mpi/collectives/all_reduce.hpp>
#include <boost/mpi/collectives/reduce.hpp>
#include <boost/mpi/operations.hpp>

#include <mpi.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <functional>
#include <numeric>
#include <stdexcept>
#include <vector>

void DipolarLayerCorrection::check_gap(Particle const &p) const {
  if (p.dipm() != 0.) {
    auto const z = p.pos()[2];
    if (z < 0. or z > dlc.box_h) {
      runtimeErrorMsg() << "Particle " << p.id() << " entered DLC gap region "
                        << "by " << z - ((z < 0.) ? 0. : dlc.box_h);
    }
  }
}

/** Calculate the maximal dipole moment in the system */
static double calc_mu_max() {
  auto const local_particles = cell_structure.local_particles();
  auto const mu_max_local = std::accumulate(
      local_particles.begin(), local_particles.end(), 0.,
      [](double mu, Particle const &p) { return std::max(mu, p.dipm()); });

  return boost::mpi::all_reduce(comm_cart, mu_max_local,
                                boost::mpi::maximum<double>());
}

static double g1_DLC_dip(double g, double x) {
  auto const c = g / x;
  auto const cc2 = c * c;
  auto const x3 = x * x * x;
  auto const a = g * g * g / x + 1.5 * cc2 + 1.5 * g / x3 + 0.75 / (x3 * x);
  return a;
}

static double g2_DLC_dip(double g, double x) {
  auto const x2 = x * x;
  auto const a = g * g / x + 2.0 * g / x2 + 2.0 / (x2 * x);
  return a;
}

/** Compute the box magnetic dipole. */
static Utils::Vector3d calc_slab_dipole(ParticleRange const &particles) {

  Utils::Vector3d local_dip{};
  for (auto const &p : particles) {
    if (p.dipm() != 0.) {
      local_dip += p.calc_dip();
    }
  }

  return boost::mpi::all_reduce(comm_cart, local_dip, std::plus<>());
}

/**
 * @brief Compute the dipolar force and torque corrections.
 * %Algorithm implemented accordingly to @cite brodka04a.
 */
static void dipolar_force_corrections(int kcut,
                                      std::vector<Utils::Vector3d> &fs,
                                      std::vector<Utils::Vector3d> &ts,
                                      ParticleRange const &particles) {
  auto const facux = 2. * Utils::pi() * box_geo.length_inv()[0];
  auto const facuy = 2. * Utils::pi() * box_geo.length_inv()[1];

  auto const n_local_particles = particles.size();
  std::vector<double> ReSjp(n_local_particles);
  std::vector<double> ReSjm(n_local_particles);
  std::vector<double> ImSjp(n_local_particles);
  std::vector<double> ImSjm(n_local_particles);
  std::vector<double> ReGrad_Mup(n_local_particles);
  std::vector<double> ImGrad_Mup(n_local_particles);
  std::vector<double> ReGrad_Mum(n_local_particles);
  std::vector<double> ImGrad_Mum(n_local_particles);

  for (int ix = -kcut; ix <= +kcut; ix++) {
    for (int iy = -kcut; iy <= +kcut; iy++) {
      if (ix == 0 and iy == 0) {
        continue;
      }
      auto const gx = static_cast<double>(ix) * facux;
      auto const gy = static_cast<double>(iy) * facuy;
      auto const gr = sqrt(gx * gx + gy * gy);

      // We assume short slab direction is in the z-direction
      auto const fa1 = 1. / (gr * (exp(gr * box_geo.length()[2]) - 1.));

      // ... Compute S+,(S+)*,S-,(S-)*, and Spj,Smj for the current g

      std::size_t ip = 0;
      double S[4] = {0., 0., 0., 0.}; // S of Brodka method, or is S[4] =
                                      // {Re(S+), Im(S+), Re(S-), Im(S-)}
      for (auto const &p : particles) {
        if (p.dipm() != 0.) {
          auto const &pos = p.pos();
          auto const dip = p.calc_dip();

          auto const a = gx * dip[0] + gy * dip[1];
          auto const b = gr * dip[2];
          auto const er = gx * pos[0] + gy * pos[1];
          auto const ez = gr * pos[2];
          auto const c = cos(er);
          auto const d = sin(er);
          auto const f = exp(ez);

          ReSjp[ip] = (b * c - a * d) * f;
          ImSjp[ip] = (c * a + b * d) * f;
          ReSjm[ip] = (-b * c - a * d) / f;
          ImSjm[ip] = (c * a - b * d) / f;
          ReGrad_Mup[ip] = c * f;
          ReGrad_Mum[ip] = c / f;
          ImGrad_Mup[ip] = d * f;
          ImGrad_Mum[ip] = d / f;

          S[0] += ReSjp[ip];
          S[1] += ImSjp[ip];
          S[2] += ReSjm[ip];
          S[3] += ImSjm[ip];
        }
        ip++;
      }

      MPI_Allreduce(MPI_IN_PLACE, S, 4, MPI_DOUBLE, MPI_SUM, comm_cart);

      // ... Now we can compute the contributions to E,Fj,Ej for the current
      // g-value
      ip = 0;
      for (auto &p : particles) {
        if (p.dipm() != 0.) {
          {
            // compute contributions to the forces
            auto const s1 = -(-ReSjp[ip] * S[3] + ImSjp[ip] * S[2]);
            auto const s2 = +(ReSjm[ip] * S[1] - ImSjm[ip] * S[0]);
            auto const s3 = -(-ReSjm[ip] * S[1] + ImSjm[ip] * S[0]);
            auto const s4 = +(ReSjp[ip] * S[3] - ImSjp[ip] * S[2]);

            auto const s1z = +(ReSjp[ip] * S[2] + ImSjp[ip] * S[3]);
            auto const s2z = -(ReSjm[ip] * S[0] + ImSjm[ip] * S[1]);
            auto const s3z = -(ReSjm[ip] * S[0] + ImSjm[ip] * S[1]);
            auto const s4z = +(ReSjp[ip] * S[2] + ImSjp[ip] * S[3]);

            auto const ss = s1 + s2 + s3 + s4;
            fs[ip][0] += fa1 * gx * ss;
            fs[ip][1] += fa1 * gy * ss;
            fs[ip][2] += fa1 * gr * (s1z + s2z + s3z + s4z);
          }
          {
            // compute contributions to the electrical field
            auto const s1 = -(-ReGrad_Mup[ip] * S[3] + ImGrad_Mup[ip] * S[2]);
            auto const s2 = +(ReGrad_Mum[ip] * S[1] - ImGrad_Mum[ip] * S[0]);
            auto const s3 = -(-ReGrad_Mum[ip] * S[1] + ImGrad_Mum[ip] * S[0]);
            auto const s4 = +(ReGrad_Mup[ip] * S[3] - ImGrad_Mup[ip] * S[2]);

            auto const s1z = +(ReGrad_Mup[ip] * S[2] + ImGrad_Mup[ip] * S[3]);
            auto const s2z = -(ReGrad_Mum[ip] * S[0] + ImGrad_Mum[ip] * S[1]);
            auto const s3z = -(ReGrad_Mum[ip] * S[0] + ImGrad_Mum[ip] * S[1]);
            auto const s4z = +(ReGrad_Mup[ip] * S[2] + ImGrad_Mup[ip] * S[3]);

            auto const ss = s1 + s2 + s3 + s4;
            ts[ip][0] += fa1 * gx * ss;
            ts[ip][1] += fa1 * gy * ss;
            ts[ip][2] += fa1 * gr * (s1z + s2z + s3z + s4z);
          }
        }
        ++ip;
      }
    }
  }

  // Convert from the corrections to the electrical field to the corrections
  // for the torques
  std::size_t ip = 0;
  for (auto const &p : particles) {
    if (p.dipm() != 0.) {
      ts[ip] = vector_product(p.calc_dip(), ts[ip]);
    }
    ++ip;
  }

  // Multiply by the factors we have left during the loops

  auto const piarea =
      Utils::pi() * box_geo.length_inv()[0] * box_geo.length_inv()[1];
  for (std::size_t j = 0; j < n_local_particles; ++j) {
    fs[j] *= piarea;
    ts[j] *= piarea;
  }
}

/**
 * @brief Compute the dipolar DLC energy correction.
 * %Algorithm implemented accordingly to @cite brodka04a.
 */
static double dipolar_energy_correction(int kcut,
                                        ParticleRange const &particles) {
  auto const facux = 2. * Utils::pi() * box_geo.length_inv()[0];
  auto const facuy = 2. * Utils::pi() * box_geo.length_inv()[1];

  double energy = 0.;
  double sum_S[4] = {0., 0., 0., 0.};
  for (int ix = -kcut; ix <= +kcut; ix++) {
    for (int iy = -kcut; iy <= +kcut; iy++) {
      if ((ix == 0 && iy == 0)) {
        continue;
      }
      auto const gx = static_cast<double>(ix) * facux;
      auto const gy = static_cast<double>(iy) * facuy;
      auto const gr = sqrt(gx * gx + gy * gy);

      // We assume short slab direction is in the z-direction
      auto const fa1 = 1. / (gr * (exp(gr * box_geo.length()[2]) - 1.));

      // ... Compute S+,(S+)*,S-,(S-)*, and Spj,Smj for the current g

      double S[4] = {0., 0., 0., 0.}; // S of Brodka method, or is S[4] =
                                      // {Re(S+), Im(S+), Re(S-), Im(S-)}
      for (auto const &p : particles) {
        if (p.dipm() != 0.) {
          auto const &pos = p.pos();
          auto const dip = p.calc_dip();

          auto const a = gx * dip[0] + gy * dip[1];
          auto const b = gr * dip[2];
          auto const er = gx * pos[0] + gy * pos[1];
          auto const ez = gr * pos[2];
          auto const c = cos(er);
          auto const d = sin(er);
          auto const f = exp(ez);

          S[0] += (b * c - a * d) * f;
          S[1] += (c * a + b * d) * f;
          S[2] += (-b * c - a * d) / f;
          S[3] += (c * a - b * d) / f;
        }
      }
      boost::mpi::reduce(comm_cart, S, 4, sum_S, std::plus<>(), 0);

      // compute contribution to the energy
      // s2=(ReSm*ReSp+ImSm*ImSp); s2=s1!!!
      energy += fa1 * 2. * (sum_S[0] * sum_S[2] + sum_S[1] * sum_S[3]);
    }
  }

  auto const piarea =
      Utils::pi() * box_geo.length_inv()[0] * box_geo.length_inv()[1];
  energy *= -piarea;
  return (this_node == 0) ? energy : 0.;
}

void DipolarLayerCorrection::add_force_corrections(
    ParticleRange const &particles) const {
  assert(dlc.far_cut > 0.);
  auto const volume = box_geo.volume();
  auto const correc = 4. * Utils::pi() / volume;

  // --- Create arrays that should contain the corrections to
  //     the forces and torques, and set them to zero.
  std::vector<Utils::Vector3d> dip_DLC_f(particles.size());
  std::vector<Utils::Vector3d> dip_DLC_t(particles.size());

  //---- Compute the corrections ----------------------------------

  // First the DLC correction
  auto const kcut = static_cast<int>(std::round(dlc.far_cut));
  dipolar_force_corrections(kcut, dip_DLC_f, dip_DLC_t, particles);

  // Now we compute the correction like Yeh and Klapp to take into account
  // the fact that you are using a 3D PBC method which uses spherical
  // summation instead of slab-wise summation.
  // Slab-wise summation is the one required to apply DLC correction.
  // This correction is often called SDC = Shape Dependent Correction.
  // See @cite brodka04a.

  auto const box_dip = calc_slab_dipole(particles);

  // --- Transfer the computed corrections to the Forces, Energy and torques
  //     of the particles

  std::size_t ip = 0;
  for (auto &p : particles) {
    check_gap(p);

    if (p.dipm() != 0.) {
      // SDC correction term is zero for the forces
      p.force() += prefactor * dip_DLC_f[ip];

      auto const dip = p.calc_dip();
      // SDC correction for the torques
      auto d = Utils::Vector3d{0., 0., -correc * box_dip[2]};
#ifdef DP3M
      if (epsilon != P3M_EPSILON_METALLIC) {
        d += correc * epsilon_correction * box_dip;
      }
#endif
      p.torque() += prefactor * (dip_DLC_t[ip] + vector_product(dip, d));
    }
    ++ip;
  }
}

double DipolarLayerCorrection::energy_correction(
    ParticleRange const &particles) const {
  assert(dlc.far_cut > 0.);
  auto const volume = box_geo.volume();
  auto const pref = prefactor * 2. * Utils::pi() / volume;

  // Check if particles aren't in the forbidden gap region
  // This loop is needed, because there is no other guaranteed
  // single pass over all particles in this function.
  for (auto const &p : particles) {
    check_gap(p);
  }

  //---- Compute the corrections ----------------------------------

  // First the DLC correction
  auto const k_cut = static_cast<int>(std::round(dlc.far_cut));
  auto dip_DLC_energy = prefactor * dipolar_energy_correction(k_cut, particles);

  // Now we compute the correction like Yeh and Klapp to take into account
  // the fact that you are using a 3D PBC method which uses spherical
  // summation instead of slab-wise summation.
  // Slab-wise summation is the one required to apply DLC correction.
  // This correction is often called SDC = Shape Dependent Correction.
  // See @cite brodka04a.

  auto const box_dip = calc_slab_dipole(particles);

  if (this_node == 0) {
    dip_DLC_energy += pref * Utils::sqr(box_dip[2]);
#ifdef DP3M
    if (epsilon != P3M_EPSILON_METALLIC) {
      dip_DLC_energy -= pref * epsilon_correction * box_dip.norm2();
    }
#endif
    return dip_DLC_energy;
  }
  return 0.;
}

static int count_magnetic_particles() {
  int local_n = 0;

  for (auto const &p : cell_structure.local_particles()) {
    if (p.dipm() != 0.) {
      local_n++;
    }
  }

  return boost::mpi::all_reduce(comm_cart, local_n, std::plus<>());
}

/** Compute the cut-off in the DLC dipolar part to get a certain accuracy.
 *  We assume particles to have all the same value of the dipolar momentum
 *  modulus, which is taken as the largest value of mu inside the system.
 *  If we assume the gap has a width @c gap_size (within
 *  which there is no particles): <tt>Lz = h + gap_size</tt>
 */
double DipolarLayerCorrection::tune_far_cut() const {
  /* we take the maximum dipole in the system, to be sure that the errors
   * in the other case will be equal or less than for this one */
  auto const mu_max_sq = Utils::sqr(calc_mu_max());
  auto const lx = box_geo.length()[0];
  auto const ly = box_geo.length()[1];
  auto const lz = box_geo.length()[2];

  if (std::abs(lx - ly) > 0.001) {
    throw std::runtime_error("DLC tuning: box size in x direction is "
                             "different from y direction. The tuning "
                             "formula requires both to be equal.");
  }

  auto constexpr limitkc = 200;
  auto const piarea = Utils::pi() / (lx * ly);
  auto const nmp = static_cast<double>(count_magnetic_particles());
  auto const h = dlc.box_h;
  auto far_cut = -1.;
  for (int kc = 1; kc < limitkc; kc++) {
    auto const gc = kc * 2. * Utils::pi() / lx;
    auto const fa0 = sqrt(9. * exp(+2. * gc * h) * g1_DLC_dip(gc, lz - h) +
                          9. * exp(-2. * gc * h) * g1_DLC_dip(gc, lz + h) +
                          22. * g1_DLC_dip(gc, lz));
    auto const fa1 = sqrt(0.125 * piarea) * fa0;
    auto const fa2 = g2_DLC_dip(gc, lz);
    auto const de = nmp * mu_max_sq / (4. * (exp(gc * lz) - 1.)) * (fa1 + fa2);
    if (de < dlc.maxPWerror) {
      far_cut = static_cast<double>(kc);
      break;
    }
  }
  if (far_cut <= 0.) {
    throw std::runtime_error("DLC tuning failed: maxPWerror too small");
  }
  return far_cut;
}

/** @brief Lock an actor and modify its parameters. */
struct AdaptSolver : public boost::static_visitor<void> {
  explicit AdaptSolver(DipolarLayerCorrection *this_ptr) : m_actor{this_ptr} {}

  void operator()(std::shared_ptr<DipolarDirectSum> const &solver) {
    m_actor->prefactor = solver->prefactor;
    m_actor->epsilon = P3M_EPSILON_METALLIC;
  }

#ifdef DP3M
  void operator()(std::shared_ptr<DipolarP3M> const &solver) {
    m_actor->prefactor = solver->prefactor;
    m_actor->epsilon = solver->dp3m.params.epsilon;
  }
#endif

private:
  DipolarLayerCorrection *m_actor;
};

void DipolarLayerCorrection::adapt_solver() {
  auto visitor = AdaptSolver{this};
  boost::apply_visitor(visitor, base_solver);
  epsilon_correction =
      (epsilon == P3M_EPSILON_METALLIC) ? 0. : 1. / (2. * epsilon + 1.);
}

void DipolarLayerCorrection::sanity_checks_node_grid() const {
  if (!box_geo.periodic(0) || !box_geo.periodic(1) || !box_geo.periodic(2)) {
    throw std::runtime_error("DLC: requires periodicity (True, True, True)");
  }
}

dlc_data::dlc_data(double maxPWerror, double gap_size, double far_cut)
    : maxPWerror{maxPWerror}, gap_size{gap_size}, box_h{box_geo.length()[2] -
                                                        gap_size},
      far_cut{far_cut}, far_calculated{far_cut == -1.} {
  if (far_cut <= 0. and not far_calculated) {
    throw std::domain_error("Parameter 'far_cut' must be > 0");
  }
  if (maxPWerror <= 0.) {
    throw std::domain_error("Parameter 'maxPWerror' must be > 0");
  }
  if (gap_size <= 0.) {
    throw std::domain_error("Parameter 'gap_size' must be > 0");
  }
}

void DipolarLayerCorrection::recalc_box_h() {
  auto const new_box_h = box_geo.length()[2] - dlc.gap_size;
  if (new_box_h < 0.) {
    throw std::runtime_error("DLC gap size (" + std::to_string(dlc.gap_size) +
                             ") larger than box length in z-direction (" +
                             std::to_string(box_geo.length()[2]) + ")");
  }
  dlc.box_h = new_box_h;
}

DipolarLayerCorrection::DipolarLayerCorrection(dlc_data &&parameters,
                                               BaseSolver &&solver)
    : dlc{parameters}, base_solver{solver} {
  adapt_solver();
}

#endif // DIPOLES
