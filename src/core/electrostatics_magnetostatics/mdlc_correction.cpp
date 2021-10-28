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
 */

#include "electrostatics_magnetostatics/mdlc_correction.hpp"

#ifdef DIPOLES
#include "Particle.hpp"
#include "cells.hpp"
#include "communication.hpp"
#include "electrostatics_magnetostatics/common.hpp"
#include "electrostatics_magnetostatics/dipole.hpp"
#include "electrostatics_magnetostatics/p3m-dipolar.hpp"
#include "errorhandling.hpp"
#include "grid.hpp"
#include "particle_data.hpp"

#include <utils/constants.hpp>
#include <utils/math/sqr.hpp>

#include <boost/mpi.hpp>
#include <boost/mpi/collectives/all_reduce.hpp>
#include <boost/mpi/operations.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <functional>
#include <numeric>
#include <stdexcept>
#include <vector>

DLC_struct dlc_params = {1e100, 0., 0., false, 0.};

/** Checks if a magnetic particle is in the forbidden gap region
 */
inline void check_gap_mdlc(const Particle &p) {
  if (p.p.dipm != 0.0) {
    if (p.r.p[2] < 0.0)
      runtimeErrorMsg() << "Particle " << p.p.identity << " entered MDLC gap "
                        << "region by " << (p.r.p[2]);
    else if (p.r.p[2] > dlc_params.h) {
      runtimeErrorMsg() << "Particle " << p.p.identity << " entered MDLC gap "
                        << "region by " << (p.r.p[2] - dlc_params.h);
    }
  }
}

/** Calculate the maximal dipole moment in the system */
double calc_mu_max() {
  auto const local_particles = cell_structure.local_particles();
  auto const mu_max_local = std::accumulate(
      local_particles.begin(), local_particles.end(), 0.0,
      [](double mu, Particle const &p) { return std::max(mu, p.p.dipm); });

  double mu_max;
  boost::mpi::reduce(comm_cart, mu_max_local, mu_max,
                     boost::mpi::maximum<double>(), 0);
  return mu_max;
}

REGISTER_CALLBACK_MAIN_RANK(calc_mu_max)

inline double g1_DLC_dip(double g, double x) {
  auto const c = g / x;
  auto const cc2 = c * c;
  auto const x3 = x * x * x;
  auto const a = g * g * g / x + 1.5 * cc2 + 1.5 * g / x3 + 0.75 / (x3 * x);
  return a;
}

inline double g2_DLC_dip(double g, double x) {
  auto const x2 = x * x;
  auto const a = g * g / x + 2.0 * g / x2 + 2.0 / (x2 * x);
  return a;
}

/** Compute the box magnetic dipole. */
inline Utils::Vector3d calc_slab_dipole(const ParticleRange &particles) {

  Utils::Vector3d local_dip{};
  for (auto const &p : particles) {
    if (p.p.dipm != 0.0) {
      local_dip += p.calc_dip();
    }
  }

  return boost::mpi::all_reduce(comm_cart, local_dip, std::plus<>());
}

/** Compute the dipolar DLC corrections for forces and torques.
 *  %Algorithm implemented accordingly to @cite brodka04a.
 */
double get_DLC_dipolar(int kcut, std::vector<Utils::Vector3d> &fs,
                       std::vector<Utils::Vector3d> &ts,
                       const ParticleRange &particles) {
  auto const n_local_particles = particles.size();

  std::vector<double> ReSjp(n_local_particles), ReSjm(n_local_particles);
  std::vector<double> ImSjp(n_local_particles), ImSjm(n_local_particles);
  std::vector<double> ReGrad_Mup(n_local_particles),
      ImGrad_Mup(n_local_particles);
  std::vector<double> ReGrad_Mum(n_local_particles),
      ImGrad_Mum(n_local_particles);
  double s1, s2, s3, s4;
  double s1z, s2z, s3z, s4z;
  double ss;

  auto const facux = 2.0 * Utils::pi() * box_geo.length_inv()[0];
  auto const facuy = 2.0 * Utils::pi() * box_geo.length_inv()[1];
  double energy = 0.0;

  for (int ix = -kcut; ix <= +kcut; ix++) {
    for (int iy = -kcut; iy <= +kcut; iy++) {
      if (!(ix == 0 && iy == 0)) {
        auto const gx = static_cast<double>(ix) * facux;
        auto const gy = static_cast<double>(iy) * facuy;

        auto const gr = sqrt(gx * gx + gy * gy);

        // We assume short slab direction is z direction
        auto const fa1 = 1. / (gr * (exp(gr * box_geo.length()[2]) - 1.0));

        // ... Compute S+,(S+)*,S-,(S-)*, and Spj,Smj for the current g

        double S[4] = {0.0, 0.0, 0.0, 0.0}; // S of Brodka method, or is S[4] =
                                            // {Re(S+), Im(S+), Re(S-), Im(S-)}
        int ip = 0;

        for (auto const &p : particles) {
          if (p.p.dipm > 0) {
            Utils::Vector3d const dip = p.calc_dip();

            auto const a = gx * dip[0] + gy * dip[1];
            auto const b = gr * dip[2];
            auto const er = gx * p.r.p[0] + gy * p.r.p[1];
            auto const ez = gr * p.r.p[2];
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

        // We compute the contribution to the energy ............

        // s2=(ReSm*ReSp+ImSm*ImSp); s2=s1!!!

        energy += fa1 * ((S[0] * S[2] + S[1] * S[3]) * 2.0);

        // ... Now we can compute the contributions to E,Fj,Ej for the current
        // g-value
        ip = 0;
        for (auto &p : particles) {
          if (p.p.dipm > 0) {
            // We compute the contributions to the forces ............

            s1 = -(-ReSjp[ip] * S[3] + ImSjp[ip] * S[2]);
            s2 = +(ReSjm[ip] * S[1] - ImSjm[ip] * S[0]);
            s3 = -(-ReSjm[ip] * S[1] + ImSjm[ip] * S[0]);
            s4 = +(ReSjp[ip] * S[3] - ImSjp[ip] * S[2]);

            s1z = +(ReSjp[ip] * S[2] + ImSjp[ip] * S[3]);
            s2z = -(ReSjm[ip] * S[0] + ImSjm[ip] * S[1]);
            s3z = -(ReSjm[ip] * S[0] + ImSjm[ip] * S[1]);
            s4z = +(ReSjp[ip] * S[2] + ImSjp[ip] * S[3]);

            ss = s1 + s2 + s3 + s4;
            fs[ip][0] += fa1 * gx * ss;
            fs[ip][1] += fa1 * gy * ss;
            fs[ip][2] += fa1 * gr * (s1z + s2z + s3z + s4z);

            // We compute the contributions to the electrical field
            // ............

            s1 = -(-ReGrad_Mup[ip] * S[3] + ImGrad_Mup[ip] * S[2]);
            s2 = +(ReGrad_Mum[ip] * S[1] - ImGrad_Mum[ip] * S[0]);
            s3 = -(-ReGrad_Mum[ip] * S[1] + ImGrad_Mum[ip] * S[0]);
            s4 = +(ReGrad_Mup[ip] * S[3] - ImGrad_Mup[ip] * S[2]);

            s1z = +(ReGrad_Mup[ip] * S[2] + ImGrad_Mup[ip] * S[3]);
            s2z = -(ReGrad_Mum[ip] * S[0] + ImGrad_Mum[ip] * S[1]);
            s3z = -(ReGrad_Mum[ip] * S[0] + ImGrad_Mum[ip] * S[1]);
            s4z = +(ReGrad_Mup[ip] * S[2] + ImGrad_Mup[ip] * S[3]);

            ss = s1 + s2 + s3 + s4;
            ts[ip][0] += fa1 * gx * ss;
            ts[ip][1] += fa1 * gy * ss;
            ts[ip][2] += fa1 * gr * (s1z + s2z + s3z + s4z);
          } // if dipm>0 ....
          ip++;
        } // loop j
      }   // end of if(ii> ...
    }
  } // end of loops for gx,gy

  // Convert from the corrections to the Electrical field to the corrections
  // for the torques ....

  int ip = 0;
  for (auto const &p : particles) {
    if (p.p.dipm > 0) {
      ts[ip] = vector_product(p.calc_dip(), ts[ip]);
    }
    ip++;
  }

  // Multiply by the factors we have left during the loops

  auto const piarea =
      Utils::pi() * box_geo.length_inv()[0] * box_geo.length_inv()[1];

  for (int j = 0; j < n_local_particles; j++) {
    fs[j] *= piarea;
    ts[j] *= piarea;
  }

  energy *= (-piarea);

  return energy;
}

/** Compute the dipolar DLC corrections
 *  %Algorithm implemented accordingly to @cite brodka04a.
 */
double get_DLC_energy_dipolar(int kcut, const ParticleRange &particles) {
  auto const facux = 2.0 * Utils::pi() * box_geo.length_inv()[0];
  auto const facuy = 2.0 * Utils::pi() * box_geo.length_inv()[1];

  double energy = 0.0;
  for (int ix = -kcut; ix <= +kcut; ix++) {
    for (int iy = -kcut; iy <= +kcut; iy++) {

      if (!(ix == 0 && iy == 0)) {
        auto const gx = static_cast<double>(ix) * facux;
        auto const gy = static_cast<double>(iy) * facuy;
        auto const gr = sqrt(gx * gx + gy * gy);
        // We assume short slab direction is z direction
        auto const fa1 = 1. / (gr * (exp(gr * box_geo.length()[2]) - 1.0));

        // ... Compute S+,(S+)*,S-,(S-)*, and Spj,Smj for the current g

        double S[4] = {0.0, 0.0, 0.0, 0.0};
        int ip = 0;
        for (auto const &p : particles) {
          if (p.p.dipm > 0) {
            const Utils::Vector3d dip = p.calc_dip();

            auto const a = gx * dip[0] + gy * dip[1];
            {
              auto const b = gr * dip[2];
              auto const er = gx * p.r.p[0] + gy * p.r.p[1];
              auto const ez = gr * p.r.p[2];
              auto const c = cos(er);
              auto const d = sin(er);
              auto const f = exp(ez);

              S[0] += (b * c - a * d) * f;
              S[1] += (c * a + b * d) * f;
              S[2] += (-b * c - a * d) / f;
              S[3] += (c * a - b * d) / f;
            }
          }
          ip++;
        }

        double global_S[4];
        MPI_Reduce(S, global_S, 4, MPI_DOUBLE, MPI_SUM, 0, comm_cart);

        // We compute the contribution to the energy ............
        auto const s1 = global_S[0] * global_S[2] + global_S[1] * global_S[3];
        // s2=(ReSm*ReSp+ImSm*ImSp); s2=s1!!!

        energy += fa1 * (s1 * 2.0);

      } // end of if(...
    }
  } // end of loops for gx,gy

  // Multiply by the factors we have left during the loops

  auto const piarea =
      Utils::pi() * box_geo.length_inv()[0] * box_geo.length_inv()[1];
  energy *= (-piarea);
  return (this_node == 0) ? energy : 0.0;
}

/** Compute and add the terms needed to correct the 3D dipolar
 *  methods when we have a slab geometry
 */
void add_mdlc_force_corrections(const ParticleRange &particles) {
  auto const volume = box_geo.volume();
  auto const correc = 4. * Utils::pi() / volume;

  // --- Create arrays that should contain the corrections to
  //     the forces and torques, and set them to zero.
  std::vector<Utils::Vector3d> dip_DLC_f(particles.size());
  std::vector<Utils::Vector3d> dip_DLC_t(particles.size());

  //---- Compute the corrections ----------------------------------

  // First the DLC correction
  get_DLC_dipolar(static_cast<int>(std::round(dlc_params.far_cut)), dip_DLC_f,
                  dip_DLC_t, particles);

  // Now we compute the correction like Yeh and Klapp to take into account
  // the fact that you are using a 3D PBC method which uses spherical
  // summation instead of slab-wise summation.
  // Slab-wise summation is the one required to apply DLC correction.
  // This correction is often called SDC = Shape Dependent Correction.
  // See @cite brodka04a.

  auto const box_dip = calc_slab_dipole(particles);

  // --- Transfer the computed corrections to the Forces, Energy and torques
  //     of the particles

  int ip = 0;
  for (auto &p : particles) {
    check_gap_mdlc(p);

    if (p.p.dipm != 0.0) {
      // SDC correction term is zero for the forces
      p.f.f += dipole.prefactor * dip_DLC_f[ip];

      auto const dip = p.calc_dip();
      // SDC correction for the torques
      Utils::Vector3d d = {0.0, 0.0, -correc * box_dip[2]};
#ifdef DP3M
      if (dipole.method == DIPOLAR_MDLC_P3M and
          dp3m.params.epsilon != P3M_EPSILON_METALLIC) {
        auto const correps = correc / (2.0 * dp3m.params.epsilon + 1.0);
        d += correps * box_dip;
      }
#endif
      p.f.torque += dipole.prefactor * (dip_DLC_t[ip] + vector_product(dip, d));
    }
    ip++;
  }
}

/** Compute and add the terms needed to correct the energy of
 *  3D dipolar methods when we have a slab geometry
 */
double add_mdlc_energy_corrections(const ParticleRange &particles) {

  auto const volume = box_geo.volume();
  auto const prefactor = dipole.prefactor * 2. * Utils::pi() / volume;

  // Check if particles aren't in the forbidden gap region
  // This loop is needed, because there is no other guaranteed
  // single pass over all particles in this function.
  for (auto const &p : particles) {
    check_gap_mdlc(p);
  }

  //---- Compute the corrections ----------------------------------

  // First the DLC correction
  auto const k_cut = static_cast<int>(std::round(dlc_params.far_cut));
  double dip_DLC_energy =
      dipole.prefactor * get_DLC_energy_dipolar(k_cut, particles);

  // Now we compute the correction like Yeh and Klapp to take into account
  // the fact that you are using a 3D PBC method which uses spherical
  // summation instead of slab-wise summation.
  // Slab-wise summation is the one required to apply DLC correction.
  // This correction is often called SDC = Shape Dependent Correction.
  // See @cite brodka04a.

  auto const box_dip = calc_slab_dipole(particles);

  if (this_node == 0) {
    dip_DLC_energy += prefactor * Utils::sqr(box_dip[2]);
#ifdef DP3M
    if (dipole.method == DIPOLAR_MDLC_P3M and
        dp3m.params.epsilon != P3M_EPSILON_METALLIC) {
      auto const correps = 1.0 / (2.0 * dp3m.params.epsilon + 1.0);
      dip_DLC_energy -= prefactor * box_dip.norm2() * correps;
    }
#endif
    return dip_DLC_energy;
  }
  return 0.0;
}

/** Compute the cut-off in the DLC dipolar part to get a certain accuracy.
 *  We assume particles to have all the same value of the dipolar momentum
 *  modulus, which is taken as the largest value of mu  inside the system.
 *  If we assume the gap has a width @c gap_size (within
 *  which there is no particles): <tt>Lz = h + gap_size</tt>
 *
 *  BE CAREFUL:
 *  1. We assume the short distance for the slab to be in the *z*-direction
 *  2. You must also tune the other 3D method to the same accuracy, otherwise
 *     it makes no sense to have an accurate result for DLC-dipolar.
 */
double mdlc_tune_far_cut(DLC_struct const &params) {
  /* we take the maximum dipole in the system, to be sure that the errors
   * in the other case will be equal or less than for this one */
  auto const mu_max = mpi_call(Communication::Result::main_rank, calc_mu_max);
  auto const mu_max_sq = mu_max * mu_max;

  const double n = get_n_part();
  auto const lx = box_geo.length()[0];
  auto const ly = box_geo.length()[1];
  auto const lz = box_geo.length()[2];
  auto const h = params.h;

  if (std::abs(lx - ly) > 0.001) {
    throw std::runtime_error("MDLC tuning: box size in x direction is "
                             "different from y direction. The tuning "
                             "formula requires both to be equal.");
  }

  constexpr int limitkc = 200;
  for (int kc = 1; kc < limitkc; kc++) {
    auto const gc = kc * 2.0 * Utils::pi() / lx;
    auto const fa0 = sqrt(9.0 * exp(+2.0 * gc * h) * g1_DLC_dip(gc, lz - h) +
                          9.0 * exp(-2.0 * gc * h) * g1_DLC_dip(gc, lz + h) +
                          22.0 * g1_DLC_dip(gc, lz));
    auto const fa1 = 0.5 * sqrt(Utils::pi() / (2.0 * lx * ly)) * fa0;
    auto const fa2 = g2_DLC_dip(gc, lz);
    auto const de = n * mu_max_sq / (4.0 * (exp(gc * lz) - 1.0)) * (fa1 + fa2);
    if (de < params.maxPWerror) {
      return static_cast<double>(kc);
    }
  }

  throw std::runtime_error("MDLC tuning failed: unable to find a proper "
                           "cut-off for the given accuracy.");
}

void mdlc_sanity_checks() {
  if (!box_geo.periodic(0) || !box_geo.periodic(1) || !box_geo.periodic(2)) {
    throw std::runtime_error("MDLC requires periodicity 1 1 1");
  }
}

void mdlc_set_params(double maxPWerror, double gap_size, double far_cut) {
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
  mdlc_sanity_checks();

  DipolarInteraction new_method;
  switch (dipole.method) {
#ifdef DP3M
  case DIPOLAR_MDLC_P3M:
  case DIPOLAR_P3M:
    new_method = DIPOLAR_MDLC_P3M;
    break;
#endif
  case DIPOLAR_MDLC_DS:
  case DIPOLAR_DS:
    new_method = DIPOLAR_MDLC_DS;
    fprintf(stderr, "You are not using the P3M method, therefore dp3m.params."
                    "epsilon unknown, I will assume metallic borders.\n");
    break;
  default:
    throw std::runtime_error(
        "MDLC cannot extend the currently active magnetostatics solver.");
  }

  DLC_struct new_dlc_params{maxPWerror, far_cut, gap_size, far_cut == -1., h};

  if (new_dlc_params.far_calculated) {
    new_dlc_params.far_cut = mdlc_tune_far_cut(new_dlc_params);
  }

  dlc_params = new_dlc_params;

  Dipole::set_method_local(new_method);
  mpi_bcast_coulomb_params();
}

#endif // DIPOLES
