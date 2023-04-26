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

/** @file
 *  P3M algorithm for long-range magnetic dipole-dipole interaction.
 *
 *  By default the magnetic epsilon is metallic = 0.
 */

#include "config/config.hpp"

#ifdef DP3M

#include "magnetostatics/dp3m.hpp"

#include "p3m/TuningAlgorithm.hpp"
#include "p3m/TuningLogger.hpp"
#include "p3m/common.hpp"
#include "p3m/fft.hpp"
#include "p3m/influence_function_dipolar.hpp"
#include "p3m/interpolation.hpp"
#include "p3m/send_mesh.hpp"

#include "Particle.hpp"
#include "ParticleRange.hpp"
#include "cell_system/CellStructureType.hpp"
#include "cells.hpp"
#include "communication.hpp"
#include "errorhandling.hpp"
#include "event.hpp"
#include "grid.hpp"
#include "integrate.hpp"
#include "npt.hpp"
#include "tuning.hpp"

#include <utils/Vector.hpp>
#include <utils/constants.hpp>
#include <utils/integral_parameter.hpp>
#include <utils/math/int_pow.hpp>
#include <utils/math/sinc.hpp>
#include <utils/math/sqr.hpp>

#include <boost/mpi/collectives/all_reduce.hpp>
#include <boost/mpi/collectives/reduce.hpp>
#include <boost/optional.hpp>

#include <algorithm>
#include <array>
#include <functional>
#include <sstream>
#include <stdexcept>
#include <vector>

void DipolarP3M::count_magnetic_particles() {
  int local_n = 0;
  double local_mu2 = 0.;

  for (auto const &p : cell_structure.local_particles()) {
    if (p.dipm() != 0.) {
      local_mu2 += p.calc_dip().norm2();
      local_n++;
    }
  }

  boost::mpi::all_reduce(comm_cart, local_mu2, dp3m.sum_mu2, std::plus<>());
  boost::mpi::all_reduce(comm_cart, local_n, dp3m.sum_dip_part, std::plus<>());
}

static double dp3m_k_space_error(double box_size, int mesh, int cao,
                                 int n_c_part, double sum_q2, double alpha_L);

static double dp3m_real_space_error(double box_size, double r_cut_iL,
                                    int n_c_part, double sum_q2,
                                    double alpha_L);
static void dp3m_tune_aliasing_sums(int nx, int ny, int nz, int mesh,
                                    double mesh_i, int cao, double alpha_L_i,
                                    double *alias1, double *alias2);

/** Compute the value of alpha through a bisection method.
 *  Based on eq. (33) @cite wang01a.
 */
double dp3m_rtbisection(double box_size, double r_cut_iL, int n_c_part,
                        double sum_q2, double x1, double x2, double xacc,
                        double tuned_accuracy);

double DipolarP3M::calc_average_self_energy_k_space() const {
  auto const start = Utils::Vector3i{dp3m.fft.plan[3].start};
  auto const size = Utils::Vector3i{dp3m.fft.plan[3].new_mesh};

  auto const node_phi = grid_influence_function_self_energy(
      dp3m.params, start, start + size, dp3m.g_energy);

  double phi = 0.;
  boost::mpi::reduce(comm_cart, node_phi, phi, std::plus<>(), 0);
  phi /= 3. * box_geo.length()[0] * Utils::int_pow<3>(dp3m.params.mesh[0]);
  return phi * Utils::pi();
}

void DipolarP3M::init() {
  assert(dp3m.params.mesh >= Utils::Vector3i::broadcast(1));
  assert(dp3m.params.cao >= 1 and dp3m.params.cao <= 7);
  assert(dp3m.params.alpha > 0.);

  dp3m.params.cao3 = Utils::int_pow<3>(dp3m.params.cao);
  dp3m.params.recalc_a_ai_cao_cut(box_geo.length());
  dp3m.local_mesh.calc_local_ca_mesh(dp3m.params, local_geo, skin, 0.);

  dp3m.sm.resize(comm_cart, dp3m.local_mesh);

  int ca_mesh_size = fft_init(dp3m.local_mesh.dim, dp3m.local_mesh.margin,
                              dp3m.params.mesh, dp3m.params.mesh_off,
                              dp3m.ks_pnum, dp3m.fft, node_grid, comm_cart);
  dp3m.rs_mesh.resize(ca_mesh_size);
  dp3m.ks_mesh.resize(ca_mesh_size);

  for (auto &val : dp3m.rs_mesh_dip) {
    val.resize(ca_mesh_size);
  }

  dp3m.calc_differential_operator();

  /* fix box length dependent constants */
  scaleby_box_l();

  count_magnetic_particles();
}

DipolarP3M::DipolarP3M(P3MParameters &&parameters, double prefactor,
                       int tune_timings, bool tune_verbose)
    : dp3m{std::move(parameters)}, prefactor{prefactor},
      tune_timings{tune_timings}, tune_verbose{tune_verbose} {

  m_is_tuned = !dp3m.params.tuning;
  dp3m.params.tuning = false;

  if (prefactor <= 0.) {
    throw std::domain_error("Parameter 'prefactor' must be > 0");
  }
  if (tune_timings <= 0) {
    throw std::domain_error("Parameter 'timings' must be > 0");
  }

  if (dp3m.params.mesh != Utils::Vector3i::broadcast(dp3m.params.mesh[0])) {
    throw std::domain_error("DipolarP3M requires a cubic mesh");
  }
}

namespace {
template <int cao> struct AssignDipole {
  void operator()(dp3m_data_struct &dp3m, Utils::Vector3d const &real_pos,
                  Utils::Vector3d const &dip) const {
    auto const weights = p3m_calculate_interpolation_weights<cao>(
        real_pos, dp3m.params.ai, dp3m.local_mesh);
    p3m_interpolate<cao>(dp3m.local_mesh, weights,
                         [&dip, &dp3m](int ind, double w) {
                           dp3m.rs_mesh_dip[0][ind] += w * dip[0];
                           dp3m.rs_mesh_dip[1][ind] += w * dip[1];
                           dp3m.rs_mesh_dip[2][ind] += w * dip[2];
                         });

    dp3m.inter_weights.store<cao>(weights);
  }
};
} // namespace

void DipolarP3M::dipole_assign(ParticleRange const &particles) {
  dp3m.inter_weights.reset(dp3m.params.cao);

  /* prepare local FFT mesh */
  for (auto &i : dp3m.rs_mesh_dip)
    for (int j = 0; j < dp3m.local_mesh.size; j++)
      i[j] = 0.;

  for (auto const &p : particles) {
    if (p.dipm() != 0.) {
      Utils::integral_parameter<int, AssignDipole, 1, 7>(dp3m.params.cao, dp3m,
                                                         p.pos(), p.calc_dip());
    }
  }
}

namespace {
template <int cao> struct AssignTorques {
  void operator()(dp3m_data_struct const &dp3m, double prefac, int d_rs,
                  ParticleRange const &particles) const {

    /* magnetic particle index */
    auto p_index = std::size_t{0ul};

    for (auto &p : particles) {
      if (p.dipm() != 0.) {
        auto const w = dp3m.inter_weights.load<cao>(p_index);

        Utils::Vector3d E{};
        p3m_interpolate(dp3m.local_mesh, w,
                        [&E, &dp3m, d_rs](int ind, double w) {
                          E[d_rs] += w * dp3m.rs_mesh[ind];
                        });

        p.torque() -= vector_product(p.calc_dip(), prefac * E);
        ++p_index;
      }
    }
  }
};

template <int cao> struct AssignForces {
  void operator()(dp3m_data_struct const &dp3m, double prefac, int d_rs,
                  ParticleRange const &particles) const {

    /* magnetic particle index */
    auto p_index = std::size_t{0ul};

    for (auto &p : particles) {
      if (p.dipm() != 0.) {
        auto const w = dp3m.inter_weights.load<cao>(p_index);

        Utils::Vector3d E{};
        p3m_interpolate(dp3m.local_mesh, w, [&E, &dp3m](int ind, double w) {
          E[0] += w * dp3m.rs_mesh_dip[0][ind];
          E[1] += w * dp3m.rs_mesh_dip[1][ind];
          E[2] += w * dp3m.rs_mesh_dip[2][ind];
        });

        p.force()[d_rs] += p.calc_dip() * prefac * E;
        ++p_index;
      }
    }
  }
};
} // namespace

double DipolarP3M::kernel(bool force_flag, bool energy_flag,
                          ParticleRange const &particles) {
  int i, ind, j[3];
  /* k-space energy */
  double k_space_energy_dip = 0.;
  double tmp0, tmp1;

  auto const dipole_prefac = prefactor / Utils::int_pow<3>(dp3m.params.mesh[0]);

  if (dp3m.sum_mu2 > 0) {
    /* Gather information for FFT grid inside the nodes domain (inner local
     * mesh) and perform forward 3D FFT (Charge Assignment Mesh). */
    std::array<double *, 3> meshes = {{dp3m.rs_mesh_dip[0].data(),
                                       dp3m.rs_mesh_dip[1].data(),
                                       dp3m.rs_mesh_dip[2].data()}};

    dp3m.sm.gather_grid(Utils::make_span(meshes), comm_cart,
                        dp3m.local_mesh.dim);

    fft_perform_forw(dp3m.rs_mesh_dip[0].data(), dp3m.fft, comm_cart);
    fft_perform_forw(dp3m.rs_mesh_dip[1].data(), dp3m.fft, comm_cart);
    fft_perform_forw(dp3m.rs_mesh_dip[2].data(), dp3m.fft, comm_cart);
    // Note: after these calls, the grids are in the order yzx and not xyz
    // anymore!!!
  }

  /* === k-space energy calculation  === */
  if (energy_flag) {
    /*********************
       Dipolar energy
    **********************/
    if (dp3m.sum_mu2 > 0) {
      /* i*k differentiation for dipolar gradients:
       * |(\Fourier{\vect{mu}}(k)\cdot \vect{k})|^2 */
      ind = 0;
      i = 0;
      double node_k_space_energy_dip = 0.0;
      for (j[0] = 0; j[0] < dp3m.fft.plan[3].new_mesh[0]; j[0]++) {
        for (j[1] = 0; j[1] < dp3m.fft.plan[3].new_mesh[1]; j[1]++) {
          for (j[2] = 0; j[2] < dp3m.fft.plan[3].new_mesh[2]; j[2]++) {
            node_k_space_energy_dip +=
                dp3m.g_energy[i] *
                (Utils::sqr(
                     dp3m.rs_mesh_dip[0][ind] *
                         dp3m.d_op[0][j[2] + dp3m.fft.plan[3].start[2]] +
                     dp3m.rs_mesh_dip[1][ind] *
                         dp3m.d_op[0][j[0] + dp3m.fft.plan[3].start[0]] +
                     dp3m.rs_mesh_dip[2][ind] *
                         dp3m.d_op[0][j[1] + dp3m.fft.plan[3].start[1]]) +
                 Utils::sqr(
                     dp3m.rs_mesh_dip[0][ind + 1] *
                         dp3m.d_op[0][j[2] + dp3m.fft.plan[3].start[2]] +
                     dp3m.rs_mesh_dip[1][ind + 1] *
                         dp3m.d_op[0][j[0] + dp3m.fft.plan[3].start[0]] +
                     dp3m.rs_mesh_dip[2][ind + 1] *
                         dp3m.d_op[0][j[1] + dp3m.fft.plan[3].start[1]]));
            ind += 2;
            i++;
          }
        }
      }
      node_k_space_energy_dip *=
          dipole_prefac * Utils::pi() * box_geo.length_inv()[0];
      boost::mpi::reduce(comm_cart, node_k_space_energy_dip, k_space_energy_dip,
                         std::plus<>(), 0);

      if (dp3m.energy_correction == 0.)
        calc_energy_correction();

      if (this_node == 0) {
        /* self energy correction */
        k_space_energy_dip -= prefactor * dp3m.sum_mu2 * 2. *
                              Utils::int_pow<3>(dp3m.params.alpha) *
                              Utils::sqrt_pi_i() / 3.;

        /* dipolar energy correction due to systematic Madelung-self effects */
        auto const volume = box_geo.volume();
        k_space_energy_dip += prefactor * dp3m.energy_correction / volume;
      }
    }
  } // if (energy_flag)

  /* === k-space force calculation  === */
  if (force_flag) {
    /****************************
     * DIPOLAR TORQUES (k-space)
     ****************************/
    if (dp3m.sum_mu2 > 0) {
      auto const two_pi_L_i = 2. * Utils::pi() * box_geo.length_inv()[0];
      /* fill in ks_mesh array for torque calculation */
      ind = 0;
      i = 0;

      for (j[0] = 0; j[0] < dp3m.fft.plan[3].new_mesh[0]; j[0]++) { // j[0]=n_y
        for (j[1] = 0; j[1] < dp3m.fft.plan[3].new_mesh[1];
             j[1]++) { // j[1]=n_z
          for (j[2] = 0; j[2] < dp3m.fft.plan[3].new_mesh[2];
               j[2]++) { // j[2]=n_x
            // tmp0 = Re(mu)*k,   tmp1 = Im(mu)*k

            tmp0 = dp3m.rs_mesh_dip[0][ind] *
                       dp3m.d_op[0][j[2] + dp3m.fft.plan[3].start[2]] +
                   dp3m.rs_mesh_dip[1][ind] *
                       dp3m.d_op[0][j[0] + dp3m.fft.plan[3].start[0]] +
                   dp3m.rs_mesh_dip[2][ind] *
                       dp3m.d_op[0][j[1] + dp3m.fft.plan[3].start[1]];

            tmp1 = dp3m.rs_mesh_dip[0][ind + 1] *
                       dp3m.d_op[0][j[2] + dp3m.fft.plan[3].start[2]] +
                   dp3m.rs_mesh_dip[1][ind + 1] *
                       dp3m.d_op[0][j[0] + dp3m.fft.plan[3].start[0]] +
                   dp3m.rs_mesh_dip[2][ind + 1] *
                       dp3m.d_op[0][j[1] + dp3m.fft.plan[3].start[1]];

            /* the optimal influence function is the same for torques
               and energy */
            dp3m.ks_mesh[ind] = tmp0 * dp3m.g_energy[i];
            dp3m.ks_mesh[ind + 1] = tmp1 * dp3m.g_energy[i];
            ind += 2;
            i++;
          }
        }
      }

      /* Force component loop */
      for (int d = 0; d < 3; d++) {
        auto const d_rs = (d + dp3m.ks_pnum) % 3;
        ind = 0;
        for (j[0] = 0; j[0] < dp3m.fft.plan[3].new_mesh[0]; j[0]++) {
          for (j[1] = 0; j[1] < dp3m.fft.plan[3].new_mesh[1]; j[1]++) {
            for (j[2] = 0; j[2] < dp3m.fft.plan[3].new_mesh[2]; j[2]++) {
              dp3m.rs_mesh[ind] =
                  dp3m.d_op[0][j[d] + dp3m.fft.plan[3].start[d]] *
                  dp3m.ks_mesh[ind];
              ind++;
              dp3m.rs_mesh[ind] =
                  dp3m.d_op[0][j[d] + dp3m.fft.plan[3].start[d]] *
                  dp3m.ks_mesh[ind];
              ind++;
            }
          }
        }

        /* Back FFT force component mesh */
        fft_perform_back(dp3m.rs_mesh.data(), false, dp3m.fft, comm_cart);
        /* redistribute force component mesh */
        dp3m.sm.spread_grid(dp3m.rs_mesh.data(), comm_cart,
                            dp3m.local_mesh.dim);
        /* Assign force component from mesh to particle */
        Utils::integral_parameter<int, AssignTorques, 1, 7>(
            dp3m.params.cao, dp3m, dipole_prefac * two_pi_L_i, d_rs, particles);
      }

      /***************************
         DIPOLAR FORCES (k-space)
      ****************************/

      // Compute forces after torques because the algorithm below overwrites the
      // grids dp3m.rs_mesh_dip !
      // Note: I'll do here 9 inverse FFTs. By symmetry, we can reduce this
      // number to 6 !
      /* fill in ks_mesh array for force calculation */
      ind = 0;
      i = 0;
      for (j[0] = 0; j[0] < dp3m.fft.plan[3].new_mesh[0]; j[0]++) { // j[0]=n_y
        for (j[1] = 0; j[1] < dp3m.fft.plan[3].new_mesh[1];
             j[1]++) { // j[1]=n_z
          for (j[2] = 0; j[2] < dp3m.fft.plan[3].new_mesh[2];
               j[2]++) { // j[2]=n_x
            // tmp0 = Im(mu)*k,   tmp1 = -Re(mu)*k
            tmp0 = dp3m.rs_mesh_dip[0][ind + 1] *
                       dp3m.d_op[0][j[2] + dp3m.fft.plan[3].start[2]] +
                   dp3m.rs_mesh_dip[1][ind + 1] *
                       dp3m.d_op[0][j[0] + dp3m.fft.plan[3].start[0]] +
                   dp3m.rs_mesh_dip[2][ind + 1] *
                       dp3m.d_op[0][j[1] + dp3m.fft.plan[3].start[1]];
            tmp1 = dp3m.rs_mesh_dip[0][ind] *
                       dp3m.d_op[0][j[2] + dp3m.fft.plan[3].start[2]] +
                   dp3m.rs_mesh_dip[1][ind] *
                       dp3m.d_op[0][j[0] + dp3m.fft.plan[3].start[0]] +
                   dp3m.rs_mesh_dip[2][ind] *
                       dp3m.d_op[0][j[1] + dp3m.fft.plan[3].start[1]];
            dp3m.ks_mesh[ind] = tmp0 * dp3m.g_force[i];
            dp3m.ks_mesh[ind + 1] = -tmp1 * dp3m.g_force[i];
            ind += 2;
            i++;
          }
        }
      }

      /* Force component loop */
      for (int d = 0; d < 3; d++) { /* direction in k-space: */
        auto const d_rs = (d + dp3m.ks_pnum) % 3;
        ind = 0;
        for (j[0] = 0; j[0] < dp3m.fft.plan[3].new_mesh[0];
             j[0]++) { // j[0]=n_y
          for (j[1] = 0; j[1] < dp3m.fft.plan[3].new_mesh[1];
               j[1]++) { // j[1]=n_z
            for (j[2] = 0; j[2] < dp3m.fft.plan[3].new_mesh[2];
                 j[2]++) { // j[2]=n_x
              tmp0 = dp3m.d_op[0][j[d] + dp3m.fft.plan[3].start[d]] *
                     dp3m.ks_mesh[ind];
              dp3m.rs_mesh_dip[0][ind] =
                  dp3m.d_op[0][j[2] + dp3m.fft.plan[3].start[2]] * tmp0;
              dp3m.rs_mesh_dip[1][ind] =
                  dp3m.d_op[0][j[0] + dp3m.fft.plan[3].start[0]] * tmp0;
              dp3m.rs_mesh_dip[2][ind] =
                  dp3m.d_op[0][j[1] + dp3m.fft.plan[3].start[1]] * tmp0;
              ind++;
              tmp0 = dp3m.d_op[0][j[d] + dp3m.fft.plan[3].start[d]] *
                     dp3m.ks_mesh[ind];
              dp3m.rs_mesh_dip[0][ind] =
                  dp3m.d_op[0][j[2] + dp3m.fft.plan[3].start[2]] * tmp0;
              dp3m.rs_mesh_dip[1][ind] =
                  dp3m.d_op[0][j[0] + dp3m.fft.plan[3].start[0]] * tmp0;
              dp3m.rs_mesh_dip[2][ind] =
                  dp3m.d_op[0][j[1] + dp3m.fft.plan[3].start[1]] * tmp0;
              ind++;
            }
          }
        }
        /* Back FFT force component mesh */
        fft_perform_back(dp3m.rs_mesh_dip[0].data(), false, dp3m.fft,
                         comm_cart);
        fft_perform_back(dp3m.rs_mesh_dip[1].data(), false, dp3m.fft,
                         comm_cart);
        fft_perform_back(dp3m.rs_mesh_dip[2].data(), false, dp3m.fft,
                         comm_cart);
        /* redistribute force component mesh */
        std::array<double *, 3> meshes = {{dp3m.rs_mesh_dip[0].data(),
                                           dp3m.rs_mesh_dip[1].data(),
                                           dp3m.rs_mesh_dip[2].data()}};

        dp3m.sm.spread_grid(Utils::make_span(meshes), comm_cart,
                            dp3m.local_mesh.dim);
        /* Assign force component from mesh to particle */
        Utils::integral_parameter<int, AssignForces, 1, 7>(
            dp3m.params.cao, dp3m, dipole_prefac * Utils::sqr(two_pi_L_i), d_rs,
            particles);
      }
    } /* if (dp3m.sum_mu2 > 0) */
  }   /* if (force_flag) */

  if (dp3m.params.epsilon != P3M_EPSILON_METALLIC) {
    auto const surface_term =
        calc_surface_term(force_flag, energy_flag, particles);
    if (this_node == 0)
      k_space_energy_dip += surface_term;
  }

  return k_space_energy_dip;
}

double DipolarP3M::calc_surface_term(bool force_flag, bool energy_flag,
                                     ParticleRange const &particles) {
  auto const pref = prefactor * 4. * Utils::pi() / box_geo.volume() /
                    (2. * dp3m.params.epsilon + 1.);
  auto const n_local_part = particles.size();

  // We put all the dipolar momenta in a the arrays mx,my,mz according to the
  // id-number of the particles
  std::vector<double> mx(n_local_part);
  std::vector<double> my(n_local_part);
  std::vector<double> mz(n_local_part);

  int ip = 0;
  for (auto const &p : particles) {
    auto const dip = p.calc_dip();
    mx[ip] = dip[0];
    my[ip] = dip[1];
    mz[ip] = dip[2];
    ip++;
  }

  // we will need the sum of all dipolar momenta vectors
  auto local_dip = Utils::Vector3d{};
  for (int i = 0; i < n_local_part; i++) {
    local_dip[0] += mx[i];
    local_dip[1] += my[i];
    local_dip[2] += mz[i];
  }
  auto const box_dip =
      boost::mpi::all_reduce(comm_cart, local_dip, std::plus<>());

  double energy = 0.;
  if (energy_flag) {
    double sum_e = 0.;
    for (int i = 0; i < n_local_part; i++) {
      sum_e += mx[i] * box_dip[0] + my[i] * box_dip[1] + mz[i] * box_dip[2];
    }
    energy =
        0.5 * pref * boost::mpi::all_reduce(comm_cart, sum_e, std::plus<>());
  }

  if (force_flag) {

    std::vector<double> sumix(n_local_part);
    std::vector<double> sumiy(n_local_part);
    std::vector<double> sumiz(n_local_part);

    for (int i = 0; i < n_local_part; i++) {
      sumix[i] = my[i] * box_dip[2] - mz[i] * box_dip[1];
      sumiy[i] = mz[i] * box_dip[0] - mx[i] * box_dip[2];
      sumiz[i] = mx[i] * box_dip[1] - my[i] * box_dip[0];
    }

    ip = 0;
    for (auto &p : particles) {
      auto &torque = p.torque();
      torque[0] -= pref * sumix[ip];
      torque[1] -= pref * sumiy[ip];
      torque[2] -= pref * sumiz[ip];
      ip++;
    }
  }

  return energy;
}

void DipolarP3M::calc_influence_function_force() {
  auto const start = Utils::Vector3i{dp3m.fft.plan[3].start};
  auto const size = Utils::Vector3i{dp3m.fft.plan[3].new_mesh};

  dp3m.g_force = grid_influence_function<3>(dp3m.params, start, start + size,
                                            box_geo.length());
}

void DipolarP3M::calc_influence_function_energy() {
  auto const start = Utils::Vector3i{dp3m.fft.plan[3].start};
  auto const size = Utils::Vector3i{dp3m.fft.plan[3].new_mesh};

  dp3m.g_energy = grid_influence_function<2>(dp3m.params, start, start + size,
                                             box_geo.length());
}

class DipolarTuningAlgorithm : public TuningAlgorithm {
  dp3m_data_struct &dp3m;
  int m_mesh_max = -1, m_mesh_min = -1;

public:
  DipolarTuningAlgorithm(dp3m_data_struct &input_dp3m, double prefactor,
                         int timings)
      : TuningAlgorithm{prefactor, timings}, dp3m{input_dp3m} {}

  P3MParameters &get_params() override { return dp3m.params; }

  void on_solver_change() const override { on_dipoles_change(); }

  boost::optional<std::string>
  layer_correction_veto_r_cut(double) const override {
    return {};
  }

  void setup_logger(bool verbose) override {
    m_logger = std::make_unique<TuningLogger>(
        verbose and this_node == 0, "DipolarP3M", TuningLogger::Mode::Dipolar);
    m_logger->tuning_goals(dp3m.params.accuracy, m_prefactor,
                           box_geo.length()[0], dp3m.sum_dip_part,
                           dp3m.sum_mu2);
    m_logger->log_tuning_start();
  }

  std::tuple<double, double, double, double>
  calculate_accuracy(Utils::Vector3i const &mesh, int cao,
                     double r_cut_iL) const override {

    double alpha_L, rs_err, ks_err;

    /* calc maximal real space error for setting */
    rs_err = dp3m_real_space_error(box_geo.length()[0], r_cut_iL,
                                   dp3m.sum_dip_part, dp3m.sum_mu2, 0.001);
    // alpha cannot be zero for dipoles because real-space formula breaks down

    if (Utils::sqrt_2() * rs_err > dp3m.params.accuracy) {
      /* assume rs_err = ks_err -> rs_err = accuracy/sqrt(2.0) -> alpha_L */
      alpha_L = dp3m_rtbisection(
          box_geo.length()[0], r_cut_iL, dp3m.sum_dip_part, dp3m.sum_mu2,
          0.0001 * box_geo.length()[0], 5. * box_geo.length()[0], 0.0001,
          dp3m.params.accuracy);
    } else {
      /* even alpha=0 is ok, however, we cannot choose it since it kills the
         k-space error formula.
         Anyways, this very likely NOT the optimal solution */
      alpha_L = 0.1;
    }

    /* calculate real-space and k-space error for this alpha_L */
    rs_err = dp3m_real_space_error(box_geo.length()[0], r_cut_iL,
                                   dp3m.sum_dip_part, dp3m.sum_mu2, alpha_L);
    ks_err = dp3m_k_space_error(box_geo.length()[0], mesh[0], cao,
                                dp3m.sum_dip_part, dp3m.sum_mu2, alpha_L);

    return {Utils::Vector2d{rs_err, ks_err}.norm(), rs_err, ks_err, alpha_L};
  }

  void determine_mesh_limits() override {
    if (dp3m.params.mesh[0] == -1) {
      /* simple heuristic to limit the tried meshes if the accuracy cannot
         be obtained with smaller meshes, but normally not all these
         meshes have to be tested */
      auto const expo = std::log(std::cbrt(dp3m.sum_dip_part)) / std::log(2.);
      /* Medium-educated guess for the minimal mesh */
      m_mesh_min = static_cast<int>(std::round(std::pow(2., std::floor(expo))));
      /* avoid using more than 1 GB of FFT arrays */
      m_mesh_max = 128;
    } else {
      m_mesh_min = m_mesh_max = dp3m.params.mesh[0];
      m_logger->report_fixed_mesh(dp3m.params.mesh);
    }
  }

  TuningAlgorithm::Parameters get_time() override {
    auto tuned_params = TuningAlgorithm::Parameters{};
    auto time_best = time_sentinel;
    for (auto tmp_mesh = m_mesh_min; tmp_mesh <= m_mesh_max; tmp_mesh += 2) {
      auto trial_params = TuningAlgorithm::Parameters{};
      trial_params.mesh = Utils::Vector3i::broadcast(tmp_mesh);
      trial_params.cao = cao_best;

      auto const trial_time =
          get_m_time(trial_params.mesh, trial_params.cao, trial_params.r_cut_iL,
                     trial_params.alpha_L, trial_params.accuracy);

      /* this mesh does not work at all */
      if (trial_time < 0.)
        continue;

      /* the optimum r_cut for this mesh is the upper limit for higher meshes,
         everything else is slower */
      m_r_cut_iL_max = trial_params.r_cut_iL;

      if (trial_time < time_best) {
        /* new optimum */
        reset_n_trials();
        tuned_params = trial_params;
        time_best = tuned_params.time = trial_time;
      } else if (trial_time > time_best + time_granularity or
                 get_n_trials() > max_n_consecutive_trials) {
        /* no hope of further optimisation */
        break;
      }
    }
    return tuned_params;
  }
};

void DipolarP3M::tune() {
  if (dp3m.params.alpha_L == 0. and dp3m.params.alpha != 0.) {
    dp3m.params.alpha_L = dp3m.params.alpha * box_geo.length()[0];
  }
  if (dp3m.params.r_cut_iL == 0. and dp3m.params.r_cut != 0.) {
    dp3m.params.r_cut_iL = dp3m.params.r_cut * box_geo.length_inv()[0];
  }
  if (not is_tuned()) {
    count_magnetic_particles();
    if (dp3m.sum_dip_part == 0) {
      throw std::runtime_error(
          "DipolarP3M: no dipolar particles in the system");
    }
    try {
      DipolarTuningAlgorithm parameters(dp3m, prefactor, tune_timings);
      parameters.setup_logger(tune_verbose);
      // parameter ranges
      parameters.determine_mesh_limits();
      parameters.determine_r_cut_limits();
      parameters.determine_cao_limits(3);
      // run tuning algorithm
      parameters.tune();
      m_is_tuned = true;
      on_dipoles_change();
    } catch (...) {
      dp3m.params.tuning = false;
      throw;
    }
  }
  init();
}

/** Calculate the k-space error of dipolar-P3M */
static double dp3m_k_space_error(double box_size, int mesh, int cao,
                                 int n_c_part, double sum_q2, double alpha_L) {
  double he_q = 0.;
  auto const mesh_i = 1. / mesh;
  auto const alpha_L_i = 1. / alpha_L;

  for (int nx = -mesh / 2; nx < mesh / 2; nx++)
    for (int ny = -mesh / 2; ny < mesh / 2; ny++)
      for (int nz = -mesh / 2; nz < mesh / 2; nz++)
        if ((nx != 0) || (ny != 0) || (nz != 0)) {
          auto const n2 = Utils::sqr(nx) + Utils::sqr(ny) + Utils::sqr(nz);
          auto const cs = p3m_analytic_cotangent_sum(nx, mesh_i, cao) *
                          p3m_analytic_cotangent_sum(ny, mesh_i, cao) *
                          p3m_analytic_cotangent_sum(nz, mesh_i, cao);
          double alias1, alias2;
          dp3m_tune_aliasing_sums(nx, ny, nz, mesh, mesh_i, cao, alpha_L_i,
                                  &alias1, &alias2);
          double d = alias1 - Utils::sqr(alias2 / cs) /
                                  Utils::int_pow<3>(static_cast<double>(n2));
          /* at high precision, d can become negative due to extinction;
             also, don't take values that have no significant digits left*/
          if (d > 0 && (fabs(d / alias1) > ROUND_ERROR_PREC))
            he_q += d;
        }

  return 8. * Utils::sqr(Utils::pi()) / 3. * sum_q2 * sqrt(he_q / n_c_part) /
         Utils::int_pow<4>(box_size);
}

/** Tuning dipolar-P3M */
void dp3m_tune_aliasing_sums(int nx, int ny, int nz, int mesh, double mesh_i,
                             int cao, double alpha_L_i, double *alias1,
                             double *alias2) {
  using Utils::sinc;

  auto const factor1 = Utils::sqr(Utils::pi() * alpha_L_i);

  *alias1 = *alias2 = 0.;
  for (int mx = -P3M_BRILLOUIN; mx <= P3M_BRILLOUIN; mx++) {
    auto const nmx = nx + mx * mesh;
    auto const fnmx = mesh_i * nmx;
    for (int my = -P3M_BRILLOUIN; my <= P3M_BRILLOUIN; my++) {
      auto const nmy = ny + my * mesh;
      auto const fnmy = mesh_i * nmy;
      for (int mz = -P3M_BRILLOUIN; mz <= P3M_BRILLOUIN; mz++) {
        auto const nmz = nz + mz * mesh;
        auto const fnmz = mesh_i * nmz;

        auto const nm2 = Utils::sqr(nmx) + Utils::sqr(nmy) + Utils::sqr(nmz);
        auto const ex = std::exp(-factor1 * nm2);

        auto const U2 = pow(sinc(fnmx) * sinc(fnmy) * sinc(fnmz), 2. * cao);

        *alias1 += Utils::sqr(ex) * nm2;
        *alias2 += U2 * ex * pow((nx * nmx + ny * nmy + nz * nmz), 3.) / nm2;
      }
    }
  }
}

/** Calculate the value of the errors for the REAL part of the force in terms
 *  of the splitting parameter alpha of Ewald. Based on eq. (33) @cite wang01a.
 *
 *  Please note that in this more refined approach we don't use
 *  eq. (37), but eq. (33) which maintains all the powers in alpha.
 */
double dp3m_real_space_error(double box_size, double r_cut_iL, int n_c_part,
                             double sum_q2, double alpha_L) {
  double d_error_f, d_cc, d_dc, d_con;

  auto const d_rcut = r_cut_iL * box_size;
  auto const d_rcut2 = Utils::sqr(d_rcut);
  auto const d_rcut4 = Utils::sqr(d_rcut2);

  auto const d_a2 = Utils::sqr(alpha_L) / Utils::sqr(box_size);

  auto const d_c = sum_q2 * exp(-d_a2 * d_rcut2);

  d_cc = 4. * Utils::sqr(d_a2) * Utils::sqr(d_rcut2) + 6. * d_a2 * d_rcut2 + 3.;

  d_dc = 8. * Utils::int_pow<3>(d_a2) * Utils::int_pow<3>(d_rcut2) +
         20. * Utils::sqr(d_a2) * d_rcut4 + 30. * d_a2 * d_rcut2 + 15.;

  d_con = 1. / sqrt(Utils::int_pow<3>(box_size) * Utils::sqr(d_a2) * d_rcut *
                    Utils::sqr(d_rcut4) * static_cast<double>(n_c_part));

  d_error_f = d_c * d_con *
              sqrt((13. / 6.) * Utils::sqr(d_cc) +
                   (2. / 15.) * Utils::sqr(d_dc) - (13. / 15.) * d_cc * d_dc);

  return d_error_f;
}

/** Using bisection, find the root of a function "func-tuned_accuracy/sqrt(2.)"
 *  known to lie between x1 and x2. The root, returned as rtbis, will be
 *  refined until its accuracy is \f$\pm\f$ @p xacc.
 */
double dp3m_rtbisection(double box_size, double r_cut_iL, int n_c_part,
                        double sum_q2, double x1, double x2, double xacc,
                        double tuned_accuracy) {
  constexpr int JJ_RTBIS_MAX = 40;

  auto const constant = tuned_accuracy / Utils::sqrt_2();

  auto const f1 =
      dp3m_real_space_error(box_size, r_cut_iL, n_c_part, sum_q2, x1) -
      constant;
  auto const f2 =
      dp3m_real_space_error(box_size, r_cut_iL, n_c_part, sum_q2, x2) -
      constant;
  if (f1 * f2 >= 0.0) {
    throw std::runtime_error(
        "Root must be bracketed for bisection in dp3m_rtbisection");
  }
  // Orient the search dx, and set rtb to x1 or x2 ...
  double dx;
  double rtb = f1 < 0.0 ? (dx = x2 - x1, x1) : (dx = x1 - x2, x2);
  for (int j = 1; j <= JJ_RTBIS_MAX; j++) {
    auto const xmid = rtb + (dx *= 0.5);
    auto const fmid =
        dp3m_real_space_error(box_size, r_cut_iL, n_c_part, sum_q2, xmid) -
        constant;
    if (fmid <= 0.0)
      rtb = xmid;
    if (fabs(dx) < xacc || fmid == 0.0)
      return rtb;
  }
  throw std::runtime_error("Too many bisections in dp3m_rtbisection");
}

void DipolarP3M::sanity_checks_boxl() const {
  for (unsigned int i = 0; i < 3; i++) {
    /* check k-space cutoff */
    if (dp3m.params.cao_cut[i] >= box_geo.length_half()[i]) {
      std::stringstream msg;
      msg << "dipolar P3M_init: k-space cutoff " << dp3m.params.cao_cut[i]
          << " is larger than half of box dimension " << box_geo.length()[i];
      throw std::runtime_error(msg.str());
    }
    if (dp3m.params.cao_cut[i] >= local_geo.length()[i]) {
      std::stringstream msg;
      msg << "dipolar P3M_init: k-space cutoff " << dp3m.params.cao_cut[i]
          << " is larger than local box dimension " << local_geo.length()[i];
      throw std::runtime_error(msg.str());
    }
  }

  if ((box_geo.length()[0] != box_geo.length()[1]) or
      (box_geo.length()[1] != box_geo.length()[2])) {
    throw std::runtime_error("DipolarP3M: requires a cubic box");
  }
}

void DipolarP3M::sanity_checks_periodicity() const {
  if (!box_geo.periodic(0) or !box_geo.periodic(1) or !box_geo.periodic(2)) {
    throw std::runtime_error(
        "DipolarP3M: requires periodicity (True, True, True)");
  }
}

void DipolarP3M::sanity_checks_cell_structure() const {
  if (local_geo.cell_structure_type() !=
          CellStructureType::CELL_STRUCTURE_REGULAR &&
      local_geo.cell_structure_type() !=
          CellStructureType::CELL_STRUCTURE_HYBRID) {
    throw std::runtime_error(
        "DipolarP3M: requires the regular or hybrid decomposition cell system");
  }
  if (n_nodes > 1 && local_geo.cell_structure_type() ==
                         CellStructureType::CELL_STRUCTURE_HYBRID) {
    throw std::runtime_error(
        "DipolarP3M: does not work with the hybrid decomposition cell system, "
        "if using more than one MPI node");
  }
}

void DipolarP3M::sanity_checks_node_grid() const {
  if (node_grid[0] < node_grid[1] || node_grid[1] < node_grid[2]) {
    throw std::runtime_error(
        "DipolarP3M: node grid must be sorted, largest first");
  }
}

void DipolarP3M::scaleby_box_l() {
  dp3m.params.r_cut = dp3m.params.r_cut_iL * box_geo.length()[0];
  dp3m.params.alpha = dp3m.params.alpha_L * box_geo.length_inv()[0];
  dp3m.params.recalc_a_ai_cao_cut(box_geo.length());
  dp3m.local_mesh.recalc_ld_pos(dp3m.params);
  sanity_checks_boxl();
  calc_influence_function_force();
  calc_influence_function_energy();
  dp3m.energy_correction = 0.0;
}

void DipolarP3M::calc_energy_correction() {
  auto const Ukp3m = calc_average_self_energy_k_space() * box_geo.volume();
  auto const Ewald_volume = Utils::int_pow<3>(dp3m.params.alpha_L);
  auto const Eself = -2. * Ewald_volume * Utils::sqrt_pi_i() / 3.;
  dp3m.energy_correction =
      -dp3m.sum_mu2 * (Ukp3m + Eself + 2. * Utils::pi() / 3.);
}

#ifdef NPT
void npt_add_virial_magnetic_contribution(double energy) {
  npt_add_virial_contribution(energy);
}
#endif // NPT
#endif // DP3M
