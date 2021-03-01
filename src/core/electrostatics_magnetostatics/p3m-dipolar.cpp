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
 *  P3M algorithm for long range magnetic dipole-dipole interaction.
 *
 *  @note
 *      In general the magnetic dipole-dipole functions bear the same name than
 *      the charge-charge, but adding a "D" in front of the name and replacing
 *      "charge" by "dipole". In this way one can recognize the similarity of
 *      the functions while avoiding nasty confusions in their use.
 *
 *  By default the magnetic epsilon is metallic = 0.
 *
 *  The corresponding header file is p3m-dipolar.hpp.
 */

#include "config.hpp"

#ifdef DP3M

#include "electrostatics_magnetostatics/p3m-dipolar.hpp"

#include "electrostatics_magnetostatics/common.hpp"
#include "electrostatics_magnetostatics/dp3m_influence_function.hpp"
#include "electrostatics_magnetostatics/fft.hpp"
#include "electrostatics_magnetostatics/p3m-common.hpp"
#include "electrostatics_magnetostatics/p3m_interpolation.hpp"
#include "electrostatics_magnetostatics/p3m_send_mesh.hpp"

#include "Particle.hpp"
#include "ParticleRange.hpp"
#include "cells.hpp"
#include "communication.hpp"
#include "errorhandling.hpp"
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
#include <boost/range/algorithm/min_element.hpp>

#include <algorithm>
#include <array>
#include <cstdio>
#include <functional>
#include <stdexcept>
#include <vector>

/************************************************
 * DEFINES
 ************************************************/

#define DP3M_RTBISECTION_ERROR 9999999

/************************************************
 * variables
 ************************************************/

dp3m_data_struct dp3m;

/** \name Private Functions */
/**@{*/

/** Initialize for magnetic dipoles the (inverse) mesh constant @ref
 *  P3MParameters::a "a" (@ref P3MParameters::ai "ai") and the
 *  cutoff for charge assignment @ref P3MParameters::cao_cut "cao_cut".
 *
 *  Function called by @ref dp3m_init() once and by @ref
 *  dp3m_scaleby_box_l() whenever the box_length changes.
 */
static void dp3m_init_a_ai_cao_cut();

/** Checks for correctness for magnetic dipoles in P3M of the cao_cut,
 *  necessary when the box length changes
 */
static bool dp3m_sanity_checks_boxl();

/** Calculate the influence function optimized for the dipolar forces. */
static void dp3m_calc_influence_function_force();

/** Calculate the influence function optimized for the dipolar energy and
 *  torques.
 */
static void dp3m_calc_influence_function_energy();

/** Calculate the constants necessary to correct the dipolar energy to minimize
 *  the error.
 */
static void dp3m_compute_constants_energy_dipolar();

static double dp3m_k_space_error(double box_size, double prefac, int mesh,
                                 int cao, int n_c_part, double sum_q2,
                                 double alpha_L);
/**@}*/

/** Compute the dipolar surface terms */
static double calc_surface_term(bool force_flag, bool energy_flag,
                                ParticleRange const &particles);

/** \name P3M Tuning Functions */
/************************************************************/
/**@{*/

double dp3m_real_space_error(double box_size, double prefac, double r_cut_iL,
                             int n_c_part, double sum_q2, double alpha_L);
static void dp3m_tune_aliasing_sums(int nx, int ny, int nz, int mesh,
                                    double mesh_i, int cao, double alpha_L_i,
                                    double *alias1, double *alias2);

/** Compute the value of alpha through a bisection method.
 *  Based on eq. (33) @cite wang01a.
 */
double dp3m_rtbisection(double box_size, double prefac, double r_cut_iL,
                        int n_c_part, double sum_q2, double x1, double x2,
                        double xacc, double tuned_accuracy);

/**@}*/

/** Correction of the dipolar p3m-energy. */
double dp3m_average_dipolar_self_energy() {
  auto const start = Utils::Vector3i{dp3m.fft.plan[3].start};
  auto const size = Utils::Vector3i{dp3m.fft.plan[3].new_mesh};

  auto const node_phi = grid_influence_function_self_energy(
      dp3m.params, start, start + size, dp3m.g_energy);

  double phi = 0.0;
  boost::mpi::reduce(comm_cart, node_phi, phi, std::plus<>(), 0);
  phi /= 3. * box_geo.length()[0] * Utils::int_pow<3>(dp3m.params.mesh[0]);
  return phi * Utils::pi();
}

/************************************************************/

dp3m_data_struct::dp3m_data_struct() {
  /* local_mesh is uninitialized */
  /* sm is uninitialized */
  sum_dip_part = 0;
  sum_mu2 = 0.0;

  pos_shift = 0.0;
  ks_pnum = 0;

  energy_correction = 0.0;
}

void dp3m_deactivate() {
  dp3m.params.alpha = 0.0;
  dp3m.params.alpha_L = 0.0;
  dp3m.params.r_cut = 0.0;
  dp3m.params.r_cut_iL = 0.0;
  dp3m.params.mesh[0] = 0;
  dp3m.params.mesh[1] = 0;
  dp3m.params.mesh[2] = 0;
  dp3m.params.cao = 0;
}

void dp3m_init() {
  if (dipole.prefactor <= 0.0) {
    // dipolar prefactor is zero: magnetostatics switched off
    dp3m.params.r_cut = 0.0;
    dp3m.params.r_cut_iL = 0.0;
    return;
  }

  if (dp3m_sanity_checks(node_grid)) {
    return;
  }

  dp3m.params.cao3 = Utils::int_pow<3>(dp3m.params.cao);

  /* initializes the (inverse) mesh constant dp3m.params.a (dp3m.params.ai)
   * and the cutoff for charge assignment dp3m.params.cao_cut */
  dp3m_init_a_ai_cao_cut();

  p3m_calc_local_ca_mesh(dp3m.local_mesh, dp3m.params, local_geo, skin, 0.0);

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
  dp3m_scaleby_box_l();

  dp3m_count_magnetic_particles();
}

/******************
 * functions related to the parsing & tuning of the dipolar parameters
 ******************/

void dp3m_set_tune_params(double r_cut, int mesh, int cao, double accuracy) {
  if (r_cut >= 0) {
    dp3m.params.r_cut = r_cut;
    dp3m.params.r_cut_iL = r_cut / box_geo.length()[0];
  }

  if (mesh >= 0)
    dp3m.params.mesh[2] = dp3m.params.mesh[1] = dp3m.params.mesh[0] = mesh;

  if (cao >= 0)
    dp3m.params.cao = cao;

  if (accuracy >= 0)
    dp3m.params.accuracy = accuracy;
}

/*****************************************************************************/

void dp3m_set_params(double r_cut, int mesh, int cao, double alpha,
                     double accuracy) {
  if (r_cut < 0)
    throw std::runtime_error("DipolarP3M: invalid r_cut");

  if (mesh < 0)
    throw std::runtime_error("DipolarP3M: invalid mesh size");

  if (cao < 1 || cao > 7)
    throw std::runtime_error("DipolarP3M: invalid cao");

  if (cao > mesh)
    throw std::runtime_error("DipolarP3M: cao larger than mesh size");

  if (alpha <= 0.0 && alpha != -1.0)
    throw std::runtime_error("DipolarP3M: invalid alpha");

  if (accuracy <= 0.0 && accuracy != -1.0)
    throw std::runtime_error("DipolarP3M: invalid accuracy");

  if (dipole.method != DIPOLAR_P3M && dipole.method != DIPOLAR_MDLC_P3M)
    Dipole::set_method_local(DIPOLAR_P3M);

  dp3m.params.r_cut = r_cut;
  dp3m.params.r_cut_iL = r_cut / box_geo.length()[0];
  dp3m.params.mesh[2] = dp3m.params.mesh[1] = dp3m.params.mesh[0] = mesh;
  dp3m.params.cao = cao;
  dp3m.params.alpha = alpha;
  dp3m.params.alpha_L = alpha * box_geo.length()[0];
  dp3m.params.accuracy = accuracy;

  mpi_bcast_coulomb_params();
}

void dp3m_set_mesh_offset(double x, double y, double z) {
  if (x == -1.0 && y == -1.0 && z == -1.0)
    return;

  if (x < 0.0 || x > 1.0 || y < 0.0 || y > 1.0 || z < 0.0 || z > 1.0)
    throw std::runtime_error("DipolarP3M: invalid mesh offset");

  dp3m.params.mesh_off[0] = x;
  dp3m.params.mesh_off[1] = y;
  dp3m.params.mesh_off[2] = z;

  mpi_bcast_coulomb_params();
}

/** We left the handling of the epsilon, due to portability reasons in
 *  the future for the electrical dipoles, or if people want to do
 *  electrical dipoles alone using the magnetic code. Currently unused.
 */
void dp3m_set_eps(double eps) {
  dp3m.params.epsilon = eps;

  mpi_bcast_coulomb_params();
}

namespace {
template <size_t cao> struct AssignDipole {
  void operator()(Utils::Vector3d const &real_pos,
                  Utils::Vector3d const &dip) const {
    auto const weights = p3m_calculate_interpolation_weights<cao>(
        real_pos, dp3m.params.ai, dp3m.local_mesh);
    p3m_interpolate<cao>(dp3m.local_mesh, weights, [&dip](int ind, double w) {
      dp3m.rs_mesh_dip[0][ind] += w * dip[0];
      dp3m.rs_mesh_dip[1][ind] += w * dip[1];
      dp3m.rs_mesh_dip[2][ind] += w * dip[2];
    });

    dp3m.inter_weights.store<cao>(weights);
  }
};
} // namespace

void dp3m_dipole_assign(const ParticleRange &particles) {
  dp3m.inter_weights.reset(dp3m.params.cao);

  /* prepare local FFT mesh */
  for (auto &i : dp3m.rs_mesh_dip)
    for (int j = 0; j < dp3m.local_mesh.size; j++)
      i[j] = 0.0;

  for (auto const &p : particles) {
    if (p.p.dipm != 0.0) {
      Utils::integral_parameter<AssignDipole, 1, 7>(dp3m.params.cao, p.r.p,
                                                    p.calc_dip());
    }
  }
}

namespace {
template <size_t cao> struct AssignTorques {
  void operator()(double prefac, int d_rs,
                  const ParticleRange &particles) const {
    /* particle counter */
    int cp_cnt = 0;
    for (auto &p : particles) {
      auto const w = dp3m.inter_weights.load<cao>(cp_cnt++);

      Utils::Vector3d E{};
      p3m_interpolate(dp3m.local_mesh, w, [&E, d_rs](int ind, double w) {
        E[d_rs] += w * dp3m.rs_mesh[ind];
      });

      p.f.torque -= vector_product(p.calc_dip(), prefac * E);
    }
  }
};

template <size_t cao> struct AssignForces {
  void operator()(double prefac, int d_rs,
                  const ParticleRange &particles) const {
    /* particle counter */
    int cp_cnt = 0;
    for (auto &p : particles) {
      auto const w = dp3m.inter_weights.load<cao>(cp_cnt++);

      Utils::Vector3d E{};

      p3m_interpolate(dp3m.local_mesh, w, [&E](int ind, double w) {
        E[0] += w * dp3m.rs_mesh_dip[0][ind];
        E[1] += w * dp3m.rs_mesh_dip[1][ind];
        E[2] += w * dp3m.rs_mesh_dip[2][ind];
      });

      p.f.f[d_rs] += p.calc_dip() * prefac * E;
    }
  }
};
} // namespace

/*****************************************************************************/

double dp3m_calc_kspace_forces(bool force_flag, bool energy_flag,
                               const ParticleRange &particles) {
  int i, d, d_rs, ind, j[3];
  /* k-space energy */
  double k_space_energy_dip = 0.0;
  double tmp0, tmp1;

  auto const dipole_prefac =
      dipole.prefactor / Utils::int_pow<3>(dp3m.params.mesh[0]);

  if (dp3m.sum_mu2 > 0) {
    /* Gather information for FFT grid inside the nodes domain (inner local
     * mesh) and perform forward 3D FFT (Charge Assignment Mesh). */
    std::array<double *, 3> meshes = {dp3m.rs_mesh_dip[0].data(),
                                      dp3m.rs_mesh_dip[1].data(),
                                      dp3m.rs_mesh_dip[2].data()};

    dp3m.sm.gather_grid(Utils::make_span(meshes), comm_cart,
                        dp3m.local_mesh.dim);

    fft_perform_forw(dp3m.rs_mesh_dip[0].data(), dp3m.fft, comm_cart);
    fft_perform_forw(dp3m.rs_mesh_dip[1].data(), dp3m.fft, comm_cart);
    fft_perform_forw(dp3m.rs_mesh_dip[2].data(), dp3m.fft, comm_cart);
    // Note: after these calls, the grids are in the order yzx and not xyz
    // anymore!!!
  }

  /* === k-space calculations === */

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
          dipole_prefac * Utils::pi() / box_geo.length()[0];
      boost::mpi::reduce(comm_cart, node_k_space_energy_dip, k_space_energy_dip,
                         std::plus<>(), 0);

      dp3m_compute_constants_energy_dipolar();

      if (this_node == 0) {
        /* self energy correction */
        k_space_energy_dip -=
            dipole.prefactor *
            (dp3m.sum_mu2 * 2 *
             pow(dp3m.params.alpha_L / box_geo.length()[0], 3) *
             Utils::sqrt_pi_i() / 3.0);

        auto const volume = box_geo.volume();
        k_space_energy_dip += dipole.prefactor * dp3m.energy_correction /
                              volume; /* add the dipolar energy correction due
                                         to systematic Madelung-Self effects */
      }
    }
  } // if (energy_flag)

  /* === k-space force calculation  === */
  if (force_flag) {
    /****************************
     * DIPOLAR TORQUES (k-space)
     ****************************/
    if (dp3m.sum_mu2 > 0) {
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
      for (d = 0; d < 3; d++) {
        d_rs = (d + dp3m.ks_pnum) % 3;
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
        Utils::integral_parameter<AssignTorques, 1, 7>(
            dp3m.params.cao,
            dipole_prefac * (2 * Utils::pi() / box_geo.length()[0]), d_rs,
            particles);
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
      for (d = 0; d < 3; d++) { /* direction in k-space: */
        d_rs = (d + dp3m.ks_pnum) % 3;
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
        std::array<double *, 3> meshes = {dp3m.rs_mesh_dip[0].data(),
                                          dp3m.rs_mesh_dip[1].data(),
                                          dp3m.rs_mesh_dip[2].data()};

        dp3m.sm.spread_grid(Utils::make_span(meshes), comm_cart,
                            dp3m.local_mesh.dim);
        /* Assign force component from mesh to particle */
        Utils::integral_parameter<AssignForces, 1, 7>(
            dp3m.params.cao,
            dipole_prefac * pow(2 * Utils::pi() / box_geo.length()[0], 2), d_rs,
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

/************************************************************/

double calc_surface_term(bool force_flag, bool energy_flag,
                         const ParticleRange &particles) {
  auto const pref = dipole.prefactor * 4 * Utils::pi() / box_geo.volume() /
                    (2 * dp3m.params.epsilon + 1);
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

  double energy = 0.0;
  if (energy_flag) {
    double sum_e = 0.0;
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
      p.f.torque[0] -= pref * sumix[ip];
      p.f.torque[1] -= pref * sumiy[ip];
      p.f.torque[2] -= pref * sumiz[ip];
      ip++;
    }
  }

  return energy;
}

/*****************************************************************************/

void dp3m_calc_influence_function_force() {
  auto const start = Utils::Vector3i{dp3m.fft.plan[3].start};
  auto const size = Utils::Vector3i{dp3m.fft.plan[3].new_mesh};

  dp3m.g_force = grid_influence_function<3>(dp3m.params, start, start + size,
                                            box_geo.length());
}

void dp3m_calc_influence_function_energy() {
  auto const start = Utils::Vector3i{dp3m.fft.plan[3].start};
  auto const size = Utils::Vector3i{dp3m.fft.plan[3].new_mesh};

  dp3m.g_energy = grid_influence_function<2>(dp3m.params, start, start + size,
                                             box_geo.length());
}

/*****************************************************************************/

/** @copybrief p3m_get_accuracy
 *
 *  The real space error is tuned such that it contributes half of the
 *  total error, and then the Fourier space error is calculated.
 *  If an optimal alpha is not found, the value 0.1 is used as fallback.
 *  @param[in]  mesh       @copybrief P3MParameters::mesh
 *  @param[in]  cao        @copybrief P3MParameters::cao
 *  @param[in]  r_cut_iL   @copybrief P3MParameters::r_cut_iL
 *  @param[out] _alpha_L   @copybrief P3MParameters::alpha_L
 *  @param[out] _rs_err    real space error
 *  @param[out] _ks_err    Fourier space error
 *  @returns Error magnitude
 */
double dp3m_get_accuracy(int mesh, int cao, double r_cut_iL, double *_alpha_L,
                         double *_rs_err, double *_ks_err) {
  double rs_err, ks_err;
  double alpha_L;

  /* calc maximal real space error for setting */

  // Alpha cannot be zero in the dipolar case because real_space formula breaks
  // down
  rs_err =
      dp3m_real_space_error(box_geo.length()[0], dipole.prefactor, r_cut_iL,
                            dp3m.sum_dip_part, dp3m.sum_mu2, 0.001);

  if (Utils::sqrt_2() * rs_err > dp3m.params.accuracy) {
    /* assume rs_err = ks_err -> rs_err = accuracy/sqrt(2.0) -> alpha_L */
    alpha_L = dp3m_rtbisection(
        box_geo.length()[0], dipole.prefactor, r_cut_iL, dp3m.sum_dip_part,
        dp3m.sum_mu2, 0.0001 * box_geo.length()[0], 5.0 * box_geo.length()[0],
        0.0001, dp3m.params.accuracy);
    if (alpha_L == -DP3M_RTBISECTION_ERROR) {
      *_rs_err = -1;
      *_ks_err = -1;
      return -DP3M_RTBISECTION_ERROR;
    }
  }

  else
    /* even alpha=0 is ok, however, we cannot choose it since it kills the
       k-space error formula.
       Anyways, this very likely NOT the optimal solution */
    alpha_L = 0.1;

  *_alpha_L = alpha_L;
  /* calculate real space and k-space error for this alpha_L */

  rs_err =
      dp3m_real_space_error(box_geo.length()[0], dipole.prefactor, r_cut_iL,
                            dp3m.sum_dip_part, dp3m.sum_mu2, alpha_L);
  ks_err = dp3m_k_space_error(box_geo.length()[0], dipole.prefactor, mesh, cao,
                              dp3m.sum_dip_part, dp3m.sum_mu2, alpha_L);

  *_rs_err = rs_err;
  *_ks_err = ks_err;
  return sqrt(Utils::sqr(rs_err) + Utils::sqr(ks_err));
}

/** @copybrief p3m_mcr_time
 *
 *  @param[in]  mesh            @copybrief P3MParameters::mesh
 *  @param[in]  cao             @copybrief P3MParameters::cao
 *  @param[in]  r_cut_iL        @copybrief P3MParameters::r_cut_iL
 *  @param[in]  alpha_L         @copybrief P3MParameters::alpha_L
 *
 *  @returns The integration time in case of success, otherwise
 *           -@ref P3M_TUNE_FAIL
 */
static double dp3m_mcr_time(int mesh, int cao, double r_cut_iL,
                            double alpha_L) {
  /* rounded up 2000/n_charges timing force evaluations */
  int int_num = (1999 + dp3m.sum_dip_part) / dp3m.sum_dip_part;

  /* broadcast p3m parameters for test run */
  if (dipole.method != DIPOLAR_P3M && dipole.method != DIPOLAR_MDLC_P3M)
    Dipole::set_method_local(DIPOLAR_P3M);
  dp3m.params.r_cut_iL = r_cut_iL;
  dp3m.params.mesh[0] = dp3m.params.mesh[1] = dp3m.params.mesh[2] = mesh;
  dp3m.params.cao = cao;
  dp3m.params.alpha_L = alpha_L;
  /* initialize p3m structures */
  mpi_bcast_coulomb_params();
  /* perform force calculation test */
  double int_time = time_force_calc(int_num);
  if (int_time == -1) {
    return -P3M_TUNE_FAIL;
  }
  return int_time;
}

/** @copybrief p3m_mc_time
 *
 *  The @p _r_cut_iL is determined via a simple bisection.
 *
 *  @param[in]  mesh            @copybrief P3MParameters::mesh
 *  @param[in]  cao             @copybrief P3MParameters::cao
 *  @param[in]  r_cut_iL_min    lower bound for @p _r_cut_iL
 *  @param[in]  r_cut_iL_max    upper bound for @p _r_cut_iL
 *  @param[out] _r_cut_iL       @copybrief P3MParameters::r_cut_iL
 *  @param[out] _alpha_L        @copybrief P3MParameters::alpha_L
 *  @param[out] _accuracy       @copybrief P3MParameters::accuracy
 *  @param[in]  verbose         printf output
 *
 *  @returns The integration time in case of success, otherwise
 *           -@ref P3M_TUNE_FAIL, -@ref P3M_TUNE_ACCURACY_TOO_LARGE,
 *           -@ref P3M_TUNE_CAO_TOO_LARGE, -@ref P3M_TUNE_ELCTEST, or
 *           -@ref P3M_TUNE_CUTOFF_TOO_LARGE
 */
static double dp3m_mc_time(int mesh, int cao, double r_cut_iL_min,
                           double r_cut_iL_max, double *_r_cut_iL,
                           double *_alpha_L, double *_accuracy, bool verbose) {
  double r_cut_iL;
  double rs_err, ks_err;

  /* initial checks. */
  auto const mesh_size = box_geo.length()[0] / static_cast<double>(mesh);
  auto const k_cut = mesh_size * cao / 2.0;

  auto const min_box_l = *boost::min_element(box_geo.length());
  auto const min_local_box_l = *boost::min_element(local_geo.length());

  if (cao >= mesh || k_cut >= std::min(min_box_l, min_local_box_l) - skin) {
    /* print result */
    if (verbose) {
      std::printf("%-4d %-3d  cao too large for this mesh\n", mesh, cao);
    }
    return -P3M_TUNE_CAO_TOO_LARGE;
  }

  /* Either low and high boundary are equal (for fixed cut), or the low border
     is initially 0 and therefore
     has infinite error estimate, as required. Therefore if the high boundary
     fails, there is no possible r_cut */
  *_accuracy =
      dp3m_get_accuracy(mesh, cao, r_cut_iL_max, _alpha_L, &rs_err, &ks_err);
  if (*_accuracy == -DP3M_RTBISECTION_ERROR)
    return *_accuracy;
  if (*_accuracy > dp3m.params.accuracy) {
    /* print result */
    if (verbose) {
      std::printf("%-4d %-3d %.5e %.5e %.5e %.3e %.3e accuracy not achieved\n",
                  mesh, cao, r_cut_iL_max, *_alpha_L, *_accuracy, rs_err,
                  ks_err);
    }
    return -P3M_TUNE_ACCURACY_TOO_LARGE;
  }

  for (;;) {
    r_cut_iL = 0.5 * (r_cut_iL_min + r_cut_iL_max);

    if (r_cut_iL_max - r_cut_iL_min < P3M_RCUT_PREC)
      break;

    /* bisection */
    auto const tmp_accuracy =
        dp3m_get_accuracy(mesh, cao, r_cut_iL, _alpha_L, &rs_err, &ks_err);
    if (tmp_accuracy == -DP3M_RTBISECTION_ERROR)
      return tmp_accuracy;
    if (tmp_accuracy > dp3m.params.accuracy)
      r_cut_iL_min = r_cut_iL;
    else
      r_cut_iL_max = r_cut_iL;
  }
  /* final result is always the upper interval boundary, since only there
     we know that the desired minimal accuracy is obtained */
  *_r_cut_iL = r_cut_iL = r_cut_iL_max;

  /* check whether we are running P3M+DLC, and whether we leave a reasonable gap
   * space */
  if (dipole.method == DIPOLAR_MDLC_P3M) {
    runtimeErrorMsg() << "dipolar P3M: tuning when dlc needs to be fixed";
  }

  double const int_time = dp3m_mcr_time(mesh, cao, r_cut_iL, *_alpha_L);
  if (int_time == -P3M_TUNE_FAIL) {
    if (verbose) {
      std::printf("tuning failed, test integration not possible\n");
    }
    return int_time;
  }

  *_accuracy =
      dp3m_get_accuracy(mesh, cao, r_cut_iL, _alpha_L, &rs_err, &ks_err);
  if (*_accuracy == -DP3M_RTBISECTION_ERROR) {
    return *_accuracy;
  }

  /* print result */
  if (verbose) {
    std::printf("%-4d %-3d %.5e %.5e %.5e %.3e %.3e %-8.0f\n", mesh, cao,
                r_cut_iL, *_alpha_L, *_accuracy, rs_err, ks_err, int_time);
  }
  return int_time;
}

/** @copybrief p3m_m_time
 *
 *  @p _cao should contain an initial guess, which is then adapted by stepping
 *  up and down.
 *
 *  @param[in]      mesh            @copybrief P3MParameters::mesh
 *  @param[in]      cao_min         lower bound for @p _cao
 *  @param[in]      cao_max         upper bound for @p _cao
 *  @param[in,out]  _cao            initial guess for the
 *                                  @copybrief P3MParameters::cao
 *  @param[in]      r_cut_iL_min    lower bound for @p _r_cut_iL
 *  @param[in]      r_cut_iL_max    upper bound for @p _r_cut_iL
 *  @param[out]     _r_cut_iL       @copybrief P3MParameters::r_cut_iL
 *  @param[out]     _alpha_L        @copybrief P3MParameters::alpha_L
 *  @param[out]     _accuracy       @copybrief P3MParameters::accuracy
 *  @param[in]      verbose         printf output
 *
 *  @returns The integration time in case of success, otherwise
 *           -@ref P3M_TUNE_FAIL or -@ref P3M_TUNE_CAO_TOO_LARGE */
static double dp3m_m_time(int mesh, int cao_min, int cao_max, int *_cao,
                          double r_cut_iL_min, double r_cut_iL_max,
                          double *_r_cut_iL, double *_alpha_L,
                          double *_accuracy, bool verbose) {
  double best_time = -1, tmp_r_cut_iL = -1., tmp_alpha_L = 0.0,
         tmp_accuracy = 0.0;
  /* in which direction improvement is possible. Initially, we don't know it
   * yet.
   */
  int final_dir = 0;
  int cao = *_cao;

  /* the initial step sets a timing mark. If there is no valid r_cut, we can
     only try
     to increase cao to increase the obtainable precision of the far formula. */
  double tmp_time;
  do {
    tmp_time =
        dp3m_mc_time(mesh, cao, r_cut_iL_min, r_cut_iL_max, &tmp_r_cut_iL,
                     &tmp_alpha_L, &tmp_accuracy, verbose);
    /* bail out if the force evaluation is not working */
    if (tmp_time == -P3M_TUNE_FAIL || tmp_time == -DP3M_RTBISECTION_ERROR)
      return tmp_time;
    /* cao is too large for this grid, but still the accuracy cannot be
     * achieved, give up */
    if (tmp_time == -P3M_TUNE_CAO_TOO_LARGE) {
      return tmp_time;
    }
    /* we have a valid time, start optimising from there */
    if (tmp_time >= 0) {
      best_time = tmp_time;
      *_r_cut_iL = tmp_r_cut_iL;
      *_alpha_L = tmp_alpha_L;
      *_accuracy = tmp_accuracy;
      *_cao = cao;
      break;
    }
    /* the required accuracy could not be obtained, try higher caos */
    cao++;
    final_dir = 1;
  } while (cao <= cao_max);
  /* with this mesh, the required accuracy cannot be obtained. */
  if (cao > cao_max)
    return -P3M_TUNE_CAO_TOO_LARGE;

  /* at the boundaries, only the opposite direction can be used for optimisation
   */
  if (cao == cao_min)
    final_dir = 1;
  else if (cao == cao_max)
    final_dir = -1;

  if (final_dir == 0) {
    /* check in which direction we can optimise. Both directions are possible */
    double dir_times[3];
    for (final_dir = -1; final_dir <= 1; final_dir += 2) {
      dir_times[final_dir + 1] = tmp_time =
          dp3m_mc_time(mesh, cao + final_dir, r_cut_iL_min, r_cut_iL_max,
                       &tmp_r_cut_iL, &tmp_alpha_L, &tmp_accuracy, verbose);
      /* bail out on errors, as usual */
      if (tmp_time == -P3M_TUNE_FAIL || tmp_time == -DP3M_RTBISECTION_ERROR)
        return tmp_time;
      /* in this direction, we cannot optimise, since we get into precision
       * trouble */
      if (tmp_time < 0)
        continue;

      if (tmp_time < best_time) {
        best_time = tmp_time;
        *_r_cut_iL = tmp_r_cut_iL;
        *_alpha_L = tmp_alpha_L;
        *_accuracy = tmp_accuracy;
        *_cao = cao + final_dir;
      }
    }
    /* choose the direction which was optimal, if any of the two */
    if (dir_times[0] == best_time) {
      final_dir = -1;
    } else if (dir_times[2] == best_time) {
      final_dir = 1;
    } else {
      /* no improvement in either direction, however if one is only marginally
       * worse, we can still try*/
      /* down is possible and not much worse, while up is either illegal or even
       * worse */
      if ((dir_times[0] >= 0 && dir_times[0] < best_time + P3M_TIME_GRAN) &&
          (dir_times[2] < 0 || dir_times[2] > dir_times[0]))
        final_dir = -1;
      /* same for up */
      else if ((dir_times[2] >= 0 &&
                dir_times[2] < best_time + P3M_TIME_GRAN) &&
               (dir_times[0] < 0 || dir_times[0] > dir_times[2]))
        final_dir = 1;
      else {
        /* really no chance for optimisation */
        return best_time;
      }
    }
    /* we already checked the initial cao and its neighbor */
    cao += 2 * final_dir;
  } else {
    /* here some constraint is active, and we only checked the initial cao
     * itself */
    cao += final_dir;
  }

  /* move cao into the optimisation direction until we do not gain anymore. */
  for (; cao >= cao_min && cao <= cao_max; cao += final_dir) {
    tmp_time =
        dp3m_mc_time(mesh, cao, r_cut_iL_min, r_cut_iL_max, &tmp_r_cut_iL,
                     &tmp_alpha_L, &tmp_accuracy, verbose);
    /* bail out on errors, as usual */
    if (tmp_time == -P3M_TUNE_FAIL || tmp_time == -DP3M_RTBISECTION_ERROR)
      return tmp_time;
    /* if we cannot meet the precision anymore, give up */
    if (tmp_time < 0)
      break;

    if (tmp_time < best_time) {
      best_time = tmp_time;
      *_r_cut_iL = tmp_r_cut_iL;
      *_alpha_L = tmp_alpha_L;
      *_accuracy = tmp_accuracy;
      *_cao = cao;
    }
    /* no hope of further optimisation */
    else if (tmp_time > best_time + P3M_TIME_GRAN)
      break;
  }
  return best_time;
}

int dp3m_adaptive_tune(bool verbose) {
  /** Tuning of dipolar P3M. The algorithm basically determines the mesh, cao
   *  and then the real space cutoff, in this nested order.
   *
   *  For each mesh, the cao optimal for the mesh tested previously is used as
   *  an initial guess, and the algorithm tries whether increasing or decreasing
   *  it leads to a better solution. This is efficient, since the optimal cao
   *  only changes little with the meshes in general.
   *
   *  The real space cutoff for a given mesh and cao is determined via a
   *  bisection on the error estimate, which determines where the error
   *  estimate equals the required accuracy. Therefore the smallest possible,
   *  i.e. fastest real space cutoff is determined.
   *
   *  Both the search over mesh and cao stop to search in a specific direction
   *  once the computation time is significantly higher than the currently
   *  known optimum.
   */
  int mesh_max, mesh = -1, tmp_mesh;
  double r_cut_iL_min, r_cut_iL_max, r_cut_iL = -1, tmp_r_cut_iL = 0.0;
  int cao_min, cao_max, cao = -1, tmp_cao;

  double alpha_L = -1, tmp_alpha_L = 0.0;
  double accuracy = -1, tmp_accuracy = 0.0;
  double time_best = 1e20, tmp_time;

  /* preparation */
  mpi_call_all(dp3m_count_magnetic_particles);

  if (dp3m.sum_dip_part == 0) {
    runtimeErrorMsg() << "no dipolar particles in the system";
    return ES_ERROR;
  }

  /* Print Status */
  if (verbose) {
    std::printf("Dipolar P3M tune parameters: Accuracy goal = %.5e prefactor "
                "= %.5e\n"
                "System: box_l = %.5e # charged part = %d Sum[q_i^2] = %.5e\n",
                dp3m.params.accuracy, dipole.prefactor, box_geo.length()[0],
                dp3m.sum_dip_part, dp3m.sum_mu2);
  }

  /* parameter ranges */
  if (dp3m.params.mesh[0] == 0) {
    double expo;
    expo = log(pow((double)dp3m.sum_dip_part, (1.0 / 3.0))) / log(2.0);

    tmp_mesh = (int)(pow(2.0, (double)((int)expo)) + 0.1);
    /* this limits the tried meshes if the accuracy cannot
       be obtained with smaller meshes, but normally not all these
       meshes have to be tested */
    mesh_max = tmp_mesh * 256;
    /* avoid using more than 1 GB of FFT arrays (per default, see config.hpp) */
    if (mesh_max > P3M_MAX_MESH)
      mesh_max = P3M_MAX_MESH;
  } else {
    tmp_mesh = mesh_max = dp3m.params.mesh[0];

    if (verbose) {
      std::printf("fixed mesh %d\n", dp3m.params.mesh[0]);
    }
  }

  if (dp3m.params.r_cut_iL == 0.0) {
    auto const min_box_l = *boost::min_element(box_geo.length());
    auto const min_local_box_l = *boost::min_element(local_geo.length());

    r_cut_iL_min = 0;
    r_cut_iL_max = std::min(min_local_box_l, min_box_l / 2) - skin;
    r_cut_iL_min *= 1. / box_geo.length()[0];
    r_cut_iL_max *= 1. / box_geo.length()[0];
  } else {
    r_cut_iL_min = r_cut_iL_max = dp3m.params.r_cut_iL;

    if (verbose) {
      std::printf("fixed r_cut_iL %f\n", dp3m.params.r_cut_iL);
    }
  }

  if (dp3m.params.cao == 0) {
    cao_min = 1;
    cao_max = 7;
    cao = 3;
  } else {
    cao_min = cao_max = cao = dp3m.params.cao;

    if (verbose) {
      std::printf("fixed cao %d\n", dp3m.params.cao);
    }
  }

  if (verbose) {
    std::printf("Dmesh cao Dr_cut_iL   Dalpha_L     Derr     "
                "    Drs_err    Dks_err    time [ms]\n");
  }

  /* mesh loop */
  for (; tmp_mesh <= mesh_max; tmp_mesh += 2) {
    tmp_cao = cao;
    tmp_time = dp3m_m_time(tmp_mesh, cao_min, cao_max, &tmp_cao, r_cut_iL_min,
                           r_cut_iL_max, &tmp_r_cut_iL, &tmp_alpha_L,
                           &tmp_accuracy, verbose);
    /* some error occurred during the tuning force evaluation */
    if (tmp_time == -P3M_TUNE_FAIL || tmp_time == -DP3M_RTBISECTION_ERROR)
      return ES_ERROR;
    /* this mesh does not work at all */
    if (tmp_time < 0)
      continue;

    /* the optimum r_cut for this mesh is the upper limit for higher meshes,
       everything else is slower */
    r_cut_iL_max = tmp_r_cut_iL;

    /* new optimum */
    if (tmp_time < time_best) {
      time_best = tmp_time;
      mesh = tmp_mesh;
      cao = tmp_cao;
      r_cut_iL = tmp_r_cut_iL;
      alpha_L = tmp_alpha_L;
      accuracy = tmp_accuracy;
    }
    /* no hope of further optimisation */
    else if (tmp_time > time_best + P3M_TIME_GRAN)
      break;
  }

  if (time_best == 1e20) {
    runtimeErrorMsg() << "failed to reach requested accuracy";
    return ES_ERROR;
  }

  /* set tuned p3m parameters */
  dp3m.params.r_cut_iL = r_cut_iL;
  dp3m.params.mesh[0] = dp3m.params.mesh[1] = dp3m.params.mesh[2] = mesh;
  dp3m.params.cao = cao;
  dp3m.params.alpha_L = alpha_L;
  dp3m.params.accuracy = accuracy;
  /* broadcast tuned p3m parameters */
  mpi_bcast_coulomb_params();
  /* Tell the user about the outcome */
  if (verbose) {
    std::printf(
        "\nresulting parameters: mesh: %d, cao: %d, r_cut_iL: %.4e,"
        "\n                      alpha_L: %.4e, accuracy: %.4e, time: %.0f\n",
        mesh, cao, r_cut_iL, alpha_L, accuracy, time_best);
  }
  return ES_OK;
}

void dp3m_count_magnetic_particles() {
  using boost::mpi::all_reduce;
  int local_n = 0;
  double local_mu2 = 0.0;

  for (auto const &p : cell_structure.local_particles()) {
    if (p.p.dipm != 0.0) {
      local_mu2 += p.calc_dip().norm2();
      local_n++;
    }
  }

  dp3m.sum_mu2 = all_reduce(comm_cart, local_mu2, std::plus<>());
  dp3m.sum_dip_part = all_reduce(comm_cart, local_n, std::plus<>());
}

REGISTER_CALLBACK(dp3m_count_magnetic_particles)

/** Calculate the k-space error of dipolar-P3M */
static double dp3m_k_space_error(double box_size, double prefac, int mesh,
                                 int cao, int n_c_part, double sum_q2,
                                 double alpha_L) {
  double he_q = 0.0;
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

  *alias1 = *alias2 = 0.0;
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

        auto const U2 = pow(sinc(fnmx) * sinc(fnmy) * sinc(fnmz), 2.0 * cao);

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
double dp3m_real_space_error(double box_size, double prefac, double r_cut_iL,
                             int n_c_part, double sum_q2, double alpha_L) {
  double d_error_f, d_cc, d_dc, d_rcut2, d_con;
  double d_a2, d_c, d_RCUT;

  d_RCUT = r_cut_iL * box_size;
  d_rcut2 = d_RCUT * d_RCUT;

  d_a2 = Utils::sqr(alpha_L) / Utils::sqr(box_size);

  d_c = sum_q2 * exp(-d_a2 * Utils::sqr(d_RCUT));

  d_cc =
      4.0 * Utils::sqr(d_a2) * Utils::sqr(d_rcut2) + 6.0 * d_a2 * d_rcut2 + 3.0;

  d_dc = 8.0 * Utils::int_pow<3>(d_a2) * Utils::int_pow<3>(d_rcut2) +
         20.0 * Utils::sqr(d_a2) * Utils::sqr(d_rcut2) + 30.0 * d_a2 * d_rcut2 +
         15.0;

  d_con = 1.0 / sqrt(Utils::int_pow<3>(box_size) * Utils::sqr(d_a2) *
                     Utils::int_pow<4>(d_rcut2) * d_RCUT * (double)n_c_part);

  d_error_f = d_c * d_con *
              sqrt((13. / 6.) * Utils::sqr(d_cc) +
                   (2. / 15.) * Utils::sqr(d_dc) - (13. / 15.) * d_cc * d_dc);

  return d_error_f;
}

/** Using bisection, find the root of a function "func-tuned_accuracy/sqrt(2.)"
 *  known to lie between x1 and x2. The root, returned as rtbis, will be
 *  refined until its accuracy is \f$\pm\f$ @p xacc.
 */
double dp3m_rtbisection(double box_size, double prefac, double r_cut_iL,
                        int n_c_part, double sum_q2, double x1, double x2,
                        double xacc, double tuned_accuracy) {
  constexpr int JJ_RTBIS_MAX = 40;

  auto const constant = tuned_accuracy / Utils::sqrt_2();

  auto const f1 =
      dp3m_real_space_error(box_size, prefac, r_cut_iL, n_c_part, sum_q2, x1) -
      constant;
  auto const f2 =
      dp3m_real_space_error(box_size, prefac, r_cut_iL, n_c_part, sum_q2, x2) -
      constant;
  if (f1 * f2 >= 0.0) {
    runtimeErrorMsg()
        << "Root must be bracketed for bisection in dp3m_rtbisection";
    return -DP3M_RTBISECTION_ERROR;
  }
  // Orient the search dx, and set rtb to x1 or x2 ...
  double dx;
  double rtb = f1 < 0.0 ? (dx = x2 - x1, x1) : (dx = x1 - x2, x2);
  for (int j = 1; j <= JJ_RTBIS_MAX; j++) {
    auto const xmid = rtb + (dx *= 0.5);
    auto const fmid = dp3m_real_space_error(box_size, prefac, r_cut_iL,
                                            n_c_part, sum_q2, xmid) -
                      constant;
    if (fmid <= 0.0)
      rtb = xmid;
    if (fabs(dx) < xacc || fmid == 0.0)
      return rtb;
  }
  runtimeErrorMsg() << "Too many bisections in dp3m_rtbisection";
  return -DP3M_RTBISECTION_ERROR;
}

/************************************************************/

void dp3m_init_a_ai_cao_cut() {
  for (int i = 0; i < 3; i++) {
    dp3m.params.ai[i] = (double)dp3m.params.mesh[i] / box_geo.length()[i];
    dp3m.params.a[i] = 1.0 / dp3m.params.ai[i];
    dp3m.params.cao_cut[i] = 0.5 * dp3m.params.a[i] * dp3m.params.cao;
  }
}

/*****************************************************************************/

bool dp3m_sanity_checks_boxl() {
  bool ret = false;
  for (int i = 0; i < 3; i++) {
    /* check k-space cutoff */
    if (dp3m.params.cao_cut[i] >= 0.5 * box_geo.length()[i]) {
      runtimeErrorMsg() << "dipolar P3M_init: k-space cutoff "
                        << dp3m.params.cao_cut[i]
                        << " is larger than half of box dimension "
                        << box_geo.length()[i];
      ret = true;
    }
    if (dp3m.params.cao_cut[i] >= local_geo.length()[i]) {
      runtimeErrorMsg() << "dipolar P3M_init: k-space cutoff "
                        << dp3m.params.cao_cut[i]
                        << " is larger than local box dimension "
                        << local_geo.length()[i];
      ret = true;
    }
  }
  return ret;
}

/*****************************************************************************/

bool dp3m_sanity_checks(const Utils::Vector3i &grid) {
  bool ret = false;

  if (!box_geo.periodic(0) || !box_geo.periodic(1) || !box_geo.periodic(2)) {
    runtimeErrorMsg() << "dipolar P3M requires periodicity 1 1 1";
    ret = true;
  }

  if (cell_structure.decomposition_type() != CELL_STRUCTURE_DOMDEC) {
    runtimeErrorMsg() << "dipolar P3M at present requires the domain "
                         "decomposition cell system";
    ret = true;
  }

  if ((box_geo.length()[0] != box_geo.length()[1]) ||
      (box_geo.length()[1] != box_geo.length()[2])) {
    runtimeErrorMsg() << "dipolar P3M requires a cubic box";
    ret = true;
  }

  if ((dp3m.params.mesh[0] != dp3m.params.mesh[1]) ||
      (dp3m.params.mesh[1] != dp3m.params.mesh[2])) {
    runtimeErrorMsg() << "dipolar P3M requires a cubic mesh";
    ret = true;
  }

  if (dp3m_sanity_checks_boxl())
    ret = true;

  if (dp3m.params.mesh[0] == 0) {
    runtimeErrorMsg() << "dipolar P3M_init: mesh size is not yet set";
    ret = true;
  }
  if (dp3m.params.cao == 0) {
    runtimeErrorMsg() << "dipolar P3M_init: cao is not yet set";
    ret = true;
  }
  if (grid[0] < grid[1] || grid[1] < grid[2]) {
    runtimeErrorMsg()
        << "dipolar P3M_init: node grid must be sorted, largest first";
    ret = true;
  }

  return ret;
}

/************************************************/

void dp3m_scaleby_box_l() {
  if (dipole.prefactor < 0.0) {
    runtimeErrorMsg() << "Dipolar prefactor has to be >=0";
    return;
  }

  dp3m.params.r_cut = dp3m.params.r_cut_iL * box_geo.length()[0];
  dp3m.params.alpha = dp3m.params.alpha_L / box_geo.length()[0];
  dp3m_init_a_ai_cao_cut();
  p3m_calc_lm_ld_pos(dp3m.local_mesh, dp3m.params);
  dp3m_sanity_checks_boxl();

  dp3m_calc_influence_function_force();
  dp3m_calc_influence_function_energy();
}

/** Calculate the dipolar-P3M energy correction */
void dp3m_compute_constants_energy_dipolar() {

  if (dp3m.energy_correction != 0.0)
    return;

  auto const Ukp3m = dp3m_average_dipolar_self_energy() * box_geo.volume();

  auto const Eself = -2 * pow(dp3m.params.alpha_L, 3) * Utils::sqrt_pi_i() / 3;

  dp3m.energy_correction =
      -dp3m.sum_mu2 * (Ukp3m + Eself + 2. * Utils::pi() / 3.);
}

#ifdef NPT
void npt_add_virial_magnetic_contribution(double energy) {
  npt_add_virial_contribution(energy);
}
#endif

#endif /* DP3M */
