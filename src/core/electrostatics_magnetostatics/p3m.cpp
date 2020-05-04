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
/** @file
 *
 *  The corresponding header file is @ref p3m.hpp.
 */
#include "p3m.hpp"

#ifdef P3M

#include "Particle.hpp"
#include "cells.hpp"
#include "communication.hpp"
#include "domain_decomposition.hpp"
#include "electrostatics_magnetostatics/coulomb.hpp"
#include "electrostatics_magnetostatics/elc.hpp"
#include "electrostatics_magnetostatics/p3m.hpp"
#include "errorhandling.hpp"
#include "event.hpp"
#include "fft.hpp"
#include "global.hpp"
#include "grid.hpp"
#include "integrate.hpp"
#include "tuning.hpp"
#ifdef CUDA
#include "p3m_gpu_error.hpp"
#endif

#include <utils/math/int_pow.hpp>
#include <utils/math/sinc.hpp>
using Utils::sinc;
#include <utils/strcat_alloc.hpp>
using Utils::strcat_alloc;
#include <utils/constants.hpp>
#include <utils/integral_parameter.hpp>
#include <utils/math/sqr.hpp>

#include <boost/optional.hpp>
#include <boost/range/algorithm/min_element.hpp>
#include <boost/range/numeric.hpp>
#include <mpi.h>

#include <complex>
#include <cstdio>
#include <cstring>

/************************************************
 * variables
 ************************************************/
p3m_data_struct p3m;

/* @name Index helpers for direct and reciprocal space
 * After the FFT the data is in order YZX, which
 * means that Y is the slowest changing index.
 * The defines are here to not get confused and
 * be able to easily change the order.
 */
/*@{*/
#define RX 0
#define RY 1
#define RZ 2
#define KY 0
#define KZ 1
#define KX 2
/*@}*/

/** \name Private Functions */
/*@{*/

/** Initialize the (inverse) mesh constant @ref P3MParameters::a "a"
 *  (@ref P3MParameters::ai "ai") and the cutoff for charge assignment
 *  @ref P3MParameters::cao_cut "cao_cut".
 *
 *  Function called by @ref p3m_init() once and by @ref p3m_scaleby_box_l()
 *  whenever the box length changes.
 */
static void p3m_init_a_ai_cao_cut();

static bool p3m_sanity_checks_system(const Utils::Vector3i &grid);

/** Checks for correctness for charges in P3M of the cao_cut,
 *  necessary when the box length changes
 */
static bool p3m_sanity_checks_boxl();

/** Shift the mesh points by mesh/2 */
static void p3m_calc_meshift();

/** Calculate the Fourier transformed differential operator.
 *  Remark: This is done on the level of n-vectors and not k-vectors,
 *          i.e. the prefactor i*2*PI/L is missing!
 */
static void p3m_calc_differential_operator();

/** Calculate the optimal influence function of @cite hockney88a.
 *  (optimised for force calculations)
 *
 *  Each node calculates only the values for its domain in k-space
 *  (see fft.plan[3].mesh and fft.plan[3].start).
 *
 *  See also: @cite hockney88a 8-22 (p275). Note the somewhat
 *  different convention for the prefactors, which is described in
 *  @cite deserno98a @cite deserno98b.
 */
static void p3m_calc_influence_function_force();

/** Calculate the influence function optimized for the energy and the
 *  self energy correction.
 */
static void p3m_calc_influence_function_energy();

/*@}*/

/** @name P3M tuning helper functions */
/*@{*/

/** Calculate the real space contribution to the rms error in the force (as
 *  described by Kolafa and Perram).
 *  \param prefac     Prefactor of Coulomb interaction.
 *  \param r_cut_iL   rescaled real space cutoff for p3m method.
 *  \param n_c_part   number of charged particles in the system.
 *  \param sum_q2     sum of square of charges in the system
 *  \param alpha_L    rescaled Ewald splitting parameter.
 *  \return real space error
 */
static double p3m_real_space_error(double prefac, double r_cut_iL, int n_c_part,
                                   double sum_q2, double alpha_L);

/** Calculate the analytic expression of the error estimate for the
 *  P3M method in @cite hockney88a (eq. (8.23)) in
 *  order to obtain the rms error in the force for a system of N
 *  randomly distributed particles in a cubic box (k-space part).
 *  \param prefac   Prefactor of Coulomb interaction.
 *  \param mesh     number of mesh points in one direction.
 *  \param cao      charge assignment order.
 *  \param n_c_part number of charged particles in the system.
 *  \param sum_q2   sum of square of charges in the system
 *  \param alpha_L  rescaled Ewald splitting parameter.
 *  \return reciprocal (k) space error
 */
static double p3m_k_space_error(double prefac, const int mesh[3], int cao,
                                int n_c_part, double sum_q2, double alpha_L);

/** Aliasing sum used by \ref p3m_k_space_error. */
static void p3m_tune_aliasing_sums(int nx, int ny, int nz, const int mesh[3],
                                   const double mesh_i[3], int cao,
                                   double alpha_L_i, double *alias1,
                                   double *alias2);

p3m_data_struct::p3m_data_struct() {
  /* local_mesh is uninitialized */
  /* sm is uninitialized */

  sum_qpart = 0;
  sum_q2 = 0.0;
  square_sum_q = 0.0;

  ks_pnum = 0;
}

void p3m_init() {
  if (coulomb.prefactor <= 0.0) {
    // prefactor is zero: electrostatics switched off
    p3m.params.r_cut = 0.0;
    p3m.params.r_cut_iL = 0.0;
  } else {
    if (p3m_sanity_checks()) {
      return;
    }
    p3m.params.cao3 = p3m.params.cao * p3m.params.cao * p3m.params.cao;

    /* initializes the (inverse) mesh constant p3m.params.a (p3m.params.ai) and
     * the cutoff for charge assignment p3m.params.cao_cut */
    p3m_init_a_ai_cao_cut();

    p3m_calc_local_ca_mesh(p3m.local_mesh, p3m.params, local_geo, skin);

    p3m.sm.resize(comm_cart, p3m.local_mesh);

    int ca_mesh_size = fft_init(p3m.local_mesh.dim, p3m.local_mesh.margin,
                                p3m.params.mesh, p3m.params.mesh_off,
                                &p3m.ks_pnum, p3m.fft, node_grid, comm_cart);
    p3m.rs_mesh.resize(ca_mesh_size);
    for (auto &e : p3m.E_mesh) {
      e.resize(ca_mesh_size);
    }

    /* k-space part: */
    p3m_calc_differential_operator();

    /* fix box length dependent constants */
    p3m_scaleby_box_l();

    p3m_count_charged_particles();
  }
}

void p3m_set_tune_params(double r_cut, const int mesh[3], int cao, double alpha,
                         double accuracy) {
  if (r_cut >= 0) {
    p3m.params.r_cut = r_cut;
    p3m.params.r_cut_iL = r_cut * (1. / box_geo.length()[0]);
  }

  if (mesh[0] >= 0) {
    p3m.params.mesh[0] = mesh[0];
    p3m.params.mesh[1] = mesh[1];
    p3m.params.mesh[2] = mesh[2];
  }

  if (cao >= 0)
    p3m.params.cao = cao;

  if (alpha >= 0) {
    p3m.params.alpha = alpha;
    p3m.params.alpha_L = alpha * box_geo.length()[0];
  }

  if (accuracy >= 0)
    p3m.params.accuracy = accuracy;
}

/*@}*/

int p3m_set_params(double r_cut, const int *mesh, int cao, double alpha,
                   double accuracy) {
  if (coulomb.method != COULOMB_P3M && coulomb.method != COULOMB_ELC_P3M &&
      coulomb.method != COULOMB_P3M_GPU)
    coulomb.method = COULOMB_P3M;

  if (r_cut < 0)
    return -1;

  if ((mesh[0] < 0) || (mesh[1] < 0) || (mesh[2] < 0))
    return -2;

  if (cao < 1 || cao > 7 || cao > mesh[0] || cao > mesh[1] || cao > mesh[2])
    return -3;

  p3m.params.r_cut = r_cut;
  p3m.params.r_cut_iL = r_cut * (1. / box_geo.length()[0]);
  p3m.params.mesh[2] = mesh[2];
  p3m.params.mesh[1] = mesh[1];
  p3m.params.mesh[0] = mesh[0];
  p3m.params.cao = cao;

  if (alpha > 0) {
    p3m.params.alpha = alpha;
    p3m.params.alpha_L = alpha * box_geo.length()[0];
  } else if (alpha != -1.0)
    return -4;

  if (accuracy >= 0)
    p3m.params.accuracy = accuracy;
  else if (accuracy != -1.0)
    return -5;

  mpi_bcast_coulomb_params();

  return 0;
}

int p3m_set_mesh_offset(double x, double y, double z) {
  if (x < 0.0 || x > 1.0 || y < 0.0 || y > 1.0 || z < 0.0 || z > 1.0)
    return ES_ERROR;

  p3m.params.mesh_off[0] = x;
  p3m.params.mesh_off[1] = y;
  p3m.params.mesh_off[2] = z;

  mpi_bcast_coulomb_params();

  return ES_OK;
}

int p3m_set_eps(double eps) {
  p3m.params.epsilon = eps;

  mpi_bcast_coulomb_params();

  return ES_OK;
}

namespace {
template <size_t cao> struct AssignCharge {
  void operator()(double q, const Utils::Vector3d &real_pos,
                  const Utils::Vector3d &ai, p3m_local_mesh const &local_mesh,
                  p3m_interpolation_cache &inter_weights) {
    auto const w =
        p3m_calculate_interpolation_weights<cao>(real_pos, ai, local_mesh);

    inter_weights.store(w);

    p3m_interpolate(local_mesh, w,
                    [q](int ind, double w) { p3m.rs_mesh[ind] += w * q; });
  }

  void operator()(double q, const Utils::Vector3d &real_pos,
                  const Utils::Vector3d &ai, p3m_local_mesh const &local_mesh) {
    p3m_interpolate(
        local_mesh,
        p3m_calculate_interpolation_weights<cao>(real_pos, ai, local_mesh),
        [q](int ind, double w) { p3m.rs_mesh[ind] += w * q; });
  }

  void operator()(const ParticleRange &particles) {
    for (auto &p : particles) {
      if (p.p.q != 0.0) {
        this->operator()(p.p.q, p.r.p, p3m.params.ai, p3m.local_mesh,
                         p3m.inter_weights);
      }
    }
  }
};
} // namespace

void p3m_charge_assign(const ParticleRange &particles) {
  p3m.inter_weights.reset(p3m.params.cao);

  /* prepare local FFT mesh */
  for (int i = 0; i < p3m.local_mesh.size; i++)
    p3m.rs_mesh[i] = 0.0;

  Utils::integral_parameter<AssignCharge, 1, 7>(p3m.params.cao, particles);
}

void p3m_assign_charge(double q, const Utils::Vector3d &real_pos,
                       p3m_interpolation_cache &inter_weights) {
  Utils::integral_parameter<AssignCharge, 1, 7>(p3m.params.cao, q, real_pos,
                                                p3m.params.ai, p3m.local_mesh,
                                                inter_weights);
}

void p3m_assign_charge(double q, const Utils::Vector3d &real_pos) {
  Utils::integral_parameter<AssignCharge, 1, 7>(p3m.params.cao, q, real_pos,
                                                p3m.params.ai, p3m.local_mesh);
}

namespace {
template <size_t cao> struct AssignForces {
  void operator()(double force_prefac, const ParticleRange &particles) const {
    using Utils::make_const_span;
    using Utils::Span;
    using Utils::Vector;

    assert(cao == p3m.inter_weights.cao());

    /* charged particle counter, charge fraction counter */
    int cp_cnt = 0;

    for (auto &p : particles) {
      auto const q = p.p.q;
      if (q != 0.0) {
        auto const pref = q * force_prefac;
        auto const w = p3m.inter_weights.load<cao>(cp_cnt++);

        Utils::Vector3d E{};
        p3m_interpolate(p3m.local_mesh, w, [&E](int ind, double w) {
          E += w * Utils::Vector3d{p3m.E_mesh[0][ind], p3m.E_mesh[1][ind],
                                   p3m.E_mesh[2][ind]};
        });

        p.f.f -= pref * E;
      }
    }
  }
};

auto dipole_moment(Particle const &p, BoxGeometry const &box) {
  return p.p.q * unfolded_position(p.r.p, p.l.i, box.length());
}

auto calc_dipole_moment(boost::mpi::communicator const &comm,
                        const ParticleRange &particles,
                        BoxGeometry const &box) {
  auto const local_dip = boost::accumulate(
      particles, Utils::Vector3d{}, [&box](Utils::Vector3d dip, auto const &p) {
        return dip + dipole_moment(p, box);
      });

  return boost::mpi::all_reduce(comm, local_dip, std::plus<>());
}

void add_dipole_correction(Utils::Vector3d const &box_dipole,
                           const ParticleRange &particles) {
  auto const pref = coulomb.prefactor * 4 * M_PI / box_geo.volume() /
                    (2 * p3m.params.epsilon + 1);

  auto const dm = pref * box_dipole;

  for (auto &p : particles) {
    p.f.f -= p.p.q * dm;
  }
}

double dipole_correction_energy(Utils::Vector3d const &box_dipole) {
  auto const pref = coulomb.prefactor * 4 * M_PI / box_geo.volume() /
                    (2 * p3m.params.epsilon + 1);

  return pref * box_dipole.norm2();
}
} // namespace

/** @details Calculate the long range electrostatics part of the stress
 *  tensor. This is part \f$\Pi_{\textrm{dir}, \alpha, \beta}\f$ eq. (2.6)
 *  in @cite essmann95a. The part \f$\Pi_{\textrm{corr}, \alpha, \beta}\f$
 *  eq. (2.8) is not present here since M is the empty set in our simulations.
 */
Utils::Vector9d p3m_calc_kspace_stress() {
  Utils::Vector9d node_k_space_stress{};

  if (p3m.sum_q2 > 0) {
    p3m.sm.gather_grid(p3m.rs_mesh.data(), comm_cart, p3m.local_mesh.dim);
    fft_perform_forw(p3m.rs_mesh.data(), p3m.fft, comm_cart);

    int ind = 0;
    int j[3];
    for (j[0] = 0; j[0] < p3m.fft.plan[3].new_mesh[RX]; j[0]++) {
      for (j[1] = 0; j[1] < p3m.fft.plan[3].new_mesh[RY]; j[1]++) {
        for (j[2] = 0; j[2] < p3m.fft.plan[3].new_mesh[RZ]; j[2]++) {
          auto const kx = 2.0 * Utils::pi() *
                          p3m.d_op[RX][j[KX] + p3m.fft.plan[3].start[KX]] /
                          box_geo.length()[RX];
          auto const ky = 2.0 * Utils::pi() *
                          p3m.d_op[RY][j[KY] + p3m.fft.plan[3].start[KY]] /
                          box_geo.length()[RY];
          auto const kz = 2.0 * Utils::pi() *
                          p3m.d_op[RZ][j[KZ] + p3m.fft.plan[3].start[KZ]] /
                          box_geo.length()[RZ];
          auto const sqk = Utils::sqr(kx) + Utils::sqr(ky) + Utils::sqr(kz);

          auto const node_k_space_energy =
              (sqk == 0)
                  ? 0.0
                  : p3m.g_energy[ind] * (Utils::sqr(p3m.rs_mesh[2 * ind]) +
                                         Utils::sqr(p3m.rs_mesh[2 * ind + 1]));
          ind++;

          auto const vterm =
              (sqk == 0)
                  ? 0.
                  : -2.0 * (1 / sqk + Utils::sqr(1.0 / 2.0 / p3m.params.alpha));

          node_k_space_stress[0] +=
              node_k_space_energy *
              (1.0 + vterm * Utils::sqr(kx)); /* sigma_xx */
          node_k_space_stress[1] +=
              node_k_space_energy * (vterm * kx * ky); /* sigma_xy */
          node_k_space_stress[2] +=
              node_k_space_energy * (vterm * kx * kz); /* sigma_xz */

          node_k_space_stress[3] +=
              node_k_space_energy * (vterm * kx * ky); /* sigma_yx */
          node_k_space_stress[4] +=
              node_k_space_energy *
              (1.0 + vterm * Utils::sqr(ky)); /* sigma_yy */
          node_k_space_stress[5] +=
              node_k_space_energy * (vterm * ky * kz); /* sigma_yz */

          node_k_space_stress[6] +=
              node_k_space_energy * (vterm * kx * kz); /* sigma_zx */
          node_k_space_stress[7] +=
              node_k_space_energy * (vterm * ky * kz); /* sigma_zy */
          node_k_space_stress[8] +=
              node_k_space_energy *
              (1.0 + vterm * Utils::sqr(kz)); /* sigma_zz */
        }
      }
    }
  }

  auto const force_prefac = coulomb.prefactor / (2.0 * box_geo.volume());
  return force_prefac * node_k_space_stress;
}

double p3m_calc_kspace_forces(bool force_flag, bool energy_flag,
                              const ParticleRange &particles) {
  /* Gather information for FFT grid inside the nodes domain (inner local mesh)
   * and perform forward 3D FFT (Charge Assignment Mesh). */
  p3m.sm.gather_grid(p3m.rs_mesh.data(), comm_cart, p3m.local_mesh.dim);
  fft_perform_forw(p3m.rs_mesh.data(), p3m.fft, comm_cart);

  // Note: after these calls, the grids are in the order yzx and not xyz
  // anymore!!!
  /* The dipole moment is only needed if we don't have metallic boundaries. */
  auto const box_dipole = (p3m.params.epsilon != P3M_EPSILON_METALLIC)
                              ? boost::make_optional(calc_dipole_moment(
                                    comm_cart, particles, box_geo))
                              : boost::none;

  /* === k-space calculations === */

  /* === k-space force calculation  === */
  if (force_flag) {
    auto const force_prefac = coulomb.prefactor / (2 * box_geo.volume());

    /* sqrt(-1)*k differentiation */
    int j[3];
    int ind = 0;
    for (j[0] = 0; j[0] < p3m.fft.plan[3].new_mesh[0]; j[0]++) {
      for (j[1] = 0; j[1] < p3m.fft.plan[3].new_mesh[1]; j[1]++) {
        for (j[2] = 0; j[2] < p3m.fft.plan[3].new_mesh[2]; j[2]++) {
          auto const rho_hat = std::complex<double>(p3m.rs_mesh[2 * ind + 0],
                                                    p3m.rs_mesh[2 * ind + 1]);
          auto const phi_hat = p3m.g_force[ind] * rho_hat;

          for (int d = 0; d < 3; d++) {
            /* direction in r-space: */
            int d_rs = (d + p3m.ks_pnum) % 3;
            /* directions */
            auto const k = 2.0 * Utils::pi() *
                           p3m.d_op[d_rs][j[d] + p3m.fft.plan[3].start[d]] /
                           box_geo.length()[d_rs];

            /* i*k*(Re+i*Im) = - Im*k + i*Re*k     (i=sqrt(-1)) */
            p3m.E_mesh[d_rs][2 * ind + 0] = -k * phi_hat.imag();
            p3m.E_mesh[d_rs][2 * ind + 1] = +k * phi_hat.real();
          }

          ind++;
        }
      }
    }

    /* Back FFT force component mesh */
    for (int d = 0; d < 3; d++) {
      fft_perform_back(p3m.E_mesh[d].data(),
                       /* check_complex */ !p3m.params.tuning, p3m.fft,
                       comm_cart);
    }

    {
      std::array<double *, 3> E_fields = {
          p3m.E_mesh[0].data(), p3m.E_mesh[1].data(), p3m.E_mesh[2].data()};
      /* redistribute force component mesh */
      p3m.sm.spread_grid(Utils::make_span(E_fields), comm_cart,
                         p3m.local_mesh.dim);
    }

    Utils::integral_parameter<AssignForces, 1, 7>(p3m.params.cao, force_prefac,
                                                  particles);

    if (p3m.params.epsilon != P3M_EPSILON_METALLIC) {
      add_dipole_correction(box_dipole.value(), particles);
    }
  } /* if(force_flag) */

  /* === k-space energy calculation  === */
  if (energy_flag) {
    double node_k_space_energy = 0.;

    for (int i = 0; i < p3m.fft.plan[3].new_size; i++) {
      // Use the energy optimized influence function for energy!
      node_k_space_energy +=
          p3m.g_energy[i] *
          (Utils::sqr(p3m.rs_mesh[2 * i]) + Utils::sqr(p3m.rs_mesh[2 * i + 1]));
    }
    node_k_space_energy *= coulomb.prefactor / (2 * box_geo.volume());

    double k_space_energy = 0.0;
    MPI_Reduce(&node_k_space_energy, &k_space_energy, 1, MPI_DOUBLE, MPI_SUM, 0,
               comm_cart);
    if (this_node == 0) {
      /* self energy correction */
      k_space_energy -= coulomb.prefactor *
                        (p3m.sum_q2 * p3m.params.alpha * Utils::sqrt_pi_i());
      /* net charge correction */
      k_space_energy -= coulomb.prefactor * p3m.square_sum_q * Utils::pi() /
                        (2.0 * box_geo.volume() * Utils::sqr(p3m.params.alpha));
      /* dipole correction */
      if (p3m.params.epsilon != P3M_EPSILON_METALLIC) {
        k_space_energy += dipole_correction_energy(box_dipole.value());
      }
    }
    return k_space_energy;
  } /* if (energy_flag) */

  return 0.0;
}

void p3m_calc_meshift() {
  p3m.meshift_x.resize(p3m.params.mesh[0]);
  p3m.meshift_y.resize(p3m.params.mesh[1]);
  p3m.meshift_z.resize(p3m.params.mesh[2]);

  p3m.meshift_x[0] = p3m.meshift_y[0] = p3m.meshift_z[0] = 0;
  for (int i = 1; i <= p3m.params.mesh[RX] / 2; i++) {
    p3m.meshift_x[i] = i;
    p3m.meshift_x[p3m.params.mesh[0] - i] = -i;
  }

  for (int i = 1; i <= p3m.params.mesh[RY] / 2; i++) {
    p3m.meshift_y[i] = i;
    p3m.meshift_y[p3m.params.mesh[1] - i] = -i;
  }

  for (int i = 1; i <= p3m.params.mesh[RZ] / 2; i++) {
    p3m.meshift_z[i] = i;
    p3m.meshift_z[p3m.params.mesh[2] - i] = -i;
  }
}

void p3m_calc_differential_operator() {

  for (int i = 0; i < 3; i++) {
    p3m.d_op[i].resize(p3m.params.mesh[i]);
    p3m.d_op[i][0] = 0;
    p3m.d_op[i][p3m.params.mesh[i] / 2] = 0.0;

    for (int j = 1; j < p3m.params.mesh[i] / 2; j++) {
      p3m.d_op[i][j] = j;
      p3m.d_op[i][p3m.params.mesh[i] - j] = -j;
    }
  }
}

namespace {

template <int cao>
inline double perform_aliasing_sums_force(int const n[3], double numerator[3]) {
  using Utils::int_pow;

  int i;
  double denominator = 0.0;
  /* lots of temporary variables... */
  double sx, sy, sz, f1, f2, mx, my, mz, nmx, nmy, nmz, nm2, expo;
  double limit = 30;

  for (i = 0; i < 3; i++)
    numerator[i] = 0.0;

  f1 = Utils::sqr(Utils::pi() / (p3m.params.alpha));

  for (mx = -P3M_BRILLOUIN; mx <= P3M_BRILLOUIN; mx++) {
    nmx = p3m.meshift_x[n[KX]] + p3m.params.mesh[RX] * mx;
    sx = int_pow<2 * cao>(sinc(nmx / (double)p3m.params.mesh[RX]));
    for (my = -P3M_BRILLOUIN; my <= P3M_BRILLOUIN; my++) {
      nmy = p3m.meshift_y[n[KY]] + p3m.params.mesh[RY] * my;
      sy = sx * int_pow<2 * cao>(sinc(nmy / (double)p3m.params.mesh[RY]));
      for (mz = -P3M_BRILLOUIN; mz <= P3M_BRILLOUIN; mz++) {
        nmz = p3m.meshift_z[n[KZ]] + p3m.params.mesh[RZ] * mz;
        sz = sy * int_pow<2 * cao>(sinc(nmz / (double)p3m.params.mesh[RZ]));

        nm2 = Utils::sqr(nmx / box_geo.length()[RX]) +
              Utils::sqr(nmy / box_geo.length()[RY]) +
              Utils::sqr(nmz / box_geo.length()[RZ]);
        expo = f1 * nm2;
        f2 = (expo < limit) ? sz * exp(-expo) / nm2 : 0.0;

        numerator[RX] += f2 * nmx / box_geo.length()[RX];
        numerator[RY] += f2 * nmy / box_geo.length()[RY];
        numerator[RZ] += f2 * nmz / box_geo.length()[RZ];

        denominator += sz;
      }
    }
  }
  return denominator;
}

template <int cao> void calc_influence_function_force() {
  int i, n[3], ind;
  int end[3];
  int size = 1;
  double fak1, fak2, fak3;
  double nominator[3] = {0.0, 0.0, 0.0};

  p3m_calc_meshift();

  for (i = 0; i < 3; i++) {
    size *= p3m.fft.plan[3].new_mesh[i];
    end[i] = p3m.fft.plan[3].start[i] + p3m.fft.plan[3].new_mesh[i];
  }

  p3m.g_force.resize(size);

  /* Skip influence function calculation in tuning mode,
     the results need not be correct for timing. */
  if (p3m.params.tuning) {
    /* If resized, fill with zeros to avoid nan forces. */
    memset(p3m.g_force.data(), 0, size * sizeof(double));

    return;
  }

  for (n[0] = p3m.fft.plan[3].start[0]; n[0] < end[0]; n[0]++) {
    for (n[1] = p3m.fft.plan[3].start[1]; n[1] < end[1]; n[1]++) {
      for (n[2] = p3m.fft.plan[3].start[2]; n[2] < end[2]; n[2]++) {
        ind =
            (n[2] - p3m.fft.plan[3].start[2]) +
            p3m.fft.plan[3].new_mesh[2] * ((n[1] - p3m.fft.plan[3].start[1]) +
                                           (p3m.fft.plan[3].new_mesh[1] *
                                            (n[0] - p3m.fft.plan[3].start[0])));

        if ((n[KX] % (p3m.params.mesh[RX] / 2) == 0) &&
            (n[KY] % (p3m.params.mesh[RY] / 2) == 0) &&
            (n[KZ] % (p3m.params.mesh[RZ] / 2) == 0)) {
          p3m.g_force[ind] = 0.0;
        } else {
          const double denominator =
              perform_aliasing_sums_force<cao>(n, nominator);

          fak1 = p3m.d_op[RX][n[KX]] * nominator[RX] / box_geo.length()[RX] +
                 p3m.d_op[RY][n[KY]] * nominator[RY] / box_geo.length()[RY] +
                 p3m.d_op[RZ][n[KZ]] * nominator[RZ] / box_geo.length()[RZ];
          fak2 = Utils::sqr(p3m.d_op[RX][n[KX]] / box_geo.length()[RX]) +
                 Utils::sqr(p3m.d_op[RY][n[KY]] / box_geo.length()[RY]) +
                 Utils::sqr(p3m.d_op[RZ][n[KZ]] / box_geo.length()[RZ]);

          if (fak2 == 0) {
            fak3 = 0;
          } else
            fak3 = fak1 / (fak2 * Utils::sqr(denominator));

          p3m.g_force[ind] = 2 * fak3 / (Utils::pi());
        }
      }
    }
  }
}

} /* namespace */

void p3m_calc_influence_function_force() {
  switch (p3m.params.cao) {
  case 1:
    calc_influence_function_force<1>();
    break;
  case 2:
    calc_influence_function_force<2>();
    break;
  case 3:
    calc_influence_function_force<3>();
    break;
  case 4:
    calc_influence_function_force<4>();
    break;
  case 5:
    calc_influence_function_force<5>();
    break;
  case 6:
    calc_influence_function_force<6>();
    break;
  case 7:
    calc_influence_function_force<7>();
    break;
  }
}

namespace {

template <int cao> inline double perform_aliasing_sums_energy(int const n[3]) {
  using Utils::int_pow;
  double numerator = 0.0, denominator = 0.0;
  /* lots of temporary variables... */
  double sx, sy, sz, f1, f2, mx, my, mz, nmx, nmy, nmz, nm2, expo;
  double limit = 30;

  f1 = Utils::sqr(Utils::pi() / (p3m.params.alpha));

  for (mx = -P3M_BRILLOUIN; mx <= P3M_BRILLOUIN; mx++) {
    nmx = p3m.meshift_x[n[KX]] + p3m.params.mesh[RX] * mx;
    sx = int_pow<2 * cao>(sinc(nmx / (double)p3m.params.mesh[RX]));
    for (my = -P3M_BRILLOUIN; my <= P3M_BRILLOUIN; my++) {
      nmy = p3m.meshift_y[n[KY]] + p3m.params.mesh[RY] * my;
      sy = sx * int_pow<2 * cao>(sinc(nmy / (double)p3m.params.mesh[RY]));
      for (mz = -P3M_BRILLOUIN; mz <= P3M_BRILLOUIN; mz++) {
        nmz = p3m.meshift_z[n[KZ]] + p3m.params.mesh[RZ] * mz;
        sz = sy * int_pow<2 * cao>(sinc(nmz / (double)p3m.params.mesh[RZ]));
        /* k = 2*pi * (nx/lx, ny/ly, nz/lz); expo = -k^2 / 4*alpha^2 */
        nm2 = Utils::sqr(nmx / box_geo.length()[RX]) +
              Utils::sqr(nmy / box_geo.length()[RY]) +
              Utils::sqr(nmz / box_geo.length()[RZ]);
        expo = f1 * nm2;
        f2 = (expo < limit) ? sz * exp(-expo) / nm2 : 0.0;

        numerator += f2;
        denominator += sz;
      }
    }
  }

  return numerator / Utils::sqr(denominator);
}

template <int cao> void calc_influence_function_energy() {
  int i, n[3], ind;
  int end[3];
  int start[3];
  int size = 1;

  p3m_calc_meshift();

  for (i = 0; i < 3; i++) {
    size *= p3m.fft.plan[3].new_mesh[i];
    end[i] = p3m.fft.plan[3].start[i] + p3m.fft.plan[3].new_mesh[i];
    start[i] = p3m.fft.plan[3].start[i];
  }

  p3m.g_energy.resize(size);

  /* Skip influence function calculation in tuning mode,
     the results need not be correct for timing. */
  if (p3m.params.tuning)
    return;

  ind = 0;

  for (n[0] = start[0]; n[0] < end[0]; n[0]++) {
    for (n[1] = start[1]; n[1] < end[1]; n[1]++) {
      for (n[2] = start[2]; n[2] < end[2]; n[2]++) {
        ind = (n[2] - start[2]) +
              p3m.fft.plan[3].new_mesh[2] * (n[1] - start[1]) +
              p3m.fft.plan[3].new_mesh[2] * p3m.fft.plan[3].new_mesh[1] *
                  (n[0] - start[0]);
        if ((n[KX] % (p3m.params.mesh[RX] / 2) == 0) &&
            (n[KY] % (p3m.params.mesh[RY] / 2) == 0) &&
            (n[KZ] % (p3m.params.mesh[RZ] / 2) == 0)) {
          p3m.g_energy[ind] = 0.0;
        }

        else
          p3m.g_energy[ind] =
              perform_aliasing_sums_energy<cao>(n) / Utils::pi();
      }
    }
  }
}

} /* namespace */

void p3m_calc_influence_function_energy() {
  switch (p3m.params.cao) {
  case 1:
    calc_influence_function_energy<1>();
    break;
  case 2:
    calc_influence_function_energy<2>();
    break;
  case 3:
    calc_influence_function_energy<3>();
    break;
  case 4:
    calc_influence_function_energy<4>();
    break;
  case 5:
    calc_influence_function_energy<5>();
    break;
  case 6:
    calc_influence_function_energy<6>();
    break;
  case 7:
    calc_influence_function_energy<7>();
    break;
  }
}

#define P3M_TUNE_MAX_CUTS 50

/** Get the minimal error for this combination of parameters.
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
static double p3m_get_accuracy(const int mesh[3], int cao, double r_cut_iL,
                               double *_alpha_L, double *_rs_err,
                               double *_ks_err) {
  double rs_err, ks_err;
  double alpha_L;

  /* calc maximal real space error for setting */
  rs_err = p3m_real_space_error(coulomb.prefactor, r_cut_iL, p3m.sum_qpart,
                                p3m.sum_q2, 0);

  if (M_SQRT2 * rs_err > p3m.params.accuracy) {
    /* assume rs_err = ks_err -> rs_err = accuracy/sqrt(2.0) -> alpha_L */
    alpha_L = sqrt(log(M_SQRT2 * rs_err / p3m.params.accuracy)) / r_cut_iL;
  } else {
    /* even alpha=0 is ok, however, we cannot choose it since it kills the
       k-space error formula.
       Anyways, this very likely NOT the optimal solution */
    alpha_L = 0.1;
  }

  *_alpha_L = alpha_L;
  /* calculate real space and k-space error for this alpha_L */
  rs_err = p3m_real_space_error(coulomb.prefactor, r_cut_iL, p3m.sum_qpart,
                                p3m.sum_q2, alpha_L);
#ifdef CUDA
  if (coulomb.method == COULOMB_P3M_GPU)
    ks_err =
        p3m_k_space_error_gpu(coulomb.prefactor, mesh, cao, p3m.sum_qpart,
                              p3m.sum_q2, alpha_L, box_geo.length().data());
  else
#endif
    ks_err = p3m_k_space_error(coulomb.prefactor, mesh, cao, p3m.sum_qpart,
                               p3m.sum_q2, alpha_L);

  *_rs_err = rs_err;
  *_ks_err = ks_err;
  return sqrt(Utils::sqr(rs_err) + Utils::sqr(ks_err));
}

/** Get the computation time for some @p mesh, @p cao, @p r_cut and @p alpha.
 *
 *  @param[in]  mesh            @copybrief P3MParameters::mesh
 *  @param[in]  cao             @copybrief P3MParameters::cao
 *  @param[in]  r_cut_iL        @copybrief P3MParameters::r_cut_iL
 *  @param[in]  alpha_L         @copybrief P3MParameters::alpha_L
 *
 *  @returns The integration time in case of success, otherwise
 *           -@ref P3M_TUNE_FAIL
 */
static double p3m_mcr_time(const int mesh[3], int cao, double r_cut_iL,
                           double alpha_L) {
  /* rounded up 5000/n_charges timing force evaluations */
  int int_num = (5000 + p3m.sum_qpart) / p3m.sum_qpart;
  double int_time;

  /* broadcast p3m parameters for test run */
  if (coulomb.method != COULOMB_P3M && coulomb.method != COULOMB_ELC_P3M &&
      coulomb.method != COULOMB_P3M_GPU)
    coulomb.method = COULOMB_P3M;

  p3m.params.r_cut = r_cut_iL * box_geo.length()[0];
  p3m.params.r_cut_iL = r_cut_iL;
  p3m.params.mesh[0] = mesh[0];
  p3m.params.mesh[1] = mesh[1];
  p3m.params.mesh[2] = mesh[2];
  p3m.params.cao = cao;
  p3m.params.alpha_L = alpha_L;
  p3m.params.alpha = p3m.params.alpha_L * (1. / box_geo.length()[0]);

  /* initialize p3m structures */
  mpi_bcast_coulomb_params();
  /* perform force calculation test */
  int_time = time_force_calc(int_num);
  if (int_time == -1) {
    return -P3M_TUNE_FAIL;
  }
  return int_time;
}

/** Get the optimal alpha and the corresponding computation time for a fixed
 *  @p mesh and @p cao.
 *
 *  The @p _r_cut_iL is determined via a simple bisection.
 *
 *  @param[out] log             log output
 *  @param[in]  mesh            @copybrief P3MParameters::mesh
 *  @param[in]  cao             @copybrief P3MParameters::cao
 *  @param[in]  r_cut_iL_min    lower bound for @p _r_cut_iL
 *  @param[in]  r_cut_iL_max    upper bound for @p _r_cut_iL
 *  @param[out] _r_cut_iL       @copybrief P3MParameters::r_cut_iL
 *  @param[out] _alpha_L        @copybrief P3MParameters::alpha_L
 *  @param[out] _accuracy       @copybrief P3MParameters::accuracy
 *
 *  @returns The integration time in case of success, otherwise
 *           -@ref P3M_TUNE_FAIL, -@ref P3M_TUNE_ACCURACY_TOO_LARGE,
 *           -@ref P3M_TUNE_CAO_TOO_LARGE, or -@ref P3M_TUNE_ELCTEST
 */
static double p3m_mc_time(char **log, const int mesh[3], int cao,
                          double r_cut_iL_min, double r_cut_iL_max,
                          double *_r_cut_iL, double *_alpha_L,
                          double *_accuracy) {
  double int_time;
  double r_cut_iL;
  double rs_err, ks_err;
  int i, n_cells;
  char b[5 * ES_DOUBLE_SPACE + 3 * ES_INTEGER_SPACE + 128];

  /* initial checks. */
  auto const k_cut =
      std::max(box_geo.length()[0] * cao / (2.0 * mesh[0]),
               std::max(box_geo.length()[1] * cao / (2.0 * mesh[1]),
                        box_geo.length()[2] * cao / (2.0 * mesh[2])));

  auto const min_box_l = *boost::min_element(box_geo.length());
  auto const min_local_box_l = *boost::min_element(local_geo.length());

  if (cao >= std::min(mesh[0], std::min(mesh[1], mesh[2])) ||
      k_cut >= (std::min(min_box_l, min_local_box_l) - skin)) {
    sprintf(b, "%-4d %-3d cao too large for this mesh\n", mesh[0], cao);
    *log = strcat_alloc(*log, b);
    return -P3M_TUNE_CAO_TOO_LARGE;
  }

  /* Either low and high boundary are equal (for fixed cut), or the low border
     is initially 0 and therefore
     has infinite error estimate, as required. Therefore if the high boundary
     fails, there is no possible r_cut */
  if ((*_accuracy = p3m_get_accuracy(mesh, cao, r_cut_iL_max, _alpha_L, &rs_err,
                                     &ks_err)) > p3m.params.accuracy) {
    /* print result */
    sprintf(b, "%-4d %-3d %.5e %.5e %.5e %.3e %.3e accuracy not achieved\n",
            mesh[0], cao, r_cut_iL_max, *_alpha_L, *_accuracy, rs_err, ks_err);
    *log = strcat_alloc(*log, b);
    return -P3M_TUNE_ACCURACY_TOO_LARGE;
  }

  for (;;) {
    r_cut_iL = 0.5 * (r_cut_iL_min + r_cut_iL_max);

    if (r_cut_iL_max - r_cut_iL_min < P3M_RCUT_PREC)
      break;

    /* bisection */
    if ((p3m_get_accuracy(mesh, cao, r_cut_iL, _alpha_L, &rs_err, &ks_err) >
         p3m.params.accuracy))
      r_cut_iL_min = r_cut_iL;
    else
      r_cut_iL_max = r_cut_iL;
  }

  /* final result is always the upper interval boundary, since only there
     we know that the desired minimal accuracy is obtained */
  *_r_cut_iL = r_cut_iL = r_cut_iL_max;

  /* check whether we are running P3M+ELC, and whether we leave a reasonable
   * gap
   * space */
  if (coulomb.method == COULOMB_ELC_P3M &&
      elc_params.gap_size <= 1.1 * r_cut_iL * box_geo.length()[0]) {
    /* print result */
    sprintf(b, "%-4d %-3d %.5e %.5e %.5e %.3e %.3e conflict with ELC\n",
            mesh[0], cao, r_cut_iL, *_alpha_L, *_accuracy, rs_err, ks_err);
    *log = strcat_alloc(*log, b);
    return -P3M_TUNE_ELCTEST;
  }

  int_time = p3m_mcr_time(mesh, cao, r_cut_iL, *_alpha_L);
  if (int_time == -P3M_TUNE_FAIL) {
    *log = strcat_alloc(*log, "tuning failed, test integration not possible\n");
    return int_time;
  }

  *_accuracy =
      p3m_get_accuracy(mesh, cao, r_cut_iL, _alpha_L, &rs_err, &ks_err);

  /* print result */
  sprintf(b, "%-4d %-3d %.5e %.5e %.5e %.3e %.3e %-8.2f\n", mesh[0], cao,
          r_cut_iL, *_alpha_L, *_accuracy, rs_err, ks_err, int_time);
  *log = strcat_alloc(*log, b);
  return int_time;
}

/** Get the optimal alpha and the corresponding computation time for a fixed
 *  @p mesh.
 *
 *  @p _cao should contain an initial guess, which is then adapted by stepping
 *  up and down.
 *
 *  @param[out]     log             log output
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
 *
 *  @returns The integration time in case of success, otherwise
 *           -@ref P3M_TUNE_FAIL or -@ref P3M_TUNE_CAO_TOO_LARGE
 */
static double p3m_m_time(char **log, const int mesh[3], int cao_min,
                         int cao_max, int *_cao, double r_cut_iL_min,
                         double r_cut_iL_max, double *_r_cut_iL,
                         double *_alpha_L, double *_accuracy) {
  double best_time = -1, tmp_time, tmp_r_cut_iL = 0.0, tmp_alpha_L = 0.0,
         tmp_accuracy = 0.0;
  /* in which direction improvement is possible. Initially, we don't know it
   * yet.
   */
  int final_dir = 0;
  int cao = *_cao;

  /* the initial step sets a timing mark. If there is no valid r_cut, we can
     only try
     to increase cao to increase the obtainable precision of the far formula.
     */
  do {
    tmp_time = p3m_mc_time(log, mesh, cao, r_cut_iL_min, r_cut_iL_max,
                           &tmp_r_cut_iL, &tmp_alpha_L, &tmp_accuracy);
    /* bail out if the force evaluation is not working */
    if (tmp_time == -P3M_TUNE_FAIL)
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

  /* at the boundaries, only the opposite direction can be used for
   * optimisation
   */
  if (cao == cao_min)
    final_dir = 1;
  else if (cao == cao_max)
    final_dir = -1;

  if (final_dir == 0) {
    /* check in which direction we can optimise. Both directions are possible
     */
    double dir_times[3];
    for (final_dir = -1; final_dir <= 1; final_dir += 2) {
      dir_times[final_dir + 1] = tmp_time =
          p3m_mc_time(log, mesh, cao + final_dir, r_cut_iL_min, r_cut_iL_max,
                      &tmp_r_cut_iL, &tmp_alpha_L, &tmp_accuracy);
      /* bail out on errors, as usual */
      if (tmp_time == -P3M_TUNE_FAIL)
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
       * worse, we can still try */
      /* down is possible and not much worse, while up is either illegal or
       * even
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
    tmp_time = p3m_mc_time(log, mesh, cao, r_cut_iL_min, r_cut_iL_max,
                           &tmp_r_cut_iL, &tmp_alpha_L, &tmp_accuracy);
    /* bail out on errors, as usual */
    if (tmp_time == -P3M_TUNE_FAIL)
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

int p3m_adaptive_tune(char **log) {
  int mesh[3] = {0, 0, 0};
  int tmp_mesh[3];
  double r_cut_iL_min, r_cut_iL_max, r_cut_iL = -1, tmp_r_cut_iL = 0.0;
  int cao_min, cao_max, cao = -1, tmp_cao;
  double alpha_L = -1, tmp_alpha_L = 0.0;
  double accuracy = -1, tmp_accuracy = 0.0;
  double time_best = 1e20, tmp_time;
  double mesh_density = 0.0, mesh_density_min, mesh_density_max;
  char b[3 * ES_INTEGER_SPACE + 3 * ES_DOUBLE_SPACE + 128];
  bool tune_mesh = false; // indicates if mesh should be tuned

  if (p3m.params.epsilon != P3M_EPSILON_METALLIC) {
    if (!((box_geo.length()[0] == box_geo.length()[1]) &&
          (box_geo.length()[1] == box_geo.length()[2]))) {
      *log = strcat_alloc(
          *log, "{049 P3M_init: Nonmetallic epsilon requires cubic box} ");
      return ES_ERROR;
    }
  }

  if (p3m_sanity_checks_system(node_grid)) {
    return ES_ERROR;
  }

  /* preparation */
  mpi_call(p3m_count_charged_particles);
  p3m_count_charged_particles();

  /* Print Status */
  sprintf(b, "P3M tune parameters: Accuracy goal = %.5e prefactor = %.5e\n",
          p3m.params.accuracy, coulomb.prefactor);
  *log = strcat_alloc(*log, b);
  sprintf(b, "System: box_l = %.5e # charged part = %d Sum[q_i^2] = %.5e\n",
          box_geo.length()[0], p3m.sum_qpart, p3m.sum_q2);
  *log = strcat_alloc(*log, b);

  if (p3m.sum_qpart == 0) {
    *log = strcat_alloc(*log,
                        "no charged particles in the system, cannot tune P3M");
    return ES_ERROR;
  }

  /* Activate tuning mode */
  p3m.params.tuning = true;

  /* parameter ranges */
  /* if at least the number of meshpoints in one direction is not set, we have
   * to tune it. */
  if (p3m.params.mesh[0] == 0 || p3m.params.mesh[1] == 0 ||
      p3m.params.mesh[2] == 0) {
    /* Medium-educated guess for the minimal mesh */
    mesh_density_min =
        pow(p3m.sum_qpart / (box_geo.length()[0] * box_geo.length()[1] *
                             box_geo.length()[2]),
            1.0 / 3.0);
    mesh_density_max = 512 / pow(box_geo.length()[0] * box_geo.length()[1] *
                                     box_geo.length()[2],
                                 1.0 / 3.0);
    tune_mesh = true;
    /* this limits the tried meshes if the accuracy cannot
       be obtained with smaller meshes, but normally not all these
       meshes have to be tested */
    /* avoid using more than 1 GB of FFT arrays (per default, see config.hpp)
     */
  } else if (p3m.params.mesh[1] == -1 && p3m.params.mesh[2] == -1) {
    mesh_density = mesh_density_min = mesh_density_max =
        p3m.params.mesh[0] / box_geo.length()[0];
    p3m.params.mesh[1] =
        static_cast<int>(std::round(mesh_density * box_geo.length()[1]));
    p3m.params.mesh[2] =
        static_cast<int>(std::round(mesh_density * box_geo.length()[2]));
    if (p3m.params.mesh[1] % 2 == 1)
      p3m.params.mesh[1]++; // Make sure that the mesh is even in all directions
    if (p3m.params.mesh[2] % 2 == 1)
      p3m.params.mesh[2]++;

    sprintf(b, "fixed mesh %d %d %d\n", p3m.params.mesh[0], p3m.params.mesh[1],
            p3m.params.mesh[2]);
    *log = strcat_alloc(*log, b);
  } else {
    mesh_density_min = mesh_density_max =
        p3m.params.mesh[0] / box_geo.length()[0];

    sprintf(b, "fixed mesh %d %d %d\n", p3m.params.mesh[0], p3m.params.mesh[1],
            p3m.params.mesh[2]);
    *log = strcat_alloc(*log, b);
  }

  if (p3m.params.r_cut_iL == 0.0) {
    auto const min_box_l = *boost::min_element(box_geo.length());
    auto const min_local_box_l = *boost::min_element(local_geo.length());

    r_cut_iL_min = 0;
    r_cut_iL_max = std::min(min_local_box_l, min_box_l / 2.0) - skin;
    r_cut_iL_min *= (1. / box_geo.length()[0]);
    r_cut_iL_max *= (1. / box_geo.length()[0]);
  } else {
    r_cut_iL_min = r_cut_iL_max = p3m.params.r_cut_iL;

    sprintf(b, "fixed r_cut_iL %f\n", p3m.params.r_cut_iL);
    *log = strcat_alloc(*log, b);
  }

  if (p3m.params.cao == 0) {
    cao_min = 1;
    cao_max = 7;
    cao = cao_max;
  } else {
    cao_min = cao_max = cao = p3m.params.cao;

    sprintf(b, "fixed cao %d\n", p3m.params.cao);
    *log = strcat_alloc(*log, b);
  }

  *log = strcat_alloc(*log, "mesh cao r_cut_iL     alpha_L      err          "
                            "rs_err     ks_err     time [ms]\n");

  /* mesh loop */
  /* we're tuning the density of mesh points, which is the same in every
   * direction. */
  for (mesh_density = mesh_density_min; mesh_density <= mesh_density_max;
       mesh_density += 0.1) {
    tmp_cao = cao;

    if (tune_mesh) {
      tmp_mesh[0] =
          static_cast<int>(std::round(box_geo.length()[0] * mesh_density));
      tmp_mesh[1] =
          static_cast<int>(std::round(box_geo.length()[1] * mesh_density));
      tmp_mesh[2] =
          static_cast<int>(std::round(box_geo.length()[2] * mesh_density));
    } else {
      tmp_mesh[0] = p3m.params.mesh[0];
      tmp_mesh[1] = p3m.params.mesh[1];
      tmp_mesh[2] = p3m.params.mesh[2];
    }

    if (tmp_mesh[0] % 2) // Make sure that the mesh is even in all directions
      tmp_mesh[0]++;
    if (tmp_mesh[1] % 2)
      tmp_mesh[1]++;
    if (tmp_mesh[2] % 2)
      tmp_mesh[2]++;

#ifdef HIP
    // When running on HIP, we don't support mesh sizes whose prime factors are
    // not 2, 3 or 5. So we skip the other supported prime factors during
    // tuning.
    if (tune_mesh && (tmp_mesh[0] % 7 == 0 || tmp_mesh[0] % 11 == 0 ||
                      tmp_mesh[0] % 13 == 0 || tmp_mesh[1] % 7 == 0 ||
                      tmp_mesh[1] % 11 == 0 || tmp_mesh[1] % 13 == 0 ||
                      tmp_mesh[2] % 7 == 0 || tmp_mesh[2] % 11 == 0 ||
                      tmp_mesh[2] % 13 == 0)) {
      continue;
    }
#endif

    tmp_time =
        p3m_m_time(log, tmp_mesh, cao_min, cao_max, &tmp_cao, r_cut_iL_min,
                   r_cut_iL_max, &tmp_r_cut_iL, &tmp_alpha_L, &tmp_accuracy);
    /* some error occurred during the tuning force evaluation */
    /* this mesh does not work at all */
    if (tmp_time < 0.0)
      continue;

    /* the optimum r_cut for this mesh is the upper limit for higher meshes,
       everything else is slower */
    if (coulomb.method == COULOMB_P3M)
      r_cut_iL_max = tmp_r_cut_iL;

    /* new optimum */
    if (tmp_time < time_best) {
      time_best = tmp_time;
      mesh[0] = tmp_mesh[0];
      mesh[1] = tmp_mesh[1];
      mesh[2] = tmp_mesh[2];
      cao = tmp_cao;
      r_cut_iL = tmp_r_cut_iL;
      alpha_L = tmp_alpha_L;
      accuracy = tmp_accuracy;
    }
    /* no hope of further optimisation */
    else if (tmp_time > time_best + P3M_TIME_GRAN) {
      break;
    }
  }

  if (time_best == 1e20) {
    *log = strcat_alloc(*log,
                        "failed to tune P3M parameters to required accuracy\n");
    return ES_ERROR;
  }

  /* set tuned p3m parameters */
  p3m.params.tuning = false;
  p3m.params.r_cut = r_cut_iL * box_geo.length()[0];
  p3m.params.r_cut_iL = r_cut_iL;
  p3m.params.mesh[0] = mesh[0];
  p3m.params.mesh[1] = mesh[1];
  p3m.params.mesh[2] = mesh[2];
  p3m.params.cao = cao;
  p3m.params.alpha_L = alpha_L;
  p3m.params.alpha = p3m.params.alpha_L * (1. / box_geo.length()[0]);
  p3m.params.accuracy = accuracy;
  /* broadcast tuned p3m parameters */
  mpi_bcast_coulomb_params();

  /* Tell the user about the outcome */
  sprintf(b,
          "\nresulting parameters: mesh: (%d %d %d), cao: %d, r_cut_iL: %.4e,"
          "\n                      alpha_L: %.4e, accuracy: %.4e, time: %.2f\n",
          mesh[0], mesh[1], mesh[2], cao, r_cut_iL, alpha_L, accuracy,
          time_best);
  *log = strcat_alloc(*log, b);
  return ES_OK;
}

void p3m_count_charged_particles() {
  double node_sums[3], tot_sums[3];

  for (int i = 0; i < 3; i++) {
    node_sums[i] = 0.0;
    tot_sums[i] = 0.0;
  }

  for (auto const &p : cell_structure.local_particles()) {
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

REGISTER_CALLBACK(p3m_count_charged_particles)

double p3m_real_space_error(double prefac, double r_cut_iL, int n_c_part,
                            double sum_q2, double alpha_L) {
  return (2.0 * prefac * sum_q2 * exp(-Utils::sqr(r_cut_iL * alpha_L))) /
         sqrt((double)n_c_part * r_cut_iL * box_geo.length()[0] *
              box_geo.length()[0] * box_geo.length()[1] * box_geo.length()[2]);
}

double p3m_k_space_error(double prefac, const int mesh[3], int cao,
                         int n_c_part, double sum_q2, double alpha_L) {
  int nx, ny, nz;
  double he_q = 0.0, mesh_i[3] = {1.0 / mesh[0], 1.0 / mesh[1], 1.0 / mesh[2]},
         alpha_L_i = 1. / alpha_L;
  double alias1, alias2, n2, cs;
  double ctan_x, ctan_y;

  for (nx = -mesh[0] / 2; nx < mesh[0] / 2; nx++) {
    ctan_x = p3m_analytic_cotangent_sum(nx, mesh_i[0], cao);
    for (ny = -mesh[1] / 2; ny < mesh[1] / 2; ny++) {
      ctan_y = ctan_x * p3m_analytic_cotangent_sum(ny, mesh_i[1], cao);
      for (nz = -mesh[2] / 2; nz < mesh[2] / 2; nz++) {
        if ((nx != 0) || (ny != 0) || (nz != 0)) {
          n2 = Utils::sqr(nx) + Utils::sqr(ny) + Utils::sqr(nz);
          cs = p3m_analytic_cotangent_sum(nz, mesh_i[2], cao) * ctan_y;
          p3m_tune_aliasing_sums(nx, ny, nz, mesh, mesh_i, cao, alpha_L_i,
                                 &alias1, &alias2);

          double d = alias1 - Utils::sqr(alias2 / cs) / n2;
          /* at high precisions, d can become negative due to extinction;
             also, don't take values that have no significant digits left*/
          if (d > 0 && (fabs(d / alias1) > ROUND_ERROR_PREC))
            he_q += d;
        }
      }
    }
  }
  return 2.0 * prefac * sum_q2 * sqrt(he_q / (double)n_c_part) /
         (box_geo.length()[1] * box_geo.length()[2]);
}

void p3m_tune_aliasing_sums(int nx, int ny, int nz, const int mesh[3],
                            const double mesh_i[3], int cao, double alpha_L_i,
                            double *alias1, double *alias2) {

  int mx, my, mz;
  double nmx, nmy, nmz;
  double fnmx, fnmy, fnmz;

  double ex, ex2, nm2, U2, factor1;

  factor1 = Utils::sqr(Utils::pi() * alpha_L_i);

  *alias1 = *alias2 = 0.0;
  for (mx = -P3M_BRILLOUIN; mx <= P3M_BRILLOUIN; mx++) {
    fnmx = mesh_i[0] * (nmx = nx + mx * mesh[0]);
    for (my = -P3M_BRILLOUIN; my <= P3M_BRILLOUIN; my++) {
      fnmy = mesh_i[1] * (nmy = ny + my * mesh[1]);
      for (mz = -P3M_BRILLOUIN; mz <= P3M_BRILLOUIN; mz++) {
        fnmz = mesh_i[2] * (nmz = nz + mz * mesh[2]);

        nm2 = Utils::sqr(nmx) + Utils::sqr(nmy) + Utils::sqr(nmz);
        ex2 = Utils::sqr(ex = exp(-factor1 * nm2));

        U2 = pow(sinc(fnmx) * sinc(fnmy) * sinc(fnmz), 2.0 * cao);

        *alias1 += ex2 / nm2;
        *alias2 += U2 * ex * (nx * nmx + ny * nmy + nz * nmz) / nm2;
      }
    }
  }
}

void p3m_init_a_ai_cao_cut() {
  int i;
  for (i = 0; i < 3; i++) {
    p3m.params.ai[i] = (double)p3m.params.mesh[i] / box_geo.length()[i];
    p3m.params.a[i] = 1.0 / p3m.params.ai[i];
    p3m.params.cao_cut[i] = 0.5 * p3m.params.a[i] * p3m.params.cao;
  }
}

bool p3m_sanity_checks_boxl() {
  int i;
  bool ret = false;
  for (i = 0; i < 3; i++) {
    /* check k-space cutoff */
    if (p3m.params.cao_cut[i] >= 0.5 * box_geo.length()[i]) {
      runtimeErrorMsg() << "P3M_init: k-space cutoff " << p3m.params.cao_cut[i]
                        << " is larger than half of box dimension "
                        << box_geo.length()[i];
      ret = true;
    }
    if (p3m.params.cao_cut[i] >= local_geo.length()[i]) {
      runtimeErrorMsg() << "P3M_init: k-space cutoff " << p3m.params.cao_cut[i]
                        << " is larger than local box dimension "
                        << local_geo.length()[i];
      ret = true;
    }
  }

  return ret;
}

/**
 * @brief General sanity checks independent of p3m parameters.
 *
 * @return false if ok, true on error.
 */
bool p3m_sanity_checks_system(const Utils::Vector3i &grid) {
  bool ret = false;

  if (!box_geo.periodic(0) || !box_geo.periodic(1) || !box_geo.periodic(2)) {
    runtimeErrorMsg() << "P3M requires periodicity 1 1 1";
    ret = true;
  }

  if (cell_structure.type != CELL_STRUCTURE_DOMDEC) {
    runtimeErrorMsg()
        << "P3M at present requires the domain decomposition cell system";
    ret = true;
  }

  if (grid[0] < grid[1] || grid[1] < grid[2]) {
    runtimeErrorMsg() << "P3M_init: node grid must be sorted, largest first";
    ret = true;
  }

  if (p3m.params.epsilon != P3M_EPSILON_METALLIC) {
    if (!((p3m.params.mesh[0] == p3m.params.mesh[1]) &&
          (p3m.params.mesh[1] == p3m.params.mesh[2]))) {
      runtimeErrorMsg() << "P3M_init: Nonmetallic epsilon requires cubic box";
      ret = true;
    }
  }

  return ret;
}

bool p3m_sanity_checks() {
  bool ret = false;

  if (p3m_sanity_checks_system(node_grid))
    ret = true;

  if (p3m_sanity_checks_boxl())
    ret = true;

  if (p3m.params.mesh[0] == 0) {
    runtimeErrorMsg() << "P3M_init: mesh size is not yet set";
    ret = true;
  }
  if (p3m.params.cao == 0) {
    runtimeErrorMsg() << "P3M_init: cao is not yet set";
    ret = true;
  }
  if (p3m.params.alpha < 0.0) {
    runtimeErrorMsg() << "P3M_init: alpha must be >0";
    ret = true;
  }

  return ret;
}

void p3m_scaleby_box_l() {
  if (coulomb.prefactor < 0.0) {
    runtimeErrorMsg() << "The Coulomb prefactor has to be >=0";
    return;
  }

  p3m.params.r_cut = p3m.params.r_cut_iL * box_geo.length()[0];
  p3m.params.alpha = p3m.params.alpha_L * (1. / box_geo.length()[0]);
  p3m_init_a_ai_cao_cut();
  p3m_calc_lm_ld_pos(p3m.local_mesh, p3m.params);
  p3m_sanity_checks_boxl();
  p3m_calc_influence_function_force();
  p3m_calc_influence_function_energy();
}

#endif /* of P3M */
