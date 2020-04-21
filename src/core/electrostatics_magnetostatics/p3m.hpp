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
#ifndef _P3M_H
#define _P3M_H
/** \file
 *  P3M algorithm for long range Coulomb interaction.
 *
 *  We use a P3M (Particle-Particle Particle-Mesh) method based on the
 *  Ewald summation. Details of the used method can be found in
 *  @cite hockney88a and @cite deserno98a @cite deserno98b.
 *
 *  Further reading: @cite ewald21a, @cite hockney88a, @cite deserno98a,
 *  @cite deserno98b, @cite deserno00e, @cite deserno00b, @cite cerda08d.
 *
 *  Implementation in p3m.cpp.
 */

#include "config.hpp"

#ifdef P3M

#include "fft.hpp"
#include "p3m-common.hpp"
#include "p3m_interpolation.hpp"
#include "p3m_send_mesh.hpp"

#include <ParticleRange.hpp>
#include <utils/constants.hpp>
#include <utils/math/AS_erfc_part.hpp>

/************************************************
 * data types
 ************************************************/

struct p3m_data_struct {
  p3m_data_struct();

  P3MParameters params;

  /** local mesh. */
  p3m_local_mesh local_mesh;
  /** real space mesh (local) for CA/FFT.*/
  fft_vector<double> rs_mesh;
  /** mesh (local) for the electric field.*/
  std::array<fft_vector<double>, 3> E_mesh;

  /** number of charged particles (only on master node). */
  int sum_qpart;
  /** Sum of square of charges (only on master node). */
  double sum_q2;
  /** square of sum of charges (only on master node). */
  double square_sum_q;

  /** help variable for calculation of aliasing sums */
  std::vector<double> meshift_x;
  std::vector<double> meshift_y;
  std::vector<double> meshift_z;

  /** Spatial differential operator in k-space. We use an i*k differentiation.
   */
  std::array<std::vector<double>, 3> d_op;
  /** Force optimised influence function (k-space) */
  std::vector<double> g_force;
  /** Energy optimised influence function (k-space) */
  std::vector<double> g_energy;

  p3m_interpolation_cache inter_weights;

  /** number of permutations in k_space */
  int ks_pnum;

  /** send/recv mesh sizes */
  p3m_send_mesh sm;

  fft_data_struct fft;
};

/** P3M parameters. */
extern p3m_data_struct p3m;

/** Tune P3M parameters to desired accuracy.
 *
 *  The parameters
 *  @ref P3MParameters::mesh "mesh",
 *  @ref P3MParameters::cao "cao",
 *  @ref P3MParameters::r_cut_iL "r_cut_iL" and
 *  @ref P3MParameters::alpha_L "alpha_L"
 *  are tuned to obtain the target accuracy (initially stored in
 *  @ref P3MParameters::accuracy "accuracy") in optimal time.
 *  These parameters are stored in the @ref p3m object.
 *
 *  The function utilizes the analytic expression of the error estimate
 *  for the P3M method in @cite hockney88a (eq. (8.23)) in
 *  order to obtain the rms error in the force for a system of N randomly
 *  distributed particles in a cubic box.
 *  For the real space error the estimate of Kolafa/Perram is used.
 *
 *  Parameter ranges if not given explicit values via p3m_set_tune_params():
 *  - @p mesh is set up such that the number of mesh points is equal to the
 *    number of charged particles
 *  - @p cao explores all possible values
 *  - @p alpha_L is tuned for each tuple (@p r_cut_iL, @p mesh, @p cao) and
 *    calculated assuming that the error contributions of real and reciprocal
 *    space should be equal
 *
 *  After checking if the total error lies below the target accuracy, the
 *  time needed for one force calculation (including Verlet list update)
 *  is measured via time_force_calc().
 *
 *  The function generates a log of the performed tuning.
 *
 *  The function is based on routines of the program HE_Q.cpp written by M.
 *  Deserno.
 *
 *  @param[out]  log  log output
 *  @retval ES_OK
 *  @retval ES_ERROR
 */
int p3m_adaptive_tune(char **log);

/** Initialize all structures, parameters and arrays needed for the
 *  P3M algorithm for charge-charge interactions.
 */
void p3m_init();

/** Update @ref P3MParameters::alpha "alpha" and
 *  @ref P3MParameters::r_cut "r_cut" if box length changed
 */
void p3m_scaleby_box_l();

/** Compute the k-space part of forces and energies for the charge-charge
 *  interaction
 */
double p3m_calc_kspace_forces(bool force_flag, bool energy_flag,
                              const ParticleRange &particles);

/** Compute the k-space part of the stress tensor **/
Utils::Vector9d p3m_calc_kspace_stress();

/** Sanity checks */
bool p3m_sanity_checks();

/** Calculate number of charged particles, the sum of the squared
 *  charges and the squared sum of the charges.
 */
void p3m_count_charged_particles();

/** Assign the physical charges using the tabulated charge assignment function.
 */
void p3m_charge_assign(const ParticleRange &particles);

/** Assign a single charge into the current charge grid.
 *
 *  @param[in] q          %Particle charge
 *  @param[in] real_pos   %Particle position in real space
 *  @param[in] inter_weights Cached interpolation weights to be used.
 */
void p3m_assign_charge(double q, const Utils::Vector3d &real_pos,
                       p3m_interpolation_cache &inter_weights);
/** @overload */
void p3m_assign_charge(double q, const Utils::Vector3d &real_pos);

/** Calculate real space contribution of Coulomb pair forces. */
inline void p3m_add_pair_force(double q1q2, Utils::Vector3d const &d,
                               double dist, Utils::Vector3d &force) {
  if (dist < p3m.params.r_cut) {
    if (dist > 0.0) {
      double adist = p3m.params.alpha * dist;
#if USE_ERFC_APPROXIMATION
      auto const erfc_part_ri = Utils::AS_erfc_part(adist) / dist;
      auto const fac1 = q1q2 * exp(-adist * adist);
      auto const fac2 =
          fac1 * (erfc_part_ri + 2.0 * p3m.params.alpha * Utils::sqrt_pi_i()) /
          (dist * dist);
#else
      auto const erfc_part_ri = erfc(adist) / dist;
      auto const fac1 = q1q2;
      auto const fac2 =
          fac1 *
          (erfc_part_ri +
           2.0 * p3m.params.alpha * Utils::sqrt_pi_i() * exp(-adist * adist)) /
          (dist * dist);
#endif
      force += fac2 * d;
    }
  }
}

/** Set initial values for p3m_adaptive_tune()
 *
 *  @param[in]  r_cut        @copybrief P3MParameters::r_cut
 *  @param[in]  mesh         @copybrief P3MParameters::mesh
 *  @param[in]  cao          @copybrief P3MParameters::cao
 *  @param[in]  alpha        @copybrief P3MParameters::alpha
 *  @param[in]  accuracy     @copybrief P3MParameters::accuracy
 */
void p3m_set_tune_params(double r_cut, const int mesh[3], int cao, double alpha,
                         double accuracy);

/** Set custom parameters
 *
 *  @param[in]  r_cut        @copybrief P3MParameters::r_cut
 *  @param[in]  mesh         @copybrief P3MParameters::mesh
 *  @param[in]  cao          @copybrief P3MParameters::cao
 *  @param[in]  alpha        @copybrief P3MParameters::alpha
 *  @param[in]  accuracy     @copybrief P3MParameters::accuracy
 *  @return Custom error code
 */
int p3m_set_params(double r_cut, const int *mesh, int cao, double alpha,
                   double accuracy);

/** Set mesh offset
 *
 *  @param[in]  x , y , z  Components of @ref P3MParameters::mesh_off
 *                         "mesh_off"
 */
int p3m_set_mesh_offset(double x, double y, double z);

/** Set @ref P3MParameters::epsilon "epsilon" parameter
 *
 *  @param[in]  eps          @copybrief P3MParameters::epsilon
 */
int p3m_set_eps(double eps);

/** Calculate real space contribution of Coulomb pair energy. */
inline double p3m_pair_energy(double chgfac, double dist) {
  if (dist < p3m.params.r_cut && dist != 0) {
    double adist = p3m.params.alpha * dist;
#if USE_ERFC_APPROXIMATION
    double erfc_part_ri = Utils::AS_erfc_part(adist) / dist;
    return chgfac * erfc_part_ri * exp(-adist * adist);
#else
    double erfc_part_ri = erfc(adist) / dist;
    return chgfac * erfc_part_ri;
#endif
  }
  return 0.0;
}

#endif /* of ifdef P3M */

#endif /*of ifndef P3M_H */
