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
#ifndef _P3M_MAGNETOSTATICS_H
#define _P3M_MAGNETOSTATICS_H
/** \file
 *  P3M algorithm for long range magnetic dipole-dipole interaction.
 *
 *  We use here a P3M (Particle-Particle Particle-Mesh) method based
 *  on the dipolar Ewald summation. Details of the used method can be found in
 *  @cite hockney88a and @cite deserno98a @cite deserno98b. The file p3m
 *  contains only the Particle-Mesh part.
 *
 *  Further reading: @cite cerda08d
 *
 *  Implementation in p3m-dipolar.cpp.
 */

#include "config.hpp"

#ifdef DP3M
#include "Particle.hpp"
#include "electrostatics_magnetostatics/dipole.hpp"
#include "fft.hpp"
#include "p3m-common.hpp"
#include "p3m_interpolation.hpp"
#include "p3m_send_mesh.hpp"

#include <ParticleRange.hpp>
#include <utils/constants.hpp>
#include <utils/math/AS_erfc_part.hpp>

struct dp3m_data_struct {
  dp3m_data_struct();

  P3MParameters params;

  /** local mesh. */
  p3m_local_mesh local_mesh;
  /** real space mesh (local) for CA/FFT.*/
  fft_vector<double> rs_mesh;
  /** real space mesh (local) for CA/FFT of the dipolar field.*/
  std::array<fft_vector<double>, 3> rs_mesh_dip;
  /** k-space mesh (local) for k-space calculation and FFT.*/
  std::vector<double> ks_mesh;

  /** number of dipolar particles (only on master node). */
  int sum_dip_part;
  /** Sum of square of magnetic dipoles (only on master node). */
  double sum_mu2;

  /** position shift for calc. of first assignment mesh point. */
  double pos_shift;
  /** help variable for calculation of aliasing sums */
  std::vector<double> meshift;

  /** Spatial differential operator in k-space. We use an i*k differentiation.
   */
  std::vector<double> d_op;
  /** Force optimised influence function (k-space) */
  std::vector<double> g_force;
  /** Energy optimised influence function (k-space) */
  std::vector<double> g_energy;

  p3m_interpolation_cache inter_weights;

  /** number of permutations in k_space */
  int ks_pnum;

  /** send/recv mesh sizes */
  p3m_send_mesh sm;

  /* Stores the value of the energy correction due to MS effects */
  double energy_correction;

  fft_data_struct fft;
};

/** dipolar P3M parameters. */
extern dp3m_data_struct dp3m;

/** \name Exported Functions */
/************************************************************/
/*@{*/

/** @copydoc p3m_set_tune_params */
void dp3m_set_tune_params(double r_cut, int mesh, int cao, double alpha,
                          double accuracy, int n_interpol);

/** @copydoc p3m_set_params */
int dp3m_set_params(double r_cut, int mesh, int cao, double alpha,
                    double accuracy);

/** @copydoc p3m_set_mesh_offset */
int dp3m_set_mesh_offset(double x, double y, double z);

/** @copydoc p3m_set_eps */
int dp3m_set_eps(double eps);

/** Initialize all structures, parameters and arrays needed for the
 *  P3M algorithm for dipole-dipole interactions.
 */
void dp3m_init();

/** @copydoc p3m_scaleby_box_l */
void dp3m_scaleby_box_l();

/** Sanity checks */
bool dp3m_sanity_checks(const Utils::Vector3i &grid);

/** Assign the physical dipoles using the tabulated assignment function.
 *  If Dstore_ca_frac is true, then the charge fractions are buffered in
 *  Dcur_ca_fmp and Dcur_ca_frac.
 */
void dp3m_dipole_assign(const ParticleRange &particles);

/** Reset @ref dp3m core parameters */
void dp3m_deactivate();

/** Tune dipolar P3M parameters to desired accuracy.
 *
 *  The parameters
 *  @ref P3MParameters::mesh "mesh",
 *  @ref P3MParameters::cao "cao",
 *  @ref P3MParameters::r_cut_iL "r_cut_iL" and
 *  @ref P3MParameters::alpha_L "alpha_L"
 *  are tuned to obtain the target accuracy (initially stored in
 *  @ref P3MParameters::accuracy "accuracy") in optimal time.
 *  These parameters are stored in the @ref dp3m object.
 *
 *  The function utilizes the analytic expression of the error estimate
 *  for the dipolar P3M method in the paper of @cite cerda08d in
 *  order to obtain the rms error in the force for a system of N randomly
 *  distributed particles in a cubic box. For the real space error, the
 *  estimate in @cite kolafa92a is used.
 *
 *  Parameter ranges if not given explicit values via dp3m_set_tune_params():
 *  - @p mesh is set up such that the number of mesh points is equal to the
 *    number of magnetic dipolar particles
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
 *  The function is based on routines of the program HE_Q.cpp for charges
 *  written by M. Deserno.
 *
 *  @param[out]  log  log output
 *  @retval ES_OK
 *  @retval ES_ERROR
 */
int dp3m_adaptive_tune(char **log);

/** Compute the k-space part of forces and energies for the magnetic
 *  dipole-dipole interaction
 */
double dp3m_calc_kspace_forces(bool force_flag, bool energy_flag,
                               ParticleRange const &particles);

/** Calculate number of magnetic particles, the sum of the squared
 *  charges and the squared sum of the charges.
 */
void dp3m_count_magnetic_particles();

/** Calculate real space contribution of p3m dipolar pair forces and torques.
 *  If NPT is compiled in, it returns the energy, which is needed for NPT.
 */
inline std::tuple<double, Utils::Vector3d, Utils::Vector3d, Utils::Vector3d>
dp3m_pair_force(Particle const &p1, Particle const &p2,
                Utils::Vector3d const &d, double dist2, double dist) {
  if ((p1.p.dipm == 0.) || (p2.p.dipm == 0.))
    return std::make_tuple(0.0, Utils::Vector3d{}, Utils::Vector3d{},
                           Utils::Vector3d{});

  if (dist < dp3m.params.r_cut && dist > 0) {
    auto const dip1 = p1.calc_dip();
    auto const dip2 = p2.calc_dip();
    auto const alpsq = dp3m.params.alpha * dp3m.params.alpha;
    auto const adist = dp3m.params.alpha * dist;
#if USE_ERFC_APPROXIMATION
    auto const erfc_part_ri = Utils::AS_erfc_part(adist) / dist;
#else
    auto const erfc_part_ri = erfc(adist) / dist;
#endif

    // Calculate scalar multiplications for vectors mi, mj, rij
    auto const mimj = dip1 * dip2;

    auto const mir = dip1 * d;
    auto const mjr = dip2 * d;

    auto const coeff = 2.0 * dp3m.params.alpha * Utils::sqrt_pi_i();
    auto const dist2i = 1 / dist2;
    auto const exp_adist2 = exp(-adist * adist);

    double B_r, C_r, D_r;
    if (dp3m.params.accuracy > 5e-06)
      B_r = (erfc_part_ri + coeff) * exp_adist2 * dist2i;
    else
      B_r = (erfc(adist) / dist + coeff * exp_adist2) * dist2i;

    C_r = (3 * B_r + 2 * alpsq * coeff * exp_adist2) * dist2i;
    D_r = (5 * C_r + 4 * coeff * alpsq * alpsq * exp_adist2) * dist2i;

    // Calculate real-space forces
    auto const force =
        dipole.prefactor *
        ((mimj * d + dip1 * mjr + dip2 * mir) * C_r - mir * mjr * D_r * d);

#ifdef ROTATION
    // Calculate vector multiplications for vectors mi, mj, rij
    auto const mixmj = vector_product(dip1, dip2);
    auto const mixr = vector_product(dip1, d);
    auto const mjxr = vector_product(dip2, d);

    // Calculate real-space torques
    auto const torque1 = dipole.prefactor * (-mixmj * B_r + mixr * (mjr * C_r));
    auto const torque2 = dipole.prefactor * (mixmj * B_r + mjxr * (mir * C_r));
#endif
#ifdef NPT
#if USE_ERFC_APPROXIMATION
    auto const fac = dipole.prefactor * p1.p.dipm * p2.p.dipm * exp_adist2;
#else
    auto const fac = dipole.prefactor * p1.p.dipm * p2.p.dipm;
#endif
    auto const _energy = fac * (mimj * B_r - mir * mjr * C_r);
#else
    auto const _energy = 0.0;
#endif
    return std::make_tuple(_energy, force, torque1, torque2);
  }
  return std::make_tuple(0.0, Utils::Vector3d{}, Utils::Vector3d{},
                         Utils::Vector3d{});
}

/** Calculate real space contribution of dipolar pair energy. */
inline double dp3m_pair_energy(Particle const &p1, Particle const &p2,
                               Utils::Vector3d const &d, double dist2,
                               double dist) {
  auto const dip1 = p1.calc_dip();
  auto const dip2 = p2.calc_dip();

  if (dist < dp3m.params.r_cut && dist > 0) {
    auto const alpsq = dp3m.params.alpha * dp3m.params.alpha;
    auto const adist = dp3m.params.alpha * dist;
    /*fac1 = dipole.prefactor;*/

#if USE_ERFC_APPROXIMATION
    auto const erfc_part_ri = Utils::AS_erfc_part(adist) / dist;
    /*  fac1 = dipole.prefactor * p1.p.dipm*p2.p.dipm; IT WAS WRONG */
    /* *exp(-adist*adist); */
#else
    auto const erfc_part_ri = erfc(adist) / dist;
    /* fac1 = dipole.prefactor * p1.p.dipm*p2.p.dipm;  IT WAS WRONG*/
#endif

    // Calculate scalar multiplications for vectors mi, mj, rij
    auto const mimj = dip1 * dip2;
    auto const mir = dip1 * d;
    auto const mjr = dip2 * d;

    auto const coeff = 2.0 * dp3m.params.alpha * Utils::sqrt_pi_i();
    auto const dist2i = 1 / dist2;
    auto const exp_adist2 = exp(-adist * adist);

    double B_r;
    if (dp3m.params.accuracy > 5e-06)
      B_r = (erfc_part_ri + coeff) * exp_adist2 * dist2i;
    else
      B_r = (erfc(adist) / dist + coeff * exp_adist2) * dist2i;

    auto const C_r = (3 * B_r + 2 * alpsq * coeff * exp_adist2) * dist2i;

    /*
      printf("(%4i %4i) pair energy = %f (B_r=%15.12f C_r=%15.12f)\n",
      p1.p.identity,p2.p.identity,fac1*(mimj*B_r-mir*mjr*C_r),B_r,C_r);
    */

    /* old line return fac1 * ( mimj*B_r - mir*mjr * C_r );*/
    return dipole.prefactor * (mimj * B_r - mir * mjr * C_r);
  }
  return 0.0;
}

/*@}*/

#endif /* DP3M */
#endif /* _P3M_DIPOLES_H */
