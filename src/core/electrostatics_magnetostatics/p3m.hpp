/*
  Copyright (C) 2010-2018 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
    Max-Planck-Institute for Polymer Research, Theory Group

  This file is part of ESPResSo.

  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef _P3M_H
#define _P3M_H
/** \file
 *  P3M algorithm for long range Coulomb interaction.
 *
 *  We use a P3M (Particle-Particle Particle-Mesh) method based on the
 *  Ewald summation. Details of the used method can be found in
 *  Hockney/Eastwood and Deserno/Holm.
 *
 *  Further reading:
 *  - P. P. Ewald,
 *    *Die Berechnung optischer und elektrostatischer Gitterpotentiale*,
 *    Ann. Phys. (64) 253-287, 1921
 *  - R. W. Hockney and J. W. Eastwood,
 *    *Computer simulation using particles*,
 *    IOP, London, 1988
 *  - M. Deserno and C. Holm,
 *    *How to mesh up Ewald sums I + II*,
 *    J. Chem. Phys. (109) 7678, 1998; (109) 7694, 1998
 *  - M. Deserno, C. Holm and H. J. Limbach,
 *    *How to mesh up Ewald sums*,
 *    in Molecular Dynamics on Parallel Computers,
 *    Ed. R. Esser et al., World Scientific, Singapore, 2000
 *  - M. Deserno,
 *    *Counterion condensation for rigid linear polyelectrolytes*,
 *    PhD Thesis, Universität Mainz, 2000
 *  - J. J. Cerda,
 *    *P3M for dipolar interactions*,
 *    J. Chem. Phys (129) 234104, 2008
 *
 *  Implementation in p3m.cpp.
 */

#include "config.hpp"

#ifdef P3M

#include "debug.hpp"
#include "fft.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "p3m-common.hpp"

#include "utils/math/AS_erfc_part.hpp"

/************************************************
 * data types
 ************************************************/

typedef struct {
  p3m_parameter_struct params;

  /** local mesh. */
  p3m_local_mesh local_mesh;
  /** real space mesh (local) for CA/FFT.*/
  double *rs_mesh;
  /** k-space mesh (local) for k-space calculation and FFT.*/
  double *ks_mesh;

  /** number of charged particles (only on master node). */
  int sum_qpart;
  /** Sum of square of charges (only on master node). */
  double sum_q2;
  /** square of sum of charges (only on master node). */
  double square_sum_q;

  /** interpolation of the charge assignment function. */
  double *int_caf[7];

  /** position shift for calc. of first assignment mesh point. */
  double pos_shift;
  /** help variable for calculation of aliasing sums */
  double *meshift_x;
  double *meshift_y;
  double *meshift_z;

  /** Spatial differential operator in k-space. We use an i*k differentiation.
   */
  double *d_op[3];
  /** Force optimised influence function (k-space) */
  double *g_force;
  /** Energy optimised influence function (k-space) */
  double *g_energy;

#ifdef P3M_STORE_CA_FRAC
  /** number of charged particles on the node. */
  int ca_num;
  /** Charge fractions for mesh assignment. */
  double *ca_frac;
  /** index of first mesh point for charge assignment. */
  int *ca_fmp;
#endif

  /** number of permutations in k_space */
  int ks_pnum;

  /** send/recv mesh sizes */
  p3m_send_mesh sm;

  /** Field to store grid points to send. */
  double *send_grid;
  /** Field to store grid points to recv */
  double *recv_grid;

  fft_data_struct fft;
} p3m_data_struct;

/** P3M parameters. */
extern p3m_data_struct p3m;

void p3m_pre_init(void);

/** Tune P3M parameters to desired accuracy.
 *
 *  The parameters
 *  @ref p3m_parameter_struct::mesh "mesh",
 *  @ref p3m_parameter_struct::cao "cao",
 *  @ref p3m_parameter_struct::r_cut_iL "r_cut_iL" and
 *  @ref p3m_parameter_struct::alpha_L "alpha_L"
 *  are tuned to obtain the target accuracy (initially stored in
 *  @ref p3m_parameter_struct::accuracy "accuracy") in optimal time.
 *  These parameters are stored in the @ref p3m object.
 *
 *  The function utilizes the analytic expression of the error estimate
 *  for the P3M method in the book of Hockney and Eastwood (Eqn. 8.23) in
 *  order to obtain the rms error in the force for a system of N randomly
 *  distributed particles in a cubic box.
 *  For the real space error the estimate of Kolafa/Perram is used.
 *
 *  Parameter ranges if not given explicit values via p3m_set_tune_params():
 *  - @p r_cut_iL starts from (@ref min_local_box_l - @ref #skin) / (
 *    n * @ref box_l), with n an integer (this implies @p r_cut_iL is the
 *    largest cutoff in the system!)
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
void p3m_init(void);

/** Update @ref p3m_parameter_struct::alpha "alpha" and
 *  @ref p3m_parameter_struct::r_cut "r_cut" if @ref box_l changed
 */
void p3m_scaleby_box_l();

/** Compute the k-space part of forces and energies for the charge-charge
 *  interaction
 */
double p3m_calc_kspace_forces(int force_flag, int energy_flag);

/** Compute the k-space part of the stress tensor **/
void p3m_calc_kspace_stress(double *stress);

/** Sanity checks */
bool p3m_sanity_checks();

/** Calculate number of charged particles, the sum of the squared
 *  charges and the squared sum of the charges.
 */
void p3m_count_charged_particles();

/** Assign the physical charges using the tabulated charge assignment function.
 *  If @ref P3M_STORE_CA_FRAC is true, then the charge fractions are buffered
 *  in @ref p3m_data_struct::ca_fmp "ca_fmp" and @ref p3m_data_struct::ca_frac
 *  "ca_frac".
 */
void p3m_charge_assign();

/** Assign a single charge into the current charge grid.
 *
 *  @param[in] q          %Particle charge
 *  @param[in] real_pos   %Particle position in real space
 *  @param[in] cp_cnt     The running index, which may be smaller than 0, in
 *                        which case the charge is assumed to be virtual and
 *                        is not stored in the @ref p3m_data_struct::ca_frac
 *                        "ca_frac" arrays
 */
void p3m_assign_charge(double q, Vector3d &real_pos, int cp_cnt);

/** Shrink wrap the charge grid */
void p3m_shrink_wrap_charge_grid(int n_charges);

/** Calculate real space contribution of Coulomb pair forces.
 *
 *  If NPT is compiled in, it returns the energy, which is needed for NPT.
 */
inline double p3m_add_pair_force(double chgfac, double const *d, double dist2,
                                 double dist, double force[3]) {
  if (dist < p3m.params.r_cut) {
    if (dist > 0.0) { // Vincent
      double adist = p3m.params.alpha * dist;
#if USE_ERFC_APPROXIMATION
      double erfc_part_ri = Utils::AS_erfc_part(adist) / dist;
      double fac1 = coulomb.prefactor * chgfac * exp(-adist * adist);
      double fac2 =
          fac1 * (erfc_part_ri + 2.0 * p3m.params.alpha * wupii) / dist2;
#else
      erfc_part_ri = erfc(adist) / dist;
      double fac1 = coulomb.prefactor * chgfac;
      double fac2 = fac1 *
                    (erfc_part_ri +
                     2.0 * p3m.params.alpha * wupii * exp(-adist * adist)) /
                    dist2;
#endif
      for (int j = 0; j < 3; j++)
        force[j] += fac2 * d[j];
      ESR_TRACE(
          fprintf(stderr, "%d: RSE: Pair dist=%.3f: force (%.3e,%.3e,%.3e)\n",
                  this_node, dist, fac2 * d[0], fac2 * d[1], fac2 * d[2]));
#ifdef NPT
      return fac1 * erfc_part_ri;
#endif
    }
  }
  return 0.0;
}

/** Set initial values for p3m_adaptive_tune()
 *
 *  @param[in]  r_cut        @copybrief p3m_parameter_struct::r_cut
 *  @param[in]  mesh         @copybrief p3m_parameter_struct::mesh
 *  @param[in]  cao          @copybrief p3m_parameter_struct::cao
 *  @param[in]  alpha        @copybrief p3m_parameter_struct::alpha
 *  @param[in]  accuracy     @copybrief p3m_parameter_struct::accuracy
 *  @param[in]  n_interpol   @copybrief p3m_parameter_struct::inter
 */
void p3m_set_tune_params(double r_cut, const int mesh[3], int cao, double alpha,
                         double accuracy, int n_interpol);

/** Set custom parameters
 *
 *  @param[in]  r_cut        @copybrief p3m_parameter_struct::r_cut
 *  @param[in]  mesh         @copybrief p3m_parameter_struct::mesh
 *  @param[in]  cao          @copybrief p3m_parameter_struct::cao
 *  @param[in]  alpha        @copybrief p3m_parameter_struct::alpha
 *  @param[in]  accuracy     @copybrief p3m_parameter_struct::accuracy
 *  @return Custom error code
 */
int p3m_set_params(double r_cut, const int *mesh, int cao, double alpha,
                   double accuracy);

/** Set mesh offset
 *
 *  @param[in]  x , y , z  Components of @ref p3m_parameter_struct::mesh_off
 *                         "mesh_off"
 */
int p3m_set_mesh_offset(double x, double y, double z);

/** Set @ref p3m_parameter_struct::epsilon "epsilon" parameter
 *
 *  @param[in]  eps          @copybrief p3m_parameter_struct::epsilon
 */
int p3m_set_eps(double eps);

/** Set @ref p3m_parameter_struct::inter "inter" parameter
 *
 *  @param[in]  n            @copybrief p3m_parameter_struct::inter
 */
int p3m_set_ninterpol(int n);

/** Calculate real space contribution of Coulomb pair energy. */
inline double p3m_pair_energy(double chgfac, double dist) {
  if (dist < p3m.params.r_cut && dist != 0) {
    double adist = p3m.params.alpha * dist;
#if USE_ERFC_APPROXIMATION
    double erfc_part_ri = Utils::AS_erfc_part(adist) / dist;
    return coulomb.prefactor * chgfac * erfc_part_ri * exp(-adist * adist);
#else
    double erfc_part_ri = erfc(adist) / dist;
    return coulomb.prefactor * chgfac * erfc_part_ri;
#endif
  }
  return 0.0;
}

/** Clean up P3M memory allocations. */
void p3m_free();

#endif /* of ifdef P3M */

#endif /*of ifndef P3M_H */
