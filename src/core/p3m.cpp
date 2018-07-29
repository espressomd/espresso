/*
  Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
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
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <mpi.h>

#include "cells.hpp"
#include "communication.hpp"
#include "domain_decomposition.hpp"
#include "elc.hpp"
#include "fft.hpp"
#include "global.hpp"
#include "grid.hpp"
#include "initialize.hpp"
#include "integrate.hpp"
#include "p3m.hpp"
#include "particle_data.hpp"
#include "thermostat.hpp"
#include "tuning.hpp"
#include "utils.hpp"
#ifdef CUDA
#include "p3m_gpu_error.hpp"
#endif

#ifdef P3M

/************************************************
 * variables
 ************************************************/
p3m_data_struct p3m;

/* MPI tags for the charge-charge p3m communications: */
/** Tag for communication in P3M_init() -> send_calc_mesh(). */
#define REQ_P3M_INIT 200
/** Tag for communication in p3m_gather_fft_grid(). */
#define REQ_P3M_GATHER 201
/** Tag for communication in p3m_spread_force_grid(). */
#define REQ_P3M_SPREAD 202

/* Index helpers for direct and reciprocal space
 * After the FFT the data is in order YZX, which
 * means that Y is the slowest changing index.
 * The defines are here to not get confused and
 * be able to easily change the order.
 */
#define RX 0
#define RY 1
#define RZ 2
#define KY 0
#define KZ 1
#define KX 2

/** \name Private Functions */
/************************************************************/
/*@{*/

#ifdef P3M_DEBUG
static void p3m_print(void) {
  fprintf(stderr,
          "general information: \n\t node: %d \n\t box_l: (%lf, %lf, %lf)\n",
          this_node, box_l[0], box_l[1], box_l[2]);

  fprintf(stderr, "p3m parameters:\n\t alpha_L: %lf\n\t r_cut_iL: %lf\n\t \
                   mesh: (%d, %d, %d)\n\t mesh_off: (%lf, %lf, %lf)\n\t \
                   cao: %d\n\t inter: %d\n\t accuracy: %lf\n\t epsilon: %lf\n\t \
                   cao_cut: (%lf, %lf, %lf)\n\t a: (%lf,%lf,%lf)\n\t \
                   ai: (%lf,%lf,%lf)\n\t alpha: %lf\n\t r_cut: %lf\n\t \
                   inter2: %d\n\t cao3: %d\n\t additional_mesh: (%lf,%lf,%lf)\n",
          p3m.params.alpha_L, p3m.params.r_cut_iL, p3m.params.mesh[0],
          p3m.params.mesh[1], p3m.params.mesh[2], p3m.params.mesh_off[0],
          p3m.params.mesh_off[1], p3m.params.mesh_off[2], p3m.params.cao,
          p3m.params.inter, p3m.params.accuracy, p3m.params.epsilon,
          p3m.params.cao_cut[0], p3m.params.cao_cut[1], p3m.params.cao_cut[2],
          p3m.params.a[0], p3m.params.a[1], p3m.params.a[2], p3m.params.ai[0],
          p3m.params.ai[1], p3m.params.ai[2], p3m.params.alpha,
          p3m.params.r_cut, p3m.params.inter2, p3m.params.cao3,
          p3m.params.additional_mesh[0], p3m.params.additional_mesh[1],
          p3m.params.additional_mesh[2]);
}

#endif

/** Calculates for charges the properties of the send/recv sub-meshes of the
 * local FFT mesh.
 *  In order to calculate the recv sub-meshes there is a communication of
 *  the margins between neighbouring nodes. */
static void p3m_calc_send_mesh();

/** Initializes the (inverse) mesh constant \ref p3m_parameter_struct::a (\ref
    p3m_parameter_struct::ai) and the cutoff for charge assignment \ref
    p3m_parameter_struct::cao_cut, which has to be done by \ref p3m_init
    once and by \ref p3m_scaleby_box_l whenever the \ref box_l
    changed.  */
static void p3m_init_a_ai_cao_cut(void);

/** Calculate the spacial position of the left down mesh point of the local
   mesh, to be
    stored in \ref p3m_local_mesh::ld_pos; function called by \ref
   p3m_calc_local_ca_mesh once
    and by \ref p3m_scaleby_box_l whenever the \ref box_l changed. */
static void p3m_calc_lm_ld_pos(void);

/** Calculates the dipole term */
static double p3m_calc_dipole_term(int force_flag, int energy_flag);

/** Gather FFT grid.
 *  After the charge assignment Each node needs to gather the
 *  information for the FFT grid in his spatial domain.
 */
static void p3m_gather_fft_grid(double *mesh);

/** Spread force grid.
 *  After the k-space calculations each node needs to get all force
 *  information to reassigne the forces from the grid to the
 *  particles.
 */
static void p3m_spread_force_grid(double *mesh);

#ifdef P3M_STORE_CA_FRAC
/** realloc charge assignment fields. */
static void p3m_realloc_ca_fields(int newsize);
#endif

static int p3m_sanity_checks_system(void);

/** checks for correctness for charges in P3M of the cao_cut, necessary when the
 * box length changes */
static int p3m_sanity_checks_boxl(void);

/** Calculate the spacial position of the left down mesh point of the local
   mesh, to be
    stored in \ref p3m_local_mesh::ld_pos; function called by \ref
   p3m_calc_local_ca_mesh once
    and by \ref p3m_scaleby_box_l whenever the \ref box_l changed. */
static void p3m_calc_lm_ld_pos(void);

/** Calculates properties of the local FFT mesh for the
    charge assignment process. */
static void p3m_calc_local_ca_mesh(void);

/** Interpolates the P-th order charge assignment function from
 * Hockney/Eastwood 5-189 (or 8-61). The following charge fractions
 * are also tabulated in Deserno/Holm. */
static void p3m_interpolate_charge_assignment_function(void);

/** shifts the mesh points by mesh/2 */
static void p3m_calc_meshift(void);

/** Calculates the Fourier transformed differential operator.
 *  Remark: This is done on the level of n-vectors and not k-vectors,
 *           i.e. the prefactor i*2*PI/L is missing! */
static void p3m_calc_differential_operator(void);

/** Calculates the optimal influence function of Hockney and Eastwood.
 * (optimised for force calculations)
 *
 *  Each node calculates only the values for its domain in k-space
 *  (see fft.plan[3].mesh and fft.plan[3].start).
 *
 *  See also: Hockney/Eastwood 8-22 (p275). Note the somewhat
 *  different convention for the prefactors, which is described in
 *  Deserno/Holm. */
static void p3m_calc_influence_function_force(void);

/** Calculates the influence function optimized for the energy and the
    self energy correction.  */
static void p3m_calc_influence_function_energy(void);

/** Calculates the aliasing sums for the optimal influence function.
 *
 * Calculates the aliasing sums in the nominator and denominator of
 * the expression for the optimal influence function (see
 * Hockney/Eastwood: 8-22, p. 275).
 *
 * \param  n           n-vector for which the aliasing sum is to be performed.
 * \param  nominator   aliasing sums in the nominator.
 * \return denominator aliasing sum in the denominator
 */
double p3m_perform_aliasing_sums_force(int n[3], double nominator[3]);
double p3m_perform_aliasing_sums_energy(int n[3]);

/*@}*/

/** \name P3M Tuning Functions (private)*/
/************************************************************/
/*@{*/

/** Calculates the real space contribution to the rms error in the force (as
   described
   by Kolafa and Perram).
   \param prefac   Prefactor of coulomb interaction.
   \param r_cut_iL rescaled real space cutoff for p3m method.
   \param n_c_part number of charged particles in the system.
   \param sum_q2   sum of square of charges in the system
   \param alpha_L  rescaled ewald splitting parameter.
   \return real space error
*/
static double p3m_real_space_error(double prefac, double r_cut_iL, int n_c_part,
                                   double sum_q2, double alpha_L);

/** Calculate the analytic expression of the error estimate for the
    P3M method in the book of Hockney and Eastwood (Eqn. 8.23) in
    order to obtain the rms error in the force for a system of N
    randomly distributed particles in a cubic box (k space part).
    \param prefac   Prefactor of coulomb interaction.
    \param mesh     number of mesh points in one direction.
    \param cao      charge assignment order.
    \param n_c_part number of charged particles in the system.
    \param sum_q2   sum of square of charges in the system
    \param alpha_L  rescaled ewald splitting parameter.
    \return reciprocal (k) space error
*/
static double p3m_k_space_error(double prefac, int mesh[3], int cao,
                                int n_c_part, double sum_q2, double alpha_L);

/** aliasing sum used by \ref p3m_k_space_error. */
static void p3m_tune_aliasing_sums(int nx, int ny, int nz, int mesh[3],
                                   double mesh_i[3], int cao, double alpha_L_i,
                                   double *alias1, double *alias2);

/** Template parameterized calculation of the charge assignment to be called by
   wrapper.
    \param cao      charge assignment order.
*/
template <int cao> static void p3m_do_charge_assign();

template <int cao>
void p3m_do_assign_charge(double q, Vector3d& real_pos, int cp_cnt);
/*@}*/

void p3m_pre_init(void) {
  p3m_common_parameter_pre_init(&p3m.params);
  /* p3m.local_mesh is uninitialized */
  /* p3m.sm is uninitialized */

  p3m.rs_mesh = nullptr;
  p3m.ks_mesh = nullptr;
  p3m.sum_qpart = 0;
  p3m.sum_q2 = 0.0;
  p3m.square_sum_q = 0.0;

  for (int i = 0; i < 7; i++)
    p3m.int_caf[i] = nullptr;
  p3m.pos_shift = 0.0;
  p3m.meshift_x = nullptr;
  p3m.meshift_y = nullptr;
  p3m.meshift_z = nullptr;

  p3m.d_op[0] = nullptr;
  p3m.d_op[1] = nullptr;
  p3m.d_op[2] = nullptr;
  p3m.g_force = nullptr;
  p3m.g_energy = nullptr;

#ifdef P3M_STORE_CA_FRAC
  p3m.ca_num = 0;
  p3m.ca_frac = nullptr;
  p3m.ca_fmp = nullptr;
#endif
  p3m.ks_pnum = 0;

  p3m.send_grid = nullptr;
  p3m.recv_grid = nullptr;

  fft_pre_init();
}

void p3m_free() {
  int i;
/* free memory */
#ifdef P3M_STORE_CA_FRAC
  free(p3m.ca_frac);
  free(p3m.ca_fmp);
#endif
  free(p3m.send_grid);
  free(p3m.recv_grid);
  free(p3m.rs_mesh);
  free(p3m.ks_mesh);
  for (i = 0; i < p3m.params.cao; i++)
    free(p3m.int_caf[i]);
}

void p3m_set_prefactor() {
  p3m.params.alpha = 0.0;
  p3m.params.alpha_L = 0.0;
  p3m.params.r_cut = 0.0;
  p3m.params.r_cut_iL = 0.0;
  p3m.params.mesh[0] = 0;
  p3m.params.mesh[1] = 0;
  p3m.params.mesh[2] = 0;
  p3m.params.cao = 0;
}

void p3m_init() {
  if (coulomb.prefactor <= 0.0) {
    p3m.params.r_cut = 0.0;
    p3m.params.r_cut_iL = 0.0;

    if (this_node == 0) {
      P3M_TRACE(fprintf(stderr, "0: P3M_init: prefactor is "
                                "zero.\nElectrostatics switched off!\n"););
    }

  } else {
    P3M_TRACE(fprintf(stderr, "%d: p3m_init: \n", this_node));

    if (p3m_sanity_checks()) {
      return;
    }

    P3M_TRACE(fprintf(stderr, "%d: p3m_init: starting\n", this_node));

    P3M_TRACE(fprintf(stderr, "%d: mesh=%d, cao=%d, mesh_off=(%f,%f,%f)\n",
                      this_node, p3m.params.mesh[0], p3m.params.cao,
                      p3m.params.mesh_off[0], p3m.params.mesh_off[1],
                      p3m.params.mesh_off[2]));
    p3m.params.cao3 = p3m.params.cao * p3m.params.cao * p3m.params.cao;

    /* initializes the (inverse) mesh constant p3m.params.a (p3m.params.ai) and
     * the cutoff for charge assignment p3m.params.cao_cut */
    p3m_init_a_ai_cao_cut();

#ifdef P3M_STORE_CA_FRAC
    /* initialize ca fields to size CA_INCREMENT: p3m.ca_frac and p3m.ca_fmp */
    p3m.ca_num = 0;
    p3m_realloc_ca_fields(CA_INCREMENT);
#endif

    p3m_calc_local_ca_mesh();

    p3m_calc_send_mesh();
    P3M_TRACE(p3m_p3m_print_local_mesh(p3m.local_mesh));
    P3M_TRACE(p3m_p3m_print_send_mesh(p3m.sm));
    p3m.send_grid = Utils::realloc(p3m.send_grid, sizeof(double) * p3m.sm.max);
    p3m.recv_grid = Utils::realloc(p3m.recv_grid, sizeof(double) * p3m.sm.max);

    /* fix box length dependent constants */
    p3m_scaleby_box_l();

    if (p3m.params.inter > 0)
      p3m_interpolate_charge_assignment_function();

    /* position offset for calc. of first meshpoint */
    p3m.pos_shift =
        (double)((p3m.params.cao - 1) / 2) - (p3m.params.cao % 2) / 2.0;
    P3M_TRACE(
        fprintf(stderr, "%d: p3m.pos_shift = %f\n", this_node, p3m.pos_shift));

    /* FFT */
    P3M_TRACE(
        fprintf(stderr, "%d: p3m.rs_mesh ADR=%p\n", this_node, (void*) p3m.rs_mesh));

    int ca_mesh_size =
        fft_init(&p3m.rs_mesh, p3m.local_mesh.dim, p3m.local_mesh.margin,
                 p3m.params.mesh, p3m.params.mesh_off, &p3m.ks_pnum);
    p3m.ks_mesh = Utils::realloc(p3m.ks_mesh, ca_mesh_size * sizeof(double));

    P3M_TRACE(
        fprintf(stderr, "%d: p3m.rs_mesh ADR=%p\n", this_node, (void*) p3m.rs_mesh));

    /* k-space part: */
    p3m_calc_differential_operator();

    p3m_count_charged_particles();

    P3M_TRACE(fprintf(stderr, "%d: p3m-charges  initialized\n", this_node));
  }
}

void p3m_set_tune_params(double r_cut, int mesh[3], int cao, double alpha,
                         double accuracy, int n_interpol) {
  if (r_cut >= 0) {
    p3m.params.r_cut = r_cut;
    p3m.params.r_cut_iL = r_cut * box_l_i[0];
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
    p3m.params.alpha_L = alpha * box_l[0];
  }

  if (accuracy >= 0)
    p3m.params.accuracy = accuracy;

  if (n_interpol != -1)
    p3m.params.inter = n_interpol;

}

/*@}*/

int p3m_set_params(double r_cut, int *mesh, int cao, double alpha,
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
  p3m.params.r_cut_iL = r_cut * box_l_i[0];
  p3m.params.mesh[2] = mesh[2];
  p3m.params.mesh[1] = mesh[1];
  p3m.params.mesh[0] = mesh[0];
  p3m.params.cao = cao;

  if (alpha > 0) {
    p3m.params.alpha = alpha;
    p3m.params.alpha_L = alpha * box_l[0];
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

int p3m_set_ninterpol(int n) {
  if (n < 0)
    return ES_ERROR;

  p3m.params.inter = n;

  mpi_bcast_coulomb_params();

  return ES_OK;
}

/************************************* method ********************************/
/*****************************************************************************/

void p3m_interpolate_charge_assignment_function() {
  double dInterpol = 0.5 / (double)p3m.params.inter;
  int i;
  long j;

  if (p3m.params.inter == 0)
    return;

  P3M_TRACE(fprintf(
      stderr,
      "%d - interpolating (%d) the order-%d charge assignment function\n",
      this_node, p3m.params.inter, p3m.params.cao));

  p3m.params.inter2 = 2 * p3m.params.inter + 1;

  for (i = 0; i < p3m.params.cao; i++) {
    /* allocate memory for interpolation array */
    p3m.int_caf[i] = Utils::realloc(
        p3m.int_caf[i], sizeof(double) * (2 * p3m.params.inter + 1));

    /* loop over all interpolation points */
    for (j = -p3m.params.inter; j <= p3m.params.inter; j++)
      p3m.int_caf[i][j + p3m.params.inter] =
          p3m_caf(i, j * dInterpol, p3m.params.cao);
  }
}

/* Template wrapper for p3m_do_charge_assign() */
void p3m_charge_assign() {
  switch (p3m.params.cao) {
  case 1:
    p3m_do_charge_assign<1>();
    break;
  case 2:
    p3m_do_charge_assign<2>();
    break;
  case 3:
    p3m_do_charge_assign<3>();
    break;
  case 4:
    p3m_do_charge_assign<4>();
    break;
  case 5:
    p3m_do_charge_assign<5>();
    break;
  case 6:
    p3m_do_charge_assign<6>();
    break;
  case 7:
    p3m_do_charge_assign<7>();
    break;
  }
}

/* assign the charges */
template <int cao> void p3m_do_charge_assign() {
  /* charged particle counter, charge fraction counter */
  int cp_cnt = 0;
  /* prepare local FFT mesh */
  for (int i = 0; i < p3m.local_mesh.size; i++)
    p3m.rs_mesh[i] = 0.0;

  for (auto &p : local_cells.particles()) {
    if (p.p.q != 0.0) {
      p3m_do_assign_charge<cao>(p.p.q, p.r.p, cp_cnt);
      cp_cnt++;
    }
  }

#ifdef P3M_STORE_CA_FRAC
  p3m_shrink_wrap_charge_grid(cp_cnt);
#endif
}

/* Template wrapper for p3m_do_assign_charge() */
void p3m_assign_charge(double q, Vector3d& real_pos, int cp_cnt) {
  switch (p3m.params.cao) {
  case 1:
    p3m_do_assign_charge<1>(q, real_pos, cp_cnt);
    break;
  case 2:
    p3m_do_assign_charge<2>(q, real_pos, cp_cnt);
    break;
  case 3:
    p3m_do_assign_charge<3>(q, real_pos, cp_cnt);
    break;
  case 4:
    p3m_do_assign_charge<4>(q, real_pos, cp_cnt);
    break;
  case 5:
    p3m_do_assign_charge<5>(q, real_pos, cp_cnt);
    break;
  case 6:
    p3m_do_assign_charge<6>(q, real_pos, cp_cnt);
    break;
  case 7:
    p3m_do_assign_charge<7>(q, real_pos, cp_cnt);
    break;
  }
}

template <int cao>
void p3m_do_assign_charge(double q, Vector3d& real_pos, int cp_cnt) {
  auto const inter = not (p3m.params.inter == 0);
  /* distance to nearest mesh point */
  double dist[3];
  /* index for caf interpolation grid */
  int arg[3];
  /* index, index jumps for rs_mesh array */
  int q_ind = 0;
#ifdef P3M_STORE_CA_FRAC
  // make sure we have enough space
  if (cp_cnt >= p3m.ca_num)
    p3m_realloc_ca_fields(cp_cnt + 1);
  // do it here, since p3m_realloc_ca_fields may change the address of
  // p3m.ca_frac
  double *cur_ca_frac = p3m.ca_frac + cao * cao * cao * cp_cnt;
#endif

  for (int d = 0; d < 3; d++) {
    /* particle position in mesh coordinates */
    auto const pos = ((real_pos[d] - p3m.local_mesh.ld_pos[d]) * p3m.params.ai[d]) -
          p3m.pos_shift;
    /* nearest mesh point */
    auto const nmp = (int)pos;
    /* 3d-array index of nearest mesh point */
    q_ind = (d == 0) ? nmp : nmp + p3m.local_mesh.dim[d] * q_ind;

    if (!inter)
      /* distance to nearest mesh point */
      dist[d] = (pos - nmp) - 0.5;
    else
      /* distance to nearest mesh point for interpolation */
      arg[d] = (int)((pos - nmp) * p3m.params.inter2);

#ifdef ADDITIONAL_CHECKS
    if (pos < -skin * p3m.params.ai[d]) {
      fprintf(stderr, "%d: rs_mesh underflow! (pos %f)\n", this_node,
              real_pos[d]);
      fprintf(stderr, "%d: allowed coordinates: %f - %f\n", this_node,
              my_left[d] - skin, my_right[d] + skin);
    }
    if ((nmp + cao) > p3m.local_mesh.dim[d]) {
      fprintf(stderr, "%d: rs_mesh overflow! (pos %f, nmp=%d)\n", this_node,
              real_pos[d], nmp);
      fprintf(stderr, "%d: allowed coordinates: %f - %f\n", this_node,
              my_left[d] - skin, my_right[d] + skin);
    }
#endif
  }

#ifdef P3M_STORE_CA_FRAC
  if (cp_cnt >= 0)
    p3m.ca_fmp[cp_cnt] = q_ind;
#endif

  if (!inter) {
    for (int i0 = 0; i0 < cao; i0++) {
      auto const tmp0 = p3m_caf(i0, dist[0], cao);
      for (int i1 = 0; i1 < cao; i1++) {
        auto const tmp1 = tmp0 * p3m_caf(i1, dist[1], cao);
        for (int i2 = 0; i2 < cao; i2++) {
          auto const cur_ca_frac_val = q * tmp1 * p3m_caf(i2, dist[2], cao);
          p3m.rs_mesh[q_ind] += cur_ca_frac_val;
#ifdef P3M_STORE_CA_FRAC
          /* store current ca frac */
          if (cp_cnt >= 0)
            *(cur_ca_frac++) = cur_ca_frac_val;
#endif
          q_ind++;
        }
        q_ind += p3m.local_mesh.q_2_off;
      }
      q_ind += p3m.local_mesh.q_21_off;
    }
  } else {
    for (int i0 = 0; i0 < cao; i0++) {
      auto const tmp0 = p3m.int_caf[i0][arg[0]];
      for (int i1 = 0; i1 < cao; i1++) {
        auto const tmp1 = tmp0 * p3m.int_caf[i1][arg[1]];
        for (int i2 = 0; i2 < cao; i2++) {
          auto const cur_ca_frac_val = q * tmp1 * p3m.int_caf[i2][arg[2]];
          p3m.rs_mesh[q_ind] += cur_ca_frac_val;
#ifdef P3M_STORE_CA_FRAC
          /* store current ca frac */
          if (cp_cnt >= 0)
            *(cur_ca_frac++) = cur_ca_frac_val;
#endif
          q_ind++;
        }
        q_ind += p3m.local_mesh.q_2_off;
      }
      q_ind += p3m.local_mesh.q_21_off;
    }
  }
}

#ifdef P3M_STORE_CA_FRAC
/** shrink wrap the charge grid */
void p3m_shrink_wrap_charge_grid(int n_charges) {
  /* we do not really want to export these */
  if (n_charges < p3m.ca_num)
    p3m_realloc_ca_fields(n_charges);
}
#endif

/* assign the forces obtained from k-space */
template <int cao>
static void P3M_assign_forces(double force_prefac, int d_rs) {
  /* charged particle counter, charge fraction counter */
  int cp_cnt = 0;
#ifdef P3M_STORE_CA_FRAC
  int cf_cnt = 0;
#else
  /* distance to nearest mesh point */
  double dist[3];
  /* index for caf interpolation grid */
  int arg[3];
#endif
  /* index, index jumps for rs_mesh array */
  int q_ind = 0;

  for (auto &p : local_cells.particles()) {
    auto const q = p.p.q;
    if (q != 0.0) {
#ifdef P3M_STORE_CA_FRAC
      q_ind = p3m.ca_fmp[cp_cnt];
      for (int i0 = 0; i0 < cao; i0++) {
        for (int i1 = 0; i1 < cao; i1++) {
          for (int i2 = 0; i2 < cao; i2++) {
            p.f.f[d_rs] -=
                force_prefac * p3m.ca_frac[cf_cnt] * p3m.rs_mesh[q_ind];
            q_ind++;
            cf_cnt++;
          }
          q_ind += p3m.local_mesh.q_2_off;
        }
        q_ind += p3m.local_mesh.q_21_off;
      }
      cp_cnt++;
#else
      double pos;
      int nmp;
      double tmp0, tmp1;
      double cur_ca_frac_val;
      for (int d = 0; d < 3; d++) {
        /* particle position in mesh coordinates */
        pos = ((p.r.p[d] - p3m.local_mesh.ld_pos[d]) * p3m.params.ai[d]) -
              p3m.pos_shift;
        /* nearest mesh point */
        nmp = (int)pos;
        /* 3d-array index of nearest mesh point */
        q_ind = (d == 0) ? nmp : nmp + p3m.local_mesh.dim[d] * q_ind;

        if (p3m.params.inter == 0)
          /* distance to nearest mesh point */
          dist[d] = (pos - nmp) - 0.5;
        else
          /* distance to nearest mesh point for interpolation */
          arg[d] = (int)((pos - nmp) * p3m.params.inter2);
      }

      if (p3m.params.inter == 0) {
        for (int i0 = 0; i0 < cao; i0++) {
          tmp0 = p3m_caf(i0, dist[0], cao);
          for (int i1 = 0; i1 < cao; i1++) {
            tmp1 = tmp0 * p3m_caf(i1, dist[1], cao);
            for (int i2 = 0; i2 < cao; i2++) {
              cur_ca_frac_val = q * tmp1 * p3m_caf(i2, dist[2], cao);
              p.f.f[d_rs] -=
                  force_prefac * cur_ca_frac_val * p3m.rs_mesh[q_ind];
              q_ind++;
            }
            q_ind += p3m.local_mesh.q_2_off;
          }
          q_ind += p3m.local_mesh.q_21_off;
        }
      } else {
        for (int i0 = 0; i0 < cao; i0++) {
          tmp0 = p3m.int_caf[i0][arg[0]];
          for (int i1 = 0; i1 < cao; i1++) {
            tmp1 = tmp0 * p3m.int_caf[i1][arg[1]];
            for (int i2 = 0; i2 < cao; i2++) {
              cur_ca_frac_val = q * tmp1 * p3m.int_caf[i2][arg[2]];
              p.f.f[d_rs] -=
                  force_prefac * cur_ca_frac_val * p3m.rs_mesh[q_ind];
              q_ind++;
            }
            q_ind += p3m.local_mesh.q_2_off;
          }
          q_ind += p3m.local_mesh.q_21_off;
        }
      }
#endif

      ONEPART_TRACE(if (p.p.identity == check_id) fprintf(
          stderr, "%d: OPT: P3M  f = (%.3e,%.3e,%.3e) in dir %d\n", this_node,
          p.f.f[0], p.f.f[1], p.f.f[2], d_rs));
    }
  }
}

double p3m_calc_kspace_forces(int force_flag, int energy_flag) {
  int i, d, d_rs, ind, j[3];
  /**************************************************************/
  /* Prefactor for force */
  double force_prefac;
  /* k space energy */
  double k_space_energy = 0.0, node_k_space_energy = 0.0;
  /* directions */
  double *d_operator = nullptr;

  P3M_TRACE(fprintf(stderr, "%d: p3m_perform: \n", this_node));
  //     fprintf(stderr, "calculating kspace forces\n");

  force_prefac = coulomb.prefactor / (2 * box_l[0] * box_l[1] * box_l[2]);

  /* Gather information for FFT grid inside the nodes domain (inner local mesh)
   */
  /* and Perform forward 3D FFT (Charge Assignment Mesh). */
  if (p3m.sum_q2 > 0) {
    p3m_gather_fft_grid(p3m.rs_mesh);
    fft_perform_forw(p3m.rs_mesh);
  }
  // Note: after these calls, the grids are in the order yzx and not xyz
  // anymore!!!

  /* === K Space Calculations === */
  P3M_TRACE(fprintf(stderr, "%d: p3m_perform: k-Space\n", this_node));

  /* === K Space Energy Calculation  === */
  //     if(energy_flag && p3m.sum_q2 > 0) {
  if (energy_flag) {
    /*********************
    Coulomb energy
    **********************/

    for (i = 0; i < fft.plan[3].new_size; i++) {
      // Use the energy optimized influence function for energy!
      node_k_space_energy += p3m.g_energy[i] * (Utils::sqr(p3m.rs_mesh[2 * i]) +
                                                Utils::sqr(p3m.rs_mesh[2 * i + 1]));
    }
    node_k_space_energy *= force_prefac;

    MPI_Reduce(&node_k_space_energy, &k_space_energy, 1, MPI_DOUBLE, MPI_SUM, 0,
               comm_cart);
    if (this_node == 0) {
      /* self energy correction */
      k_space_energy -=
          coulomb.prefactor * (p3m.sum_q2 * p3m.params.alpha * wupii);
      /* net charge correction */
      k_space_energy -=
          coulomb.prefactor * p3m.square_sum_q * PI /
          (2.0 * box_l[0] * box_l[1] * box_l[2] * Utils::sqr(p3m.params.alpha));
    }

  } /* if (energy_flag) */

  /* === K Space Force Calculation  === */
  if (force_flag && p3m.sum_q2 > 0) {
    /***************************
     COULOMB FORCES (k-space)
     ****************************/
    /* Force preparation */
    ind = 0;
    /* apply the influence function */
    for (i = 0; i < fft.plan[3].new_size; i++) {
      p3m.ks_mesh[ind] = p3m.g_force[i] * p3m.rs_mesh[ind];
      ind++;
      p3m.ks_mesh[ind] = p3m.g_force[i] * p3m.rs_mesh[ind];
      ind++;
    }

    /* === 3 Fold backward 3D FFT (Force Component Meshs) === */

    /* Force component loop */
    for (d = 0; d < 3; d++) {
      if (d == KX)
        d_operator = p3m.d_op[RX];
      else if (d == KY)
        d_operator = p3m.d_op[RY];
      else if (d == KZ)
        d_operator = p3m.d_op[RZ];

      /* direction in k space: */
      d_rs = (d + p3m.ks_pnum) % 3;
      /* srqt(-1)*k differentiation */
      ind = 0;
      for (j[0] = 0; j[0] < fft.plan[3].new_mesh[0]; j[0]++) {
        for (j[1] = 0; j[1] < fft.plan[3].new_mesh[1]; j[1]++) {
          for (j[2] = 0; j[2] < fft.plan[3].new_mesh[2]; j[2]++) {
            /* i*k*(Re+i*Im) = - Im*k + i*Re*k     (i=sqrt(-1)) */
            p3m.rs_mesh[ind] = -2.0 * PI *
                               (p3m.ks_mesh[ind + 1] *
                                d_operator[j[d] + fft.plan[3].start[d]]) /
                               box_l[d_rs];
            ind++;
            p3m.rs_mesh[ind] = 2.0 * PI * p3m.ks_mesh[ind - 1] *
                               d_operator[j[d] + fft.plan[3].start[d]] /
                               box_l[d_rs];
            ind++;
          }
        }
      }
      /* Back FFT force component mesh */
      fft_perform_back(p3m.rs_mesh, /* check_complex */ !p3m.params.tuning);
      /* redistribute force component mesh */
      p3m_spread_force_grid(p3m.rs_mesh);
      /* Assign force component from mesh to particle */
      switch (p3m.params.cao) {
      case 1:
        P3M_assign_forces<1>(force_prefac, d_rs);
        break;
      case 2:
        P3M_assign_forces<2>(force_prefac, d_rs);
        break;
      case 3:
        P3M_assign_forces<3>(force_prefac, d_rs);
        break;
      case 4:
        P3M_assign_forces<4>(force_prefac, d_rs);
        break;
      case 5:
        P3M_assign_forces<5>(force_prefac, d_rs);
        break;
      case 6:
        P3M_assign_forces<6>(force_prefac, d_rs);
        break;
      case 7:
        P3M_assign_forces<7>(force_prefac, d_rs);
        break;
      }
    }
  } /* if(force_flag) */

  if (p3m.params.epsilon != P3M_EPSILON_METALLIC) {
    k_space_energy += p3m_calc_dipole_term(force_flag, energy_flag);
  }

  return k_space_energy;
}

/************************************************************/

double p3m_calc_dipole_term(int force_flag, int energy_flag) {
  double pref = coulomb.prefactor * 4 * M_PI * box_l_i[0] * box_l_i[1] *
                box_l_i[2] / (2 * p3m.params.epsilon + 1);
  double lcl_dm[3], gbl_dm[3];
  double en;

  for (int j = 0; j < 3; j++)
    lcl_dm[j] = 0;

  for (auto const &p : local_cells.particles()) {
    for (int j = 0; j < 3; j++)
      /* dipole moment with unfolded coordinates */
      lcl_dm[j] += p.p.q * (p.r.p[j] + p.l.i[j] * box_l[j]);
  }

  MPI_Allreduce(lcl_dm, gbl_dm, 3, MPI_DOUBLE, MPI_SUM, comm_cart);

  if (energy_flag)
    en = 0.5 * pref * (Utils::sqr(gbl_dm[0]) + Utils::sqr(gbl_dm[1]) + Utils::sqr(gbl_dm[2]));
  else
    en = 0;
  if (force_flag) {
    for (int j = 0; j < 3; j++)
      gbl_dm[j] *= pref;
    for (auto &p : local_cells.particles()) {
      for (int j = 0; j < 3; j++)
        p.f.f[j] -= gbl_dm[j] * p.p.q;
    }
    return en;
  }
  return 0;
}

/************************************************************/

void p3m_gather_fft_grid(double *themesh) {
  int s_dir, r_dir, evenodd;
  MPI_Status status;
  double *tmp_ptr;

  P3M_TRACE(fprintf(stderr, "%d: p3m_gather_fft_grid:\n", this_node));

  /* direction loop */
  for (s_dir = 0; s_dir < 6; s_dir++) {
    if (s_dir % 2 == 0)
      r_dir = s_dir + 1;
    else
      r_dir = s_dir - 1;
    /* pack send block */
    if (p3m.sm.s_size[s_dir] > 0)
      fft_pack_block(themesh, p3m.send_grid, p3m.sm.s_ld[s_dir],
                     p3m.sm.s_dim[s_dir], p3m.local_mesh.dim, 1);

    /* communication */
    if (node_neighbors[s_dir] != this_node) {
      for (evenodd = 0; evenodd < 2; evenodd++) {
        if ((node_pos[s_dir / 2] + evenodd) % 2 == 0) {
          if (p3m.sm.s_size[s_dir] > 0)
            MPI_Send(p3m.send_grid, p3m.sm.s_size[s_dir], MPI_DOUBLE,
                     node_neighbors[s_dir], REQ_P3M_GATHER, comm_cart);
        } else {
          if (p3m.sm.r_size[r_dir] > 0)
            MPI_Recv(p3m.recv_grid, p3m.sm.r_size[r_dir], MPI_DOUBLE,
                     node_neighbors[r_dir], REQ_P3M_GATHER, comm_cart, &status);
        }
      }
    } else {
      tmp_ptr = p3m.recv_grid;
      p3m.recv_grid = p3m.send_grid;
      p3m.send_grid = tmp_ptr;
    }
    /* add recv block */
    if (p3m.sm.r_size[r_dir] > 0) {
      p3m_add_block(p3m.recv_grid, themesh, p3m.sm.r_ld[r_dir],
                    p3m.sm.r_dim[r_dir], p3m.local_mesh.dim);
    }
  }
}

void p3m_spread_force_grid(double *themesh) {
  int s_dir, r_dir, evenodd;
  MPI_Status status;
  double *tmp_ptr;
  P3M_TRACE(fprintf(stderr, "%d: p3m_spread_force_grid:\n", this_node));

  /* direction loop */
  for (s_dir = 5; s_dir >= 0; s_dir--) {
    if (s_dir % 2 == 0)
      r_dir = s_dir + 1;
    else
      r_dir = s_dir - 1;
    /* pack send block */
    if (p3m.sm.s_size[s_dir] > 0)
      fft_pack_block(themesh, p3m.send_grid, p3m.sm.r_ld[r_dir],
                     p3m.sm.r_dim[r_dir], p3m.local_mesh.dim, 1);
    /* communication */
    if (node_neighbors[r_dir] != this_node) {
      for (evenodd = 0; evenodd < 2; evenodd++) {
        if ((node_pos[r_dir / 2] + evenodd) % 2 == 0) {
          if (p3m.sm.r_size[r_dir] > 0)
            MPI_Send(p3m.send_grid, p3m.sm.r_size[r_dir], MPI_DOUBLE,
                     node_neighbors[r_dir], REQ_P3M_SPREAD, comm_cart);
        } else {
          if (p3m.sm.s_size[s_dir] > 0)
            MPI_Recv(p3m.recv_grid, p3m.sm.s_size[s_dir], MPI_DOUBLE,
                     node_neighbors[s_dir], REQ_P3M_SPREAD, comm_cart, &status);
        }
      }
    } else {
      tmp_ptr = p3m.recv_grid;
      p3m.recv_grid = p3m.send_grid;
      p3m.send_grid = tmp_ptr;
    }
    /* un pack recv block */
    if (p3m.sm.s_size[s_dir] > 0) {
      fft_unpack_block(p3m.recv_grid, themesh, p3m.sm.s_ld[s_dir],
                       p3m.sm.s_dim[s_dir], p3m.local_mesh.dim, 1);
    }
  }
}

#ifdef P3M_STORE_CA_FRAC
void p3m_realloc_ca_fields(int newsize) {
  newsize = ((newsize + CA_INCREMENT - 1) / CA_INCREMENT) * CA_INCREMENT;
  if (newsize == p3m.ca_num)
    return;
  if (newsize < CA_INCREMENT)
    newsize = CA_INCREMENT;

  P3M_TRACE(fprintf(stderr,
                    "%d: p3m_realloc_ca_fields: old_size=%d -> new_size=%d\n",
                    this_node, p3m.ca_num, newsize));
  p3m.ca_num = newsize;
  p3m.ca_frac = Utils::realloc(p3m.ca_frac,
                               p3m.params.cao3 * p3m.ca_num * sizeof(double));
  p3m.ca_fmp = Utils::realloc(p3m.ca_fmp, p3m.ca_num * sizeof(int));
}
#endif

void p3m_calc_meshift(void) {
  int i;

  p3m.meshift_x =
      Utils::realloc(p3m.meshift_x, p3m.params.mesh[0] * sizeof(double));
  p3m.meshift_y =
      Utils::realloc(p3m.meshift_y, p3m.params.mesh[1] * sizeof(double));
  p3m.meshift_z =
      Utils::realloc(p3m.meshift_z, p3m.params.mesh[2] * sizeof(double));

  p3m.meshift_x[0] = p3m.meshift_y[0] = p3m.meshift_z[0] = 0;
  for (i = 1; i <= p3m.params.mesh[RX] / 2; i++) {
    p3m.meshift_x[i] = i;
    p3m.meshift_x[p3m.params.mesh[0] - i] = -i;
  }

  for (i = 1; i <= p3m.params.mesh[RY] / 2; i++) {
    p3m.meshift_y[i] = i;
    p3m.meshift_y[p3m.params.mesh[1] - i] = -i;
  }

  for (i = 1; i <= p3m.params.mesh[RZ] / 2; i++) {
    p3m.meshift_z[i] = i;
    p3m.meshift_z[p3m.params.mesh[2] - i] = -i;
  }
}

void p3m_calc_differential_operator() {
  int i, j;

  for (i = 0; i < 3; i++) {
    p3m.d_op[i] =
        Utils::realloc(p3m.d_op[i], p3m.params.mesh[i] * sizeof(double));
    p3m.d_op[i][0] = 0;
    p3m.d_op[i][p3m.params.mesh[i] / 2] = 0.0;

    for (j = 1; j < p3m.params.mesh[i] / 2; j++) {
      p3m.d_op[i][j] = j;
      p3m.d_op[i][p3m.params.mesh[i] - j] = -j;
    }
  }
}

namespace {

template <int cao>
inline double perform_aliasing_sums_force(int n[3], double numerator[3]) {
  using Utils::int_pow;

  int i;
  double denominator = 0.0;
  /* lots of temporary variables... */
  double sx, sy, sz, f1, f2, mx, my, mz, nmx, nmy, nmz, nm2, expo;
  double limit = 30;

  for (i = 0; i < 3; i++)
    numerator[i] = 0.0;

  f1 = Utils::sqr(PI / (p3m.params.alpha));

  for (mx = -P3M_BRILLOUIN; mx <= P3M_BRILLOUIN; mx++) {
    nmx = p3m.meshift_x[n[KX]] + p3m.params.mesh[RX] * mx;
    sx = int_pow<2 * cao>(sinc(nmx / (double)p3m.params.mesh[RX]));
    for (my = -P3M_BRILLOUIN; my <= P3M_BRILLOUIN; my++) {
      nmy = p3m.meshift_y[n[KY]] + p3m.params.mesh[RY] * my;
      sy = sx * int_pow<2 * cao>(sinc(nmy / (double)p3m.params.mesh[RY]));
      for (mz = -P3M_BRILLOUIN; mz <= P3M_BRILLOUIN; mz++) {
        nmz = p3m.meshift_z[n[KZ]] + p3m.params.mesh[RZ] * mz;
        sz = sy * int_pow<2 * cao>(sinc(nmz / (double)p3m.params.mesh[RZ]));

        nm2 =
            Utils::sqr(nmx / box_l[RX]) + Utils::sqr(nmy / box_l[RY]) + Utils::sqr(nmz / box_l[RZ]);
        expo = f1 * nm2;
        f2 = (expo < limit) ? sz * exp(-expo) / nm2 : 0.0;

        numerator[RX] += f2 * nmx / box_l[RX];
        numerator[RY] += f2 * nmy / box_l[RY];
        numerator[RZ] += f2 * nmz / box_l[RZ];

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
    size *= fft.plan[3].new_mesh[i];
    end[i] = fft.plan[3].start[i] + fft.plan[3].new_mesh[i];
  }

  p3m.g_force = Utils::realloc(p3m.g_force, size * sizeof(double));

  /* Skip influence function calculation in tuning mode,
     the results need not be correct for timing. */
  if(p3m.params.tuning) {
    /* If resized, fill with zeros to avoid nan forces. */
      memset(p3m.g_force, 0, size * sizeof(double));

    return;
  }

  for (n[0] = fft.plan[3].start[0]; n[0] < end[0]; n[0]++) {
    for (n[1] = fft.plan[3].start[1]; n[1] < end[1]; n[1]++) {
      for (n[2] = fft.plan[3].start[2]; n[2] < end[2]; n[2]++) {
        ind = (n[2] - fft.plan[3].start[2]) +
              fft.plan[3].new_mesh[2] *
                  ((n[1] - fft.plan[3].start[1]) +
                   (fft.plan[3].new_mesh[1] * (n[0] - fft.plan[3].start[0])));

        if ((n[KX] % (p3m.params.mesh[RX] / 2) == 0) &&
            (n[KY] % (p3m.params.mesh[RY] / 2) == 0) &&
            (n[KZ] % (p3m.params.mesh[RZ] / 2) == 0)) {
          p3m.g_force[ind] = 0.0;
        } else {
          const double denominator =
              perform_aliasing_sums_force<cao>(n, nominator);

          fak1 = p3m.d_op[RX][n[KX]] * nominator[RX] / box_l[RX] +
                 p3m.d_op[RY][n[KY]] * nominator[RY] / box_l[RY] +
                 p3m.d_op[RZ][n[KZ]] * nominator[RZ] / box_l[RZ];
          fak2 = Utils::sqr(p3m.d_op[RX][n[KX]] / box_l[RX]) +
                 Utils::sqr(p3m.d_op[RY][n[KY]] / box_l[RY]) +
                 Utils::sqr(p3m.d_op[RZ][n[KZ]] / box_l[RZ]);

          fak3 = fak1 / (fak2 * Utils::sqr(denominator));
          p3m.g_force[ind] = 2 * fak3 / (PI);
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

template <int cao> inline double perform_aliasing_sums_energy(int n[3]) {
  using Utils::int_pow;
  double numerator = 0.0, denominator = 0.0;
  /* lots of temporary variables... */
  double sx, sy, sz, f1, f2, mx, my, mz, nmx, nmy, nmz, nm2, expo;
  double limit = 30;

  f1 = Utils::sqr(PI / (p3m.params.alpha));

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
        nm2 =
            Utils::sqr(nmx / box_l[RX]) + Utils::sqr(nmy / box_l[RY]) + Utils::sqr(nmz / box_l[RZ]);
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
    size *= fft.plan[3].new_mesh[i];
    end[i] = fft.plan[3].start[i] + fft.plan[3].new_mesh[i];
    start[i] = fft.plan[3].start[i];
  }

  p3m.g_energy = Utils::realloc(p3m.g_energy, size * sizeof(double));

  /* Skip influence function calculation in tuning mode,
     the results need not be correct for timing. */
  if(p3m.params.tuning)
    return;

  ind = 0;

  for (n[0] = start[0]; n[0] < end[0]; n[0]++) {
    for (n[1] = start[1]; n[1] < end[1]; n[1]++) {
      for (n[2] = start[2]; n[2] < end[2]; n[2]++) {
        ind = (n[2] - start[2]) + fft.plan[3].new_mesh[2] * (n[1] - start[1]) +
              fft.plan[3].new_mesh[2] * fft.plan[3].new_mesh[1] *
                  (n[0] - start[0]);
        if ((n[KX] % (p3m.params.mesh[RX] / 2) == 0) &&
            (n[KY] % (p3m.params.mesh[RY] / 2) == 0) &&
            (n[KZ] % (p3m.params.mesh[RZ] / 2) == 0)) {
          p3m.g_energy[ind] = 0.0;
        }

        else
          p3m.g_energy[ind] = perform_aliasing_sums_energy<cao>(n) / PI;
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

/************************************************
 * Functions for P3M Parameter tuning
 * This tuning is based on P3M_tune by M. Deserno
 ************************************************/

#define P3M_TUNE_MAX_CUTS 50

/** get the minimal error for this combination of parameters. In fact,
    the real space error is tuned such that it contributes half of the
    total error, and then the Fourier space error is
    calculated. Returns the error and the optimal alpha, or 0 if this
    combination does not work at all */
static double p3m_get_accuracy(int mesh[3], int cao, double r_cut_iL,
                               double *_alpha_L, double *_rs_err,
                               double *_ks_err) {
  double rs_err, ks_err;
  double alpha_L;
  P3M_TRACE(fprintf(stderr,
                    "p3m_get_accuracy: mesh (%d, %d, %d), cao %d, r_cut %f ",
                    mesh[0], mesh[1], mesh[2], cao, r_cut_iL));

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
  /* calculate real space and k space error for this alpha_L */
  rs_err = p3m_real_space_error(coulomb.prefactor, r_cut_iL, p3m.sum_qpart,
                                p3m.sum_q2, alpha_L);
#ifdef CUDA
  if (coulomb.method == COULOMB_P3M_GPU)
    ks_err = p3m_k_space_error_gpu(coulomb.prefactor, mesh, cao, p3m.sum_qpart,
                                   p3m.sum_q2, alpha_L, box_l);
  else
#endif
    ks_err = p3m_k_space_error(coulomb.prefactor, mesh, cao, p3m.sum_qpart,
                               p3m.sum_q2, alpha_L);

  *_rs_err = rs_err;
  *_ks_err = ks_err;
  P3M_TRACE(fprintf(
      stderr, "resulting: alpha_L %g -> rs_err: %g, ks_err %g, total_err %g\n",
      alpha_L, rs_err, ks_err, sqrt(Utils::sqr(rs_err) + Utils::sqr(ks_err))));
  return sqrt(Utils::sqr(rs_err) + Utils::sqr(ks_err));
}

/** get the optimal alpha and the corresponding computation time for fixed
 * mesh,
 * cao, r_cut and alpha */
static double p3m_mcr_time(int mesh[3], int cao, double r_cut_iL,
                           double alpha_L) {
  /* rounded up 5000/n_charges timing force evaluations */
  int int_num = (5000 + p3m.sum_qpart) / p3m.sum_qpart;
  double int_time;

  /* broadcast p3m parameters for test run */
  if (coulomb.method != COULOMB_P3M && coulomb.method != COULOMB_ELC_P3M &&
      coulomb.method != COULOMB_P3M_GPU)
    coulomb.method = COULOMB_P3M;

  p3m.params.r_cut = r_cut_iL * box_l[0];
  p3m.params.r_cut_iL = r_cut_iL;
  p3m.params.mesh[0] = mesh[0];
  p3m.params.mesh[1] = mesh[1];
  p3m.params.mesh[2] = mesh[2];
  p3m.params.cao = cao;
  p3m.params.alpha_L = alpha_L;
  p3m.params.alpha = p3m.params.alpha_L * box_l_i[0];

  /* initialize p3m structures */
  mpi_bcast_coulomb_params();
  /* perform force calculation test */
  int_time = time_force_calc(int_num);
  P3M_TRACE(fprintf(stderr, "%d: test integration with mesh (%d %d %d), "
                            "r_cut_iL %lf, cao %d, alpha_L %lf returned %lf.\n",
                    this_node, mesh[0], mesh[1], mesh[2], r_cut_iL, cao,
                    alpha_L, int_time));
  return int_time;
}

/** get the optimal alpha and the corresponding computation time for fixed
   mesh,
   cao. The r_cut is determined via
    a simple bisection. Returns -1 if the force evaluation does not work, -2
   if
   there is no valid r_cut, and -3 if
    the charge assigment order is to large for this grid */
static double p3m_mc_time(char **log, int mesh[3], int cao, double r_cut_iL_min,
                          double r_cut_iL_max, double *_r_cut_iL,
                          double *_alpha_L, double *_accuracy) {
  double int_time;
  double r_cut_iL;
  double rs_err, ks_err;
  int i, n_cells;
  char b[5 * ES_DOUBLE_SPACE + 3 * ES_INTEGER_SPACE + 128];

  /* initial checks. */
  auto const k_cut = std::max(box_l[0] * cao / (2.0 * mesh[0]),
                              std::max(box_l[1] * cao / (2.0 * mesh[1]),
                                       box_l[2] * cao / (2.0 * mesh[2])));

  P3M_TRACE(fprintf(
      stderr, "p3m_mc_time: mesh=(%d, %d, %d), cao=%d, rmin=%f, rmax=%f\n",
      mesh[0], mesh[1], mesh[2], cao, r_cut_iL_min, r_cut_iL_max));
  if (cao >= std::min(mesh[0], std::min(mesh[1], mesh[2])) ||
      k_cut >= (std::min(min_box_l, min_local_box_l) - skin)) {
    sprintf(b, "%-4d %-3d cao too large for this mesh\n", mesh[0], cao);
    *log = strcat_alloc(*log, b);
    return -3;
  }

  /* Either low and high boundary are equal (for fixed cut), or the low border
     is initially 0 and therefore
     has infinite error estimate, as required. Therefore if the high boundary
     fails, there is no possible r_cut */
  if ((*_accuracy = p3m_get_accuracy(mesh, cao, r_cut_iL_max, _alpha_L, &rs_err,
                                     &ks_err)) > p3m.params.accuracy) {
    /* print result */
    P3M_TRACE(puts("p3m_mc_time: accuracy not achieved."));
    sprintf(b, "%-4d %-3d %.5e %.5e %.5e %.3e %.3e accuracy not achieved\n",
            mesh[0], cao, r_cut_iL_max, *_alpha_L, *_accuracy, rs_err, ks_err);
    *log = strcat_alloc(*log, b);
    return -2;
  }

  for (;;) {
    P3M_TRACE(fprintf(stderr, "p3m_mc_time: interval [%f,%f]\n", r_cut_iL_min,
                      r_cut_iL_max));
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
      elc_params.gap_size <= 1.1 * r_cut_iL * box_l[0]) {
    P3M_TRACE(fprintf(stderr, "p3m_mc_time: mesh (%d, %d, %d) cao %d r_cut %f "
                              "reject r_cut %f > gap %f\n",
                      mesh[0], mesh[1], mesh[2], cao, r_cut_iL,
                      2 * r_cut_iL * box_l[0], elc_params.gap_size));
    /* print result */
    sprintf(b, "%-4d %-3d %.5e %.5e %.5e %.3e %.3e conflict with ELC\n",
            mesh[0], cao, r_cut_iL, *_alpha_L, *_accuracy, rs_err, ks_err);
    *log = strcat_alloc(*log, b);
    return -P3M_TUNE_ELCTEST;
  }

  /* check whether this radius is too large, so that we would use less cells
   * than allowed */
  n_cells = 1;
  for (i = 0; i < 3; i++)
    n_cells *= (int)(floor(local_box_l[i] / (r_cut_iL * box_l[0] + skin)));
  if (n_cells < min_num_cells) {
    P3M_TRACE(fprintf(
        stderr,
        "p3m_mc_time: mesh (%d, %d, %d) cao %d r_cut %f reject n_cells %d\n",
        mesh[0], mesh[1], mesh[2], cao, r_cut_iL, n_cells));
    /* print result */
    sprintf(b, "%-4d %-3d %.5e %.5e %.5e %.3e %.3e radius dangerously high\n\n",
            mesh[0], cao, r_cut_iL, *_alpha_L, *_accuracy, rs_err, ks_err);
    *log = strcat_alloc(*log, b);
  }

  int_time = p3m_mcr_time(mesh, cao, r_cut_iL, *_alpha_L);
  if (int_time == -1) {
    *log = strcat_alloc(*log, "tuning failed, test integration not possible\n");
    return -P3M_TUNE_FAIL;
  }

  *_accuracy =
      p3m_get_accuracy(mesh, cao, r_cut_iL, _alpha_L, &rs_err, &ks_err);

  P3M_TRACE(fprintf(stderr,
                    "p3m_mc_time: mesh (%d, %d, %d) cao %d r_cut %f time %f\n",
                    mesh[0], mesh[1], mesh[2], cao, r_cut_iL, int_time));
  /* print result */
  sprintf(b, "%-4d %-3d %.5e %.5e %.5e %.3e %.3e %-8.2f\n", mesh[0], cao,
          r_cut_iL, *_alpha_L, *_accuracy, rs_err, ks_err, int_time);
  *log = strcat_alloc(*log, b);
  return int_time;
}

/** get the optimal alpha and the corresponding computation time for fixed
   mesh.
   *cao
    should contain an initial guess, which is then adapted by stepping up and
   down. Returns the time
    upon completion, -1 if the force evaluation does not work, and -2 if the
   accuracy cannot be met */
static double p3m_m_time(char **log, int mesh[3], int cao_min, int cao_max,
                         int *_cao, double r_cut_iL_min, double r_cut_iL_max,
                         double *_r_cut_iL, double *_alpha_L,
                         double *_accuracy) {
  double best_time = -1, tmp_time, tmp_r_cut_iL = 0.0, tmp_alpha_L = 0.0,
         tmp_accuracy = 0.0;
  /* in which direction improvement is possible. Initially, we dont know it
   * yet.
   */
  int final_dir = 0;
  int cao = *_cao;

  P3M_TRACE(fprintf(stderr, "p3m_m_time: mesh=(%d, %d %d), cao_min=%d, "
                            "cao_max=%d, rmin=%f, rmax=%f\n",
                    mesh[0], mesh[1], mesh[2], cao_min, cao_max, r_cut_iL_min,
                    r_cut_iL_max));
  /* the initial step sets a timing mark. If there is no valid r_cut, we can
     only try
     to increase cao to increase the obtainable precision of the far formula.
     */
  do {
    tmp_time = p3m_mc_time(log, mesh, cao, r_cut_iL_min, r_cut_iL_max,
                           &tmp_r_cut_iL, &tmp_alpha_L, &tmp_accuracy);
    /* bail out if the force evaluation is not working */
    if (tmp_time == -1)
      return -1;
    /* cao is too large for this grid, but still the accuracy cannot be
     * achieved, give up */
    if (tmp_time == -3) {
      P3M_TRACE(fprintf(stderr, "p3m_m_time: no possible cao found\n"));
      return -2;
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
    /* the required accuracy could not be obtained, try higher caos. Therefore
       optimisation can only be
       obtained with even higher caos, but not lower ones */
    P3M_TRACE(fprintf(stderr, "p3m_m_time: doesn't give precision, step up\n"));
    cao++;
    final_dir = 1;
  } while (cao <= cao_max);
  /* with this mesh, the required accuracy cannot be obtained. */
  if (cao > cao_max)
    return -2;

  /* at the boundaries, only the opposite direction can be used for
   * optimisation
   */
  if (cao == cao_min)
    final_dir = 1;
  else if (cao == cao_max)
    final_dir = -1;

  P3M_TRACE(
      fprintf(stderr, "p3m_m_time: final constraints dir %d\n", final_dir));

  if (final_dir == 0) {
    /* check in which direction we can optimise. Both directions are possible
     */
    double dir_times[3];
    for (final_dir = -1; final_dir <= 1; final_dir += 2) {
      dir_times[final_dir + 1] = tmp_time =
          p3m_mc_time(log, mesh, cao + final_dir, r_cut_iL_min, r_cut_iL_max,
                      &tmp_r_cut_iL, &tmp_alpha_L, &tmp_accuracy);
      /* bail out on errors, as usual */
      if (tmp_time == -1)
        return -P3M_TUNE_FAIL;
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
        P3M_TRACE(fprintf(
            stderr, "p3m_m_time: mesh=(%d, %d, %d) final cao=%d time=%f\n",
            mesh[0], mesh[1], mesh[2], cao, best_time));
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

  P3M_TRACE(
      fprintf(stderr, "p3m_m_time: optimise in direction %d\n", final_dir));

  /* move cao into the optimisation direction until we do not gain anymore. */
  for (; cao >= cao_min && cao <= cao_max; cao += final_dir) {
    tmp_time = p3m_mc_time(log, mesh, cao, r_cut_iL_min, r_cut_iL_max,
                           &tmp_r_cut_iL, &tmp_alpha_L, &tmp_accuracy);
    /* bail out on errors, as usual */
    if (tmp_time == -1)
      return -1;
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
  P3M_TRACE(fprintf(
      stderr, "p3m_m_time: mesh=(%d, %d, %d) final cao=%d r_cut=%f time=%f\n",
      mesh[0], mesh[1], mesh[2], *_cao, *_r_cut_iL, best_time));
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
  int tune_mesh = 0; // boolean to indicate if mesh should be tuned

  if (p3m.params.epsilon != P3M_EPSILON_METALLIC) {
    if (!((box_l[0] == box_l[1]) && (box_l[1] == box_l[2]))) {
      *log = strcat_alloc(
          *log, "{049 P3M_init: Nonmetallic epsilon requires cubic box} ");
      return ES_ERROR;
    }
  }

  if (p3m_sanity_checks_system()) {
    return ES_ERROR;
  }

  /* preparation */
  mpi_bcast_event(P3M_COUNT_CHARGES);
  /* Print Status */
  sprintf(b,
          "P3M tune parameters: Accuracy goal = %.5e prefactor = %.5e \n",
          p3m.params.accuracy, coulomb.prefactor);
  *log = strcat_alloc(*log, b);
  sprintf(b, "System: box_l = %.5e # charged part = %d Sum[q_i^2] = %.5e\n",
          box_l[0], p3m.sum_qpart, p3m.sum_q2);
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
        pow(p3m.sum_qpart / (box_l[0] * box_l[1] * box_l[2]), 1.0 / 3.0);
    mesh_density_max = 512 / pow(box_l[0] * box_l[1] * box_l[2], 1.0 / 3.0);
    tune_mesh = 1;
    /* this limits the tried meshes if the accuracy cannot
       be obtained with smaller meshes, but normally not all these
       meshes have to be tested */
    /* avoid using more than 1 GB of FFT arrays (per default, see config.hpp)
     */

    P3M_TRACE(fprintf(stderr,
                      "%d: starting with meshdensity %lf, using at most %lf.\n",
                      this_node, mesh_density_min, mesh_density_max));

  } else if (p3m.params.mesh[1] == -1 && p3m.params.mesh[2] == -1) {
    mesh_density = mesh_density_min = mesh_density_max =
        p3m.params.mesh[0] / box_l[0];
    p3m.params.mesh[1] = lround(mesh_density * box_l[1]);
    p3m.params.mesh[2] = lround(mesh_density * box_l[2]);
    if (p3m.params.mesh[1] % 2 == 1)
      p3m.params.mesh[1]++; // Make sure that the mesh is even in all directions
    if (p3m.params.mesh[2] % 2 == 1)
      p3m.params.mesh[2]++;

    sprintf(b, "fixed mesh %d %d %d\n", p3m.params.mesh[0], p3m.params.mesh[1],
            p3m.params.mesh[2]);
    *log = strcat_alloc(*log, b);
  } else {
    mesh_density = mesh_density_min = mesh_density_max =
        p3m.params.mesh[0] / box_l[0];

    sprintf(b, "fixed mesh %d %d %d\n", p3m.params.mesh[0], p3m.params.mesh[1],
            p3m.params.mesh[2]);
    *log = strcat_alloc(*log, b);
  }

  if (p3m.params.r_cut_iL == 0.0) {
    r_cut_iL_min = 0;
    r_cut_iL_max = std::min(min_local_box_l, min_box_l / 2.0) - skin;
    r_cut_iL_min *= box_l_i[0];
    r_cut_iL_max *= box_l_i[0];
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

    P3M_TRACE(fprintf(stderr, "%d: trying meshdensity %lf.\n", this_node,
                      mesh_density));

    if (tune_mesh) {
      tmp_mesh[0] = lround(box_l[0] * mesh_density);
      tmp_mesh[1] = lround(box_l[1] * mesh_density);
      tmp_mesh[2] = lround(box_l[2] * mesh_density);
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

    tmp_time =
        p3m_m_time(log, tmp_mesh, cao_min, cao_max, &tmp_cao, r_cut_iL_min,
                   r_cut_iL_max, &tmp_r_cut_iL, &tmp_alpha_L, &tmp_accuracy);
    /* some error occured during the tuning force evaluation */
    P3M_TRACE(fprintf(stderr, "delta_accuracy: %lf tune time: %lf\n",
                      p3m.params.accuracy - tmp_accuracy, tmp_time));
    //    if (tmp_time == -1) con;
    /* this mesh does not work at all */
    if (tmp_time < 0.0)
      continue;

    /* the optimum r_cut for this mesh is the upper limit for higher meshes,
       everything else is slower */
    if (coulomb.method == COULOMB_P3M)
      r_cut_iL_max = tmp_r_cut_iL;

    /* new optimum */
    if (tmp_time < time_best) {
      P3M_TRACE(fprintf(stderr,
                        "Found new optimum: time %lf, mesh (%d %d %d)\n",
                        tmp_time, tmp_mesh[0], tmp_mesh[1], tmp_mesh[2]));
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
      P3M_TRACE(fprintf(stderr,
                        "%d: %lf is mush slower then best time, aborting.\n",
                        this_node, tmp_time));
      break;
    }
  }

  P3M_TRACE(fprintf(stderr, "%d: finished tuning, best time: %lf\n", this_node,
                    time_best));
  if (time_best == 1e20) {
    *log = strcat_alloc(*log,
                        "failed to tune P3M parameters to required accuracy\n");
    return ES_ERROR;
  }

  /* set tuned p3m parameters */
  p3m.params.tuning = false;
  p3m.params.r_cut = r_cut_iL * box_l[0];
  p3m.params.r_cut_iL = r_cut_iL;
  p3m.params.mesh[0] = mesh[0];
  p3m.params.mesh[1] = mesh[1];
  p3m.params.mesh[2] = mesh[2];
  p3m.params.cao = cao;
  p3m.params.alpha_L = alpha_L;
  p3m.params.alpha = p3m.params.alpha_L * box_l_i[0];
  p3m.params.accuracy = accuracy;
  /* broadcast tuned p3m parameters */
  P3M_TRACE(fprintf(stderr, "%d: Broadcasting P3M parameters: mesh: (%d %d "
                            "%d), cao: %d, alpha_L: %lf, acccuracy: %lf\n",
                    this_node, p3m.params.mesh[0], p3m.params.mesh[1],
                    p3m.params.mesh[2], p3m.params.cao, p3m.params.alpha_L,
                    p3m.params.accuracy));
  mpi_bcast_coulomb_params();

  P3M_TRACE(p3m_print());

  /* Tell the user about the outcome */
  sprintf(
      b, "\nresulting parameters:\n%-4d %-4d %-4d %-3d %.5e %.5e %.5e %-8.2f\n",
      mesh[0], mesh[1], mesh[2], cao, r_cut_iL, alpha_L, accuracy, time_best);
  *log = strcat_alloc(*log, b);
  return ES_OK;
}

void p3m_count_charged_particles() {
  double node_sums[3], tot_sums[3];

  P3M_TRACE(fprintf(stderr, "%d: p3m_count_charged_particles\n", this_node));

  for (int i = 0; i < 3; i++) {
    node_sums[i] = 0.0;
    tot_sums[i] = 0.0;
  }

  for (auto const &p : local_cells.particles()) {
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

  P3M_TRACE(fprintf(
      stderr, "%d: p3m.sum_qpart: %d, p3m.sum_q2: %lf, total_charge %lf\n",
      this_node, p3m.sum_qpart, p3m.sum_q2, sqrt(p3m.square_sum_q)));
}

double p3m_real_space_error(double prefac, double r_cut_iL, int n_c_part,
                            double sum_q2, double alpha_L) {
  return (2.0 * prefac * sum_q2 * exp(-Utils::sqr(r_cut_iL * alpha_L))) /
         sqrt((double)n_c_part * r_cut_iL * box_l[0] * box_l[0] * box_l[1] *
              box_l[2]);
}

double p3m_k_space_error(double prefac, int mesh[3], int cao, int n_c_part,
                         double sum_q2, double alpha_L) {
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
         (box_l[1] * box_l[2]);
}

void p3m_tune_aliasing_sums(int nx, int ny, int nz, int mesh[3],
                            double mesh_i[3], int cao, double alpha_L_i,
                            double *alias1, double *alias2) {

  int mx, my, mz;
  double nmx, nmy, nmz;
  double fnmx, fnmy, fnmz;

  double ex, ex2, nm2, U2, factor1;

  factor1 = Utils::sqr(PI * alpha_L_i);

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

/************************************************************/

void p3m_calc_local_ca_mesh() {
  int i;
  int ind[3];
  /* total skin size */
  double full_skin[3];

  for (i = 0; i < 3; i++)
    full_skin[i] = p3m.params.cao_cut[i] + skin + p3m.params.additional_mesh[i];

  /* inner left down grid point (global index) */
  for (i = 0; i < 3; i++)
    p3m.local_mesh.in_ld[i] =
        (int)ceil(my_left[i] * p3m.params.ai[i] - p3m.params.mesh_off[i]);
  /* inner up right grid point (global index) */
  for (i = 0; i < 3; i++)
    p3m.local_mesh.in_ur[i] =
        (int)floor(my_right[i] * p3m.params.ai[i] - p3m.params.mesh_off[i]);

  /* correct roundof errors at boundary */
  for (i = 0; i < 3; i++) {
    if ((my_right[i] * p3m.params.ai[i] - p3m.params.mesh_off[i]) -
            p3m.local_mesh.in_ur[i] <
        ROUND_ERROR_PREC)
      p3m.local_mesh.in_ur[i]--;
    if (1.0 + (my_left[i] * p3m.params.ai[i] - p3m.params.mesh_off[i]) -
            p3m.local_mesh.in_ld[i] <
        ROUND_ERROR_PREC)
      p3m.local_mesh.in_ld[i]--;
  }
  /* inner grid dimensions */
  for (i = 0; i < 3; i++)
    p3m.local_mesh.inner[i] =
        p3m.local_mesh.in_ur[i] - p3m.local_mesh.in_ld[i] + 1;
  /* index of left down grid point in global mesh */
  for (i = 0; i < 3; i++)
    p3m.local_mesh.ld_ind[i] =
        (int)ceil((my_left[i] - full_skin[i]) * p3m.params.ai[i] -
                  p3m.params.mesh_off[i]);
  /* left down margin */
  for (i = 0; i < 3; i++)
    p3m.local_mesh.margin[i * 2] =
        p3m.local_mesh.in_ld[i] - p3m.local_mesh.ld_ind[i];
  /* up right grid point */
  for (i = 0; i < 3; i++)
    ind[i] = (int)floor((my_right[i] + full_skin[i]) * p3m.params.ai[i] -
                        p3m.params.mesh_off[i]);
  /* correct roundof errors at up right boundary */
  for (i = 0; i < 3; i++)
    if (((my_right[i] + full_skin[i]) * p3m.params.ai[i] -
         p3m.params.mesh_off[i]) -
            ind[i] ==
        0)
      ind[i]--;
  /* up right margin */
  for (i = 0; i < 3; i++)
    p3m.local_mesh.margin[(i * 2) + 1] = ind[i] - p3m.local_mesh.in_ur[i];

  /* grid dimension */
  p3m.local_mesh.size = 1;
  for (i = 0; i < 3; i++) {
    p3m.local_mesh.dim[i] = ind[i] - p3m.local_mesh.ld_ind[i] + 1;
    p3m.local_mesh.size *= p3m.local_mesh.dim[i];
  }
  /* reduce inner grid indices from global to local */
  for (i = 0; i < 3; i++)
    p3m.local_mesh.in_ld[i] = p3m.local_mesh.margin[i * 2];
  for (i = 0; i < 3; i++)
    p3m.local_mesh.in_ur[i] =
        p3m.local_mesh.margin[i * 2] + p3m.local_mesh.inner[i];

  p3m.local_mesh.q_2_off = p3m.local_mesh.dim[2] - p3m.params.cao;
  p3m.local_mesh.q_21_off =
      p3m.local_mesh.dim[2] * (p3m.local_mesh.dim[1] - p3m.params.cao);
}

void p3m_calc_lm_ld_pos() {
  int i;
  /* spacial position of left down mesh point */
  for (i = 0; i < 3; i++) {
    p3m.local_mesh.ld_pos[i] =
        (p3m.local_mesh.ld_ind[i] + p3m.params.mesh_off[i]) * p3m.params.a[i];
  }
}

void p3m_init_a_ai_cao_cut() {
  int i;
  for (i = 0; i < 3; i++) {
    p3m.params.ai[i] = (double)p3m.params.mesh[i] / box_l[i];
    p3m.params.a[i] = 1.0 / p3m.params.ai[i];
    p3m.params.cao_cut[i] = 0.5 * p3m.params.a[i] * p3m.params.cao;
  }
}

int p3m_sanity_checks_boxl() {
  // char *errtxt;
  int i, ret = 0;
  for (i = 0; i < 3; i++) {
    /* check k-space cutoff */
    if (p3m.params.cao_cut[i] >= 0.5 * box_l[i]) {
      runtimeErrorMsg() << "P3M_init: k-space cutoff " << p3m.params.cao_cut[i]
                        << " is larger than half of box dimension " << box_l[i];
      ret = 1;
    }
    if (p3m.params.cao_cut[i] >= local_box_l[i]) {
      runtimeErrorMsg() << "P3M_init: k-space cutoff " << p3m.params.cao_cut[i]
                        << " is larger than local box dimension "
                        << local_box_l[i];
      ret = 1;
    }
  }

  return ret;
}

/**
 * @brief General sanity checks independent of p3m parameters.
 *
 * @return 0 if ok, 1 on error.
 */
int p3m_sanity_checks_system() {
  int ret = 0;

  if (!PERIODIC(0) || !PERIODIC(1) || !PERIODIC(2)) {
    runtimeErrorMsg() << "P3M requires periodicity 1 1 1";
    ret = 1;
  }

  if (cell_structure.type != CELL_STRUCTURE_DOMDEC) {
    runtimeErrorMsg()
        << "P3M at present requires the domain decomposition cell system";
    ret = 1;
  }

  if (node_grid[0] < node_grid[1] || node_grid[1] < node_grid[2]) {
    runtimeErrorMsg() << "P3M_init: node grid must be sorted, largest first";
    ret = 1;
  }

  if (p3m.params.epsilon != P3M_EPSILON_METALLIC) {
    if (!((p3m.params.mesh[0] == p3m.params.mesh[1]) &&
          (p3m.params.mesh[1] == p3m.params.mesh[2]))) {
      runtimeErrorMsg() << "P3M_init: Nonmetallic epsilon requires cubic box";
      ret = 1;
    }
  }

  return ret;
}

int p3m_sanity_checks() {
  int ret = 0;

  if (p3m_sanity_checks_system())
    ret = 1;

  if (p3m_sanity_checks_boxl())
    ret = 1;

  if (p3m.params.mesh[0] == 0) {
    runtimeErrorMsg() << "P3M_init: mesh size is not yet set";
    ret = 1;
  }
  if (p3m.params.cao == 0) {
    runtimeErrorMsg() << "P3M_init: cao is not yet set";
    ret = 1;
  }
  if (p3m.params.alpha < 0.0) {
    runtimeErrorMsg() << "P3M_init: alpha must be >0";
    ret = 1;
  }

  return ret;
}

void p3m_calc_send_mesh() {
  int i, j, evenodd;
  int done[3] = {0, 0, 0};
  MPI_Status status;
  /* send grids */
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      /* left */
      p3m.sm.s_ld[i * 2][j] = 0 + done[j] * p3m.local_mesh.margin[j * 2];
      if (j == i)
        p3m.sm.s_ur[i * 2][j] = p3m.local_mesh.margin[j * 2];
      else
        p3m.sm.s_ur[i * 2][j] = p3m.local_mesh.dim[j] -
                                done[j] * p3m.local_mesh.margin[(j * 2) + 1];
      /* right */
      if (j == i)
        p3m.sm.s_ld[(i * 2) + 1][j] = p3m.local_mesh.in_ur[j];
      else
        p3m.sm.s_ld[(i * 2) + 1][j] =
            0 + done[j] * p3m.local_mesh.margin[j * 2];
      p3m.sm.s_ur[(i * 2) + 1][j] =
          p3m.local_mesh.dim[j] - done[j] * p3m.local_mesh.margin[(j * 2) + 1];
    }
    done[i] = 1;
  }
  p3m.sm.max = 0;
  for (i = 0; i < 6; i++) {
    p3m.sm.s_size[i] = 1;
    for (j = 0; j < 3; j++) {
      p3m.sm.s_dim[i][j] = p3m.sm.s_ur[i][j] - p3m.sm.s_ld[i][j];
      p3m.sm.s_size[i] *= p3m.sm.s_dim[i][j];
    }
    if (p3m.sm.s_size[i] > p3m.sm.max)
      p3m.sm.max = p3m.sm.s_size[i];
  }
  /* communication */
  for (i = 0; i < 6; i++) {
    if (i % 2 == 0)
      j = i + 1;
    else
      j = i - 1;
    if (node_neighbors[i] != this_node) {
      /* two step communication: first all even positions than all odd */
      for (evenodd = 0; evenodd < 2; evenodd++) {
        if ((node_pos[i / 2] + evenodd) % 2 == 0)
          MPI_Send(&(p3m.local_mesh.margin[i]), 1, MPI_INT, node_neighbors[i],
                   REQ_P3M_INIT, comm_cart);
        else
          MPI_Recv(&(p3m.local_mesh.r_margin[j]), 1, MPI_INT, node_neighbors[j],
                   REQ_P3M_INIT, comm_cart, &status);
      }
    } else {
      p3m.local_mesh.r_margin[j] = p3m.local_mesh.margin[i];
    }
  }
  /* recv grids */
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++) {
      if (j == i) {
        p3m.sm.r_ld[i * 2][j] =
            p3m.sm.s_ld[i * 2][j] + p3m.local_mesh.margin[2 * j];
        p3m.sm.r_ur[i * 2][j] =
            p3m.sm.s_ur[i * 2][j] + p3m.local_mesh.r_margin[2 * j];
        p3m.sm.r_ld[(i * 2) + 1][j] =
            p3m.sm.s_ld[(i * 2) + 1][j] - p3m.local_mesh.r_margin[(2 * j) + 1];
        p3m.sm.r_ur[(i * 2) + 1][j] =
            p3m.sm.s_ur[(i * 2) + 1][j] - p3m.local_mesh.margin[(2 * j) + 1];
      } else {
        p3m.sm.r_ld[i * 2][j] = p3m.sm.s_ld[i * 2][j];
        p3m.sm.r_ur[i * 2][j] = p3m.sm.s_ur[i * 2][j];
        p3m.sm.r_ld[(i * 2) + 1][j] = p3m.sm.s_ld[(i * 2) + 1][j];
        p3m.sm.r_ur[(i * 2) + 1][j] = p3m.sm.s_ur[(i * 2) + 1][j];
      }
    }
  for (i = 0; i < 6; i++) {
    p3m.sm.r_size[i] = 1;
    for (j = 0; j < 3; j++) {
      p3m.sm.r_dim[i][j] = p3m.sm.r_ur[i][j] - p3m.sm.r_ld[i][j];
      p3m.sm.r_size[i] *= p3m.sm.r_dim[i][j];
    }
    if (p3m.sm.r_size[i] > p3m.sm.max)
      p3m.sm.max = p3m.sm.r_size[i];
  }
}

/************************************************/

void p3m_scaleby_box_l() {
  if (coulomb.prefactor < 0.0) {
    runtimeErrorMsg() << "The Coulomb prefactor has to be >=0";
    return;
  }

  p3m.params.r_cut = p3m.params.r_cut_iL * box_l[0];
  p3m.params.alpha = p3m.params.alpha_L * box_l_i[0];
  p3m_init_a_ai_cao_cut();
  p3m_calc_lm_ld_pos();
  p3m_sanity_checks_boxl();
  p3m_calc_influence_function_force();
  p3m_calc_influence_function_energy();
}

/************************************************/

void p3m_calc_kspace_stress(double *stress) {
  /**
  Calculates the long range electrostatics part of the stress tensor. This is part Pi_{dir, alpha,beta} in the paper by Essmann et al "A smooth particle mesh Ewald method", The Journal of Chemical Physics 103, 8577 (1995); doi: 10.1063/1.470117. The part Pi_{corr, alpha, beta} in the Essmann paper is not present here since M is the empty set in our simulations.
  */
  if (p3m.sum_q2 > 0) {
    double *node_k_space_stress;
    double *k_space_stress;
    double force_prefac, node_k_space_energy, sqk, vterm, kx, ky, kz;

    int j[3], i, ind = 0;
    // ordering after fourier transform
    node_k_space_stress = (double *)Utils::malloc(9 * sizeof(double));
    k_space_stress = (double *)Utils::malloc(9 * sizeof(double));

    for (i = 0; i < 9; i++) {
      node_k_space_stress[i] = 0.0;
      k_space_stress[i] = 0.0;
    }

    p3m_gather_fft_grid(p3m.rs_mesh);
    fft_perform_forw(p3m.rs_mesh);
    force_prefac = coulomb.prefactor / (2.0 * box_l[0] * box_l[1] * box_l[2]);

    for (j[0] = 0; j[0] < fft.plan[3].new_mesh[RX]; j[0]++) {
      for (j[1] = 0; j[1] < fft.plan[3].new_mesh[RY]; j[1]++) {
        for (j[2] = 0; j[2] < fft.plan[3].new_mesh[RZ]; j[2]++) {
          kx = 2.0 * PI * p3m.d_op[RX][j[KX] + fft.plan[3].start[KX]] /
               box_l[RX];
          ky = 2.0 * PI * p3m.d_op[RY][j[KY] + fft.plan[3].start[KY]] /
               box_l[RY];
          kz = 2.0 * PI * p3m.d_op[RZ][j[KZ] + fft.plan[3].start[KZ]] /
               box_l[RZ];
          sqk = Utils::sqr(kx) + Utils::sqr(ky) + Utils::sqr(kz);
          if (sqk == 0) {
            node_k_space_energy = 0.0;
            vterm = 0.0;
          } else {
            vterm = -2.0 * (1 / sqk + Utils::sqr(1.0 / 2.0 / p3m.params.alpha));
            node_k_space_energy =
                p3m.g_energy[ind] *
                (Utils::sqr(p3m.rs_mesh[2 * ind]) + Utils::sqr(p3m.rs_mesh[2 * ind + 1]));
          }
          ind++;
          node_k_space_stress[0] +=
              node_k_space_energy * (1.0 + vterm * Utils::sqr(kx)); /* sigma_xx */
          node_k_space_stress[1] +=
              node_k_space_energy * (vterm * kx * ky); /* sigma_xy */
          node_k_space_stress[2] +=
              node_k_space_energy * (vterm * kx * kz); /* sigma_xz */

          node_k_space_stress[3] +=
              node_k_space_energy * (vterm * kx * ky); /* sigma_yx */
          node_k_space_stress[4] +=
              node_k_space_energy * (1.0 + vterm * Utils::sqr(ky)); /* sigma_yy */
          node_k_space_stress[5] +=
              node_k_space_energy * (vterm * ky * kz); /* sigma_yz */

          node_k_space_stress[6] +=
              node_k_space_energy * (vterm * kx * kz); /* sigma_zx */
          node_k_space_stress[7] +=
              node_k_space_energy * (vterm * ky * kz); /* sigma_zy */
          node_k_space_stress[8] +=
              node_k_space_energy * (1.0 + vterm * Utils::sqr(kz)); /* sigma_zz */
        }
      }
    }

    MPI_Reduce(node_k_space_stress, k_space_stress, 9, MPI_DOUBLE, MPI_SUM, 0,
               comm_cart);
    if (this_node == 0) {
      for (i = 0; i < 9; i++) {
        stress[i] = k_space_stress[i] * force_prefac;
      }
    }
    free(node_k_space_stress);
    free(k_space_stress);
  }
}

/************************************************/

/*********************** miscelanea of functions
 * *************************************/

/************************************************
 * Debug functions printing p3m structures
 ************************************************/

void p3m_p3m_print_struct(p3m_parameter_struct ps) {
  fprintf(stderr, "%d: p3m_parameter_struct: \n", this_node);
  fprintf(stderr, "   alpha_L=%f, r_cut_iL=%f \n", ps.alpha_L, ps.r_cut_iL);
  fprintf(stderr, "   mesh=(%d,%d,%d), mesh_off=(%.4f,%.4f,%.4f)\n", ps.mesh[0],
          ps.mesh[1], ps.mesh[2], ps.mesh_off[0], ps.mesh_off[1],
          ps.mesh_off[2]);
  fprintf(stderr, "   cao=%d, inter=%d, epsilon=%f\n", ps.cao, ps.inter,
          ps.epsilon);
  fprintf(stderr, "   cao_cut=(%f,%f,%f)\n", ps.cao_cut[0], ps.cao_cut[1],
          ps.cao_cut[2]);
  fprintf(stderr, "   a=(%f,%f,%f), ai=(%f,%f,%f)\n", ps.a[0], ps.a[1], ps.a[2],
          ps.ai[0], ps.ai[1], ps.ai[2]);
}

#endif /* of P3M */
