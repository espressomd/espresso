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

/** \file iccp3m.cpp
  Detailed Information about the method is included in the corresponding header
  file \ref iccp3m.hpp.
*/

#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>

#include "elc.hpp"
#include "iccp3m.hpp"
#include "mmm1d.hpp"
#include "mmm2d.hpp"
#include "p3m.hpp"
#include "p3m_gpu.hpp"

#include "communication.hpp"

#include "cells.hpp"
#include "config.hpp"
#include "forces.hpp"
#include "global.hpp"
#include "interaction_data.hpp"
#include "particle_data.hpp"
#include "utils.hpp"
#include "initialize.hpp"

#include "short_range_loop.hpp"
#include "utils/NoOp.hpp"

using std::ostringstream;

#ifdef ELECTROSTATICS

iccp3m_struct iccp3m_cfg;

int iccp3m_initialized = 0;
/* functions that are used in icc* to compute the electric field acting on the
 * induced charges, excluding forces other than the electrostatic ones */
void init_forces_iccp3m();
void calc_long_range_forces_iccp3m();

inline void init_local_particle_force_iccp3m(Particle *part);
inline void init_ghost_force_iccp3m(Particle *part);

/** Calculation of the electrostatic forces between source charges (= real
 * charges) and wall charges. For each electrostatic method the proper functions
 * for short and long range parts are called. Long Range Parts are calculated
 * directly, short range parts need helper functions according to the particle
 * data organisation. A modified version of \ref force_calc in \ref forces.hpp.
 */
void force_calc_iccp3m();

/** Variant of add_non_bonded_pair_force where only coulomb
 *  contributions are calculated   */
inline void add_non_bonded_pair_force_iccp3m(Particle *p1, Particle *p2,
                                             double d[3], double dist,
                                             double dist2) {
  /* IA_parameters *ia_params = get_ia_param(p1->p.type,p2->p.type);*/
  double force[3] = {0, 0, 0};

  FORCE_TRACE(fprintf(stderr, "%d: interaction %d<->%d dist %f\n", this_node,
                      p1->p.identity, p2->p.identity, dist));

  /***********************************************/
  /* long range electrostatics                   */
  /***********************************************/

  /* real space coulomb */
  auto const q1q2 = p1->p.q * p2->p.q;
  switch (coulomb.method) {
#ifdef P3M
  case COULOMB_ELC_P3M:
    if (q1q2)
      p3m_add_pair_force(q1q2, d, dist2, dist, force);
    break;
  case COULOMB_P3M_GPU:
  case COULOMB_P3M:
    if (q1q2)
      p3m_add_pair_force(q1q2, d, dist2, dist, force);
    break;
#endif /* P3M */
  case COULOMB_MMM1D:
    if (q1q2)
      add_mmm1d_coulomb_pair_force(q1q2, d, dist2, dist, force);
    break;
  case COULOMB_MMM2D:
    if (q1q2)
      add_mmm2d_coulomb_pair_force(q1q2, d, dist2, dist, force);
    break;
  default:
    break;
  }

  /***********************************************/
  /* add total nonbonded forces to particle      */
  /***********************************************/
  for (int j = 0; j < 3; j++) {
    p1->f.f[j] += force[j];
    p2->f.f[j] -= force[j];
  }
  /***********************************************/
}

void iccp3m_set_initialized() { iccp3m_initialized = 1; }

void iccp3m_init(void) {
  iccp3m_cfg.set_flag = 0;
  iccp3m_cfg.areas = nullptr;
  iccp3m_cfg.ein = nullptr;
  iccp3m_cfg.nvectorx = nullptr;
  iccp3m_cfg.nvectory = nullptr;
  iccp3m_cfg.nvectorz = nullptr;
  iccp3m_cfg.extx = 0;
  iccp3m_cfg.exty = 0;
  iccp3m_cfg.extz = 0;
  iccp3m_cfg.first_id = 0;
  iccp3m_cfg.num_iteration = 30;
  iccp3m_cfg.convergence = 1e-2;
  iccp3m_cfg.relax = 0.7;
  iccp3m_cfg.eout = 1;
  iccp3m_cfg.citeration = 0;
}

void iccp3m_alloc_lists() {
  iccp3m_cfg.areas = (double *)Utils::realloc(
      iccp3m_cfg.areas, (iccp3m_cfg.n_ic) * sizeof(double));
  iccp3m_cfg.ein = (double *)Utils::realloc(iccp3m_cfg.ein,
                                            (iccp3m_cfg.n_ic) * sizeof(double));
  iccp3m_cfg.nvectorx = (double *)Utils::realloc(
      iccp3m_cfg.nvectorx, (iccp3m_cfg.n_ic) * sizeof(double));
  iccp3m_cfg.nvectory = (double *)Utils::realloc(
      iccp3m_cfg.nvectory, (iccp3m_cfg.n_ic) * sizeof(double));
  iccp3m_cfg.nvectorz = (double *)Utils::realloc(
      iccp3m_cfg.nvectorz, (iccp3m_cfg.n_ic) * sizeof(double));
  iccp3m_cfg.sigma = (double *)Utils::realloc(
      iccp3m_cfg.sigma, (iccp3m_cfg.n_ic) * sizeof(double));
}

int bcast_iccp3m_cfg(void) {
  int i;
  MPI_Bcast((int *)&iccp3m_cfg.n_ic, 1, MPI_INT, 0, comm_cart);
  /* allocates Memory on slave nodes
   * */
  if (this_node != 0) {
    iccp3m_cfg.areas = (double *)Utils::realloc(
        iccp3m_cfg.areas, (iccp3m_cfg.n_ic) * sizeof(double));
    iccp3m_cfg.ein = (double *)Utils::realloc(
        iccp3m_cfg.ein, (iccp3m_cfg.n_ic) * sizeof(double));
    iccp3m_cfg.nvectorx = (double *)Utils::realloc(
        iccp3m_cfg.nvectorx, (iccp3m_cfg.n_ic) * sizeof(double));
    iccp3m_cfg.nvectory = (double *)Utils::realloc(
        iccp3m_cfg.nvectory, (iccp3m_cfg.n_ic) * sizeof(double));
    iccp3m_cfg.nvectorz = (double *)Utils::realloc(
        iccp3m_cfg.nvectorz, (iccp3m_cfg.n_ic) * sizeof(double));
    iccp3m_cfg.sigma = (double *)Utils::realloc(
        iccp3m_cfg.sigma, (iccp3m_cfg.n_ic) * sizeof(double));
  }

  MPI_Bcast((int *)&iccp3m_cfg.num_iteration, 1, MPI_INT, 0, comm_cart);
  MPI_Bcast((int *)&iccp3m_cfg.first_id, 1, MPI_INT, 0, comm_cart);
  MPI_Bcast((double *)&iccp3m_cfg.convergence, 1, MPI_DOUBLE, 0, comm_cart);
  MPI_Bcast((double *)&iccp3m_cfg.eout, 1, MPI_DOUBLE, 0, comm_cart);
  MPI_Bcast((double *)&iccp3m_cfg.relax, 1, MPI_DOUBLE, 0, comm_cart);

  /* broadcast the vectors element by element. This is slow
   * but safe and only performed at the beginning of each simulation*/
  for (i = 0; i < iccp3m_cfg.n_ic; i++) {
    MPI_Bcast((double *)&iccp3m_cfg.areas[i], 1, MPI_DOUBLE, 0, comm_cart);
    MPI_Bcast((double *)&iccp3m_cfg.ein[i], 1, MPI_DOUBLE, 0, comm_cart);
    MPI_Bcast((double *)&iccp3m_cfg.nvectorx[i], 1, MPI_DOUBLE, 0, comm_cart);
    MPI_Bcast((double *)&iccp3m_cfg.nvectory[i], 1, MPI_DOUBLE, 0, comm_cart);
    MPI_Bcast((double *)&iccp3m_cfg.nvectorz[i], 1, MPI_DOUBLE, 0, comm_cart);
    MPI_Bcast((double *)&iccp3m_cfg.sigma[i], 1, MPI_DOUBLE, 0, comm_cart);
  }
  MPI_Bcast((double *)&iccp3m_cfg.extx, 1, MPI_DOUBLE, 0, comm_cart);
  MPI_Bcast((double *)&iccp3m_cfg.exty, 1, MPI_DOUBLE, 0, comm_cart);
  MPI_Bcast((double *)&iccp3m_cfg.extz, 1, MPI_DOUBLE, 0, comm_cart);

  MPI_Bcast(&iccp3m_cfg.citeration, 1, MPI_INT, 0, comm_cart);
  MPI_Bcast(&iccp3m_cfg.set_flag, 1, MPI_INT, 0, comm_cart);

  return 0;
}

int iccp3m_iteration() {
  double fdot, hold, hnew, del_eps, diff = 0.0, difftemp = 0.0, ex, ey,
                                          ez, pref;
  Cell *cell;
  int c, np;
  Particle *part;
  int i, j, id;
  double globalmax = 0;
  double f1, f2 = 0;

  iccp3m_sanity_check();

  if ((iccp3m_cfg.eout <= 0)) {
    ostringstream msg;
    msg << "ICCP3M: nonpositive dielectric constant is not allowed. Put a "
           "decent exception here\n";
    runtimeError(msg);
  }

  pref = 1.0 / (coulomb.prefactor * 6.283185307);
  iccp3m_cfg.citeration = 0;

  for (j = 0; j < iccp3m_cfg.num_iteration; j++) {
    double hmax = 0.;

    force_calc_iccp3m(); /* Calculate electrostatic forces (SR+LR) excluding
                            source source interaction*/
    ghost_communicator(&cell_structure.collect_ghost_force_comm);

    diff = 0;
    for (c = 0; c < local_cells.n; c++) {
      cell = local_cells.cell[c];
      part = cell->part;
      np = cell->n;
      for (i = 0; i < np; i++) {
        if (part[i].p.identity < iccp3m_cfg.n_ic + iccp3m_cfg.first_id &&
            part[i].p.identity >= iccp3m_cfg.first_id) {
          id = part[i].p.identity - iccp3m_cfg.first_id;
          /* the dielectric-related prefactor: */
          del_eps = (iccp3m_cfg.ein[id] - iccp3m_cfg.eout) /
                    (iccp3m_cfg.ein[id] + iccp3m_cfg.eout);
          /* calculate the electric field at the certain position */
          ex = part[i].f.f[0] / part[i].p.q;
          ey = part[i].f.f[1] / part[i].p.q;
          ez = part[i].f.f[2] / part[i].p.q;
          /* let's add the contribution coming from the external field */
          ex += iccp3m_cfg.extx;
          ey += iccp3m_cfg.exty;
          ez += iccp3m_cfg.extz;
          if (ex == 0 && ey == 0 && ez == 0) {
            ostringstream msg;
            msg << "ICCP3M found zero electric field on a charge. This must "
                   "never happen";
            runtimeError(msg);

          }
          /* the dot product   */
          fdot = ex * iccp3m_cfg.nvectorx[id] + ey * iccp3m_cfg.nvectory[id] +
                 ez * iccp3m_cfg.nvectorz[id];
          /* recalculate the old charge density */
          hold = part[i].p.q / iccp3m_cfg.areas[id];
          /* determine if it is higher than the previously highest charge
           * density */
          if (fabs(hold) > hmax)
            hmax = fabs(hold);
          f1 = del_eps * fdot * pref;
          if (iccp3m_cfg.sigma != 0)
            f2 = (2 * iccp3m_cfg.eout) /
                 (iccp3m_cfg.eout + iccp3m_cfg.ein[id]) *
                 (iccp3m_cfg.sigma[id]);
          /* relative variation: never use an estimator which can be negative
           * here */
          hnew =
              (1. - iccp3m_cfg.relax) * hold + (iccp3m_cfg.relax) * (f1 + f2);
          difftemp = fabs(1 * (hnew - hold) / (hmax + fabs(hnew + hold)));
          // if(difftemp > diff && part[i].p.q > 1e-5)
          if (difftemp > diff)
            diff =
                difftemp; /* Take the largest error to check for convergence */
          part[i].p.q = hnew * iccp3m_cfg.areas[id];

          /* check if the charge now is more than 1e6, to determine if ICC still
           * leads to reasonable results */
          /* this is kind a arbitrary measure but, does a good job spotting
           * divergence !*/
          if (fabs(part[i].p.q) > 1e6) {
            ostringstream msg;
            msg << "too big charge assignment in iccp3m! q >1e6 , assigned "
                   "charge= "
                << part[i].p.q << "\n";
            runtimeError(msg);
            diff = 1e90; /* A very high value is used as error code */
            break;
          }
        }
      } /* cell particles */
    }   /* local cells */
    iccp3m_cfg.citeration++;
    MPI_Allreduce(&diff, &globalmax, 1, MPI_DOUBLE, MPI_MAX, comm_cart);

    if (globalmax < iccp3m_cfg.convergence)
      break;
    if (diff > 1e89) {
      return iccp3m_cfg.citeration++;
    }

    /* Update charges on ghosts. */
    ghost_communicator(&cell_structure.exchange_ghosts_comm);
  } /* iteration */

  if (globalmax > iccp3m_cfg.convergence) {
    runtimeErrorMsg()
        << "ICC failed to converge in the given number of maximal steps.";
  }

  on_particle_charge_change();

  return iccp3m_cfg.citeration;
}

void force_calc_iccp3m() {
  init_forces_iccp3m();

  short_range_loop(Utils::NoOp{}, [](Particle &p1, Particle &p2, Distance &d) {
    /* calc non bonded interactions */
    add_non_bonded_pair_force_iccp3m(&(p1), &(p2), d.vec21.data(), sqrt(d.dist2),
                                     d.dist2);
  });

  calc_long_range_forces_iccp3m();
}

void init_forces_iccp3m() {
  for(auto & p : local_cells.particles()) {
    p.f = ParticleForce{};
  }

  for(auto & p : ghost_cells.particles()) {
    p.f = ParticleForce{};
  }
}

void calc_long_range_forces_iccp3m() {
#ifdef ELECTROSTATICS
  /* calculate k-space part of electrostatic interaction. */
  if (!(coulomb.method == COULOMB_ELC_P3M ||
        coulomb.method == COULOMB_P3M_GPU || coulomb.method == COULOMB_P3M ||
        coulomb.method == COULOMB_MMM2D || coulomb.method == COULOMB_MMM1D)) {
    ostringstream msg;
    msg << "ICCP3M implemented only for MMM1D,MMM2D,ELC or P3M ";
    runtimeError(msg);
  }
  switch (coulomb.method) {
#ifdef P3M
  case COULOMB_ELC_P3M:
    if (elc_params.dielectric_contrast_on) {
      ostringstream msg;
      msg << "ICCP3M conflicts with ELC dielectric constrast";
      runtimeError(msg);
    }
    p3m_charge_assign();
    p3m_calc_kspace_forces(1, 0);
    ELC_add_force();
    break;

#ifdef CUDA
  case COULOMB_P3M_GPU:
    if (this_node == 0) {
      FORCE_TRACE(printf("Computing GPU P3M forces.\n"));
      p3m_gpu_add_farfield_force();
    }
    break;
#endif
  case COULOMB_P3M:
    p3m_charge_assign();
    p3m_calc_kspace_forces(1, 0);
    break;
#endif
  case COULOMB_MMM2D:
    MMM2D_add_far_force();
    MMM2D_dielectric_layers_force_contribution();
    break;
  default:
    break;
  }

#endif
}

/** \name Private Functions */
/************************************************************/
/*@{*/

int iccp3m_sanity_check() {
  switch (coulomb.method) {
#ifdef P3M
  case COULOMB_ELC_P3M: {
    if (elc_params.dielectric_contrast_on) {
      ostringstream msg;
      msg << "ICCP3M conflicts with ELC dielectric constrast";
      runtimeError(msg);
      return 1;
    }
    break;
  }
#endif
  case COULOMB_DH: {
    ostringstream msg;
    msg << "ICCP3M does not work with Debye-Hueckel iccp3m.h";
    runtimeError(msg);
    return 1;
  }
  case COULOMB_RF: {
    ostringstream msg;
    msg << "ICCP3M does not work with COULOMB_RF iccp3m.h";
    runtimeError(msg);
    return 1;
  }
  default:
    break;
  }

#ifdef NPT
  if (integ_switch == INTEG_METHOD_NPT_ISO) {
    ostringstream msg;
    msg << "ICCP3M does not work in the NPT ensemble";
    runtimeError(msg);
    return 1;
  }
#endif

  return 0;
}

#endif
