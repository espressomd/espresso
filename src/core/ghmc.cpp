/*
  Copyright (C) 2010,2012,2013,2014,2015,2016 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
  Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021
  Mainz, Germany

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
/** \file ghmc.cpp

    For more information see \ref ghmc.hpp
 */

#include "integrate.hpp"
#include "thermostat.hpp"
#include "utils.hpp"
#include <cmath>

#include "cells.hpp"
#include "communication.hpp"
#include "energy.hpp"
#include "ghmc.hpp"
#include "particle_data.hpp"
#include "random.hpp"
#include "statistics.hpp"
#include "virtual_sites.hpp"
#include "global.hpp"

/************************************************************/

// momentum flip flag
int ghmc_mflip = GHMC_MFLIP_OFF;
// temperature scaling flag
int ghmc_tscale = GHMC_TSCALE_OFF;
// result of mc decision
int ghmc_mc_res;

#ifdef GHMC
// MC statistics variables
int ghmc_att, ghmc_acc;
// partial momentum update parameters
double cosp, sinp;
// inverse temperature
double beta;
// thermostat data struct
Ghmc ghmcdata = {0, 0, 0.0, 0.0};
#endif

/************************************************************/
/* local prototypes                                         */
/************************************************************/

void hamiltonian_calc(int ekin_update_flag);

double calc_local_temp();

void calc_kinetic(double *ek_trans, double *ek_rot);

void simple_momentum_update();

void partial_momentum_update();

void tscale_momentum_update();

void momentum_flip();
/************************************************************/

/************************************************************/
/** \name Privat Functions */
/************************************************************/
/*@{*/

#ifdef GHMC

void save_last_state() {
  for (auto &p : local_cells.particles()) {
    memmove(&p.l.r_ls, &p.r, sizeof(ParticlePosition));
    memmove(&p.l.m_ls, &p.m, sizeof(ParticleMomentum));
  }
}

void load_last_state() {
  for (auto &p : local_cells.particles()) {
    memmove(&p.r, &p.l.r_ls, sizeof(ParticlePosition));
    memmove(&p.m, &p.l.m_ls, sizeof(ParticleMomentum));
  }
}

void hamiltonian_calc(int ekin_update_flag) {

  /* if ekin_update_flag = 0, we calculate all energies with \ref
   master_energy_calc().
   if ekin_update_flag = 1, we only updated momenta, so there we only need to
   recalculate
   kinetic energy with \ref calc_kinetic(). */

  double ekt, ekr;

  INTEG_TRACE(fprintf(stderr, "%d: hamiltonian_calc:\n", this_node));

  // initialize energy struct
  if (total_energy.init_status == 0) {
    init_energies(&total_energy);
    // if we are initializing energy we have to calculate all energies anyway
    ekin_update_flag = 0;
  }

  // calculate energies
  if (ekin_update_flag == 0)
    master_energy_calc();
  else
    calc_kinetic(&ekt, &ekr);

  // sum up energies on master node, and update ghmcdata struct
  if (this_node == 0) {
    double result = 0.0;
    ghmcdata.hmlt_old = ghmcdata.hmlt_new;
    for (int i = ekin_update_flag; i < total_energy.data.n; i++) {
      result += total_energy.data.e[i];
    }
    if (ekin_update_flag == 1)
      result += ekt + ekr;
    ghmcdata.hmlt_new = result;
  }
}

// get local temperature - here for debbuging purposes
double calc_local_temp() {
  int tot_np = 0;
  double temp = 0.0;

  for (auto &p : local_cells.particles()) {
    for (int j = 0; j < 3; j++) {
      temp += p.p.mass * Utils::sqr(p.m.v[j] / time_step);
#ifdef ROTATION
      temp += Utils::sqr(p.m.omega[j]);
#endif
    }
  }
#ifdef ROTATION
  temp /= 6 * tot_np;
#else
  temp /= 3 * tot_np;
#endif
  return temp;
}

void calc_kinetic(double *ek_trans, double *ek_rot) {
  double et = 0.0, er = 0.0;

  for (auto &p : local_cells.particles()) {
#ifdef VIRTUAL_SITES
    if (p.p.is_virtual)
      continue;
#endif

    /* kinetic energy */
    et += (Utils::sqr(p.m.v[0]) + Utils::sqr(p.m.v[1]) + Utils::sqr(p.m.v[2])) * (p).p.mass;

/* rotational energy */
#ifdef ROTATION
    er += Utils::sqr(p.m.omega[0]) * p.p.rinertia[0] +
          Utils::sqr(p.m.omega[1]) * p.p.rinertia[1] +
          Utils::sqr(p.m.omega[2]) * p.p.rinertia[2];
#endif
  }

  MPI_Allreduce(MPI_IN_PLACE, &et, 1, MPI_DOUBLE, MPI_SUM, comm_cart);
  MPI_Allreduce(MPI_IN_PLACE, &er, 1, MPI_DOUBLE, MPI_SUM, comm_cart);

  et /= (2.0 * time_step * time_step);
#ifdef ROTATION
  er /= 2.0;
#endif

  *ek_trans = et;
  *ek_rot = er;
}

/* momentum update step of ghmc */
void ghmc_momentum_update() {

  INTEG_TRACE(fprintf(stderr, "%d: ghmc_momentum_update:\n", this_node));

  /* currently, when temperature scaling is enabled the procedure
   tscale_momentum_update is used,
   although it may seem like a code duplication, this separation is maintained
   for now until this feature is tested thoroughly*/

  if (ghmc_tscale == GHMC_TSCALE_ON)
    tscale_momentum_update();
  else
    partial_momentum_update();

  int ekin_update_flag = 1;
  hamiltonian_calc(ekin_update_flag);
}

/* momentum update step of ghmc with temperature scaling */
void tscale_momentum_update() {
  save_last_state();

  simple_momentum_update();

  double tempt, tempr;
  calc_kinetic(&tempt, &tempr);
  tempt /= (1.5 * n_part);
  tempr /= (1.5 * n_part);

  double scalet = sqrt(temperature / tempt);
#ifdef ROTATION
  double scaler = sqrt(temperature / tempr);
#endif

  for (auto &p : local_cells.particles()) {
    for (int j = 0; j < 3; j++) {
      p.m.v[j] *= scalet;
#ifdef ROTATION
      p.m.omega[j] *= scaler;
#endif
    }
  }

  for (auto &p : local_cells.particles()) {
#ifdef VIRTUAL_SITES
    if (p.p.is_virtual)
      continue;
#endif

    for (int j = 0; j < 3; j++) {
      p.m.v[j] = cosp * (p.l.m_ls.v[j]) + sinp * (p.m.v[j]);
#ifdef ROTATION
      p.m.omega[j] = cosp * (p.l.m_ls.omega[j]) + sinp * (p.m.omega[j]);
#endif
    }
  }
}

/* momentum update step of ghmc */
void simple_momentum_update() {
  double sigmat, sigmar;

  sigmat = sqrt(temperature);
  sigmar = sqrt(temperature);

  for (auto &p : local_cells.particles()) {
#ifdef VIRTUAL_SITES
    if (p.p.is_virtual)
      continue;
#endif

#ifdef MASS
    sigmat = sqrt(temperature / (p).p.mass);
#endif
    for (int j = 0; j < 3; j++) {
      if (sigmat > 0.0) {
        p.m.v[j] = sigmat * gaussian_random() * time_step;
      }
#ifdef ROTATION
#ifdef ROTATIONAL_INERTIA
      sigmar = sqrt(temperature / p.p.rinertia[j]);
#endif
      if (sigmar > 0.0) {
        p.m.omega[j] = sigmar * gaussian_random();
      }
#endif
    }
  }
}

/* partial momentum update step for ghmc */
void partial_momentum_update() {
  double sigmat, sigmar;

  sigmat = sqrt(temperature);
  sigmar = sqrt(temperature);

  for (auto &p : local_cells.particles()) {
#ifdef MASS
    sigmat = sqrt(temperature / (p).p.mass);
#endif
    for (int j = 0; j < 3; j++) {
      if (sigmat > 0.0) {
        p.m.v[j] =
          cosp * (p.m.v[j]) + sinp * (sigmat * gaussian_random() * time_step);
      } else {
        p.m.v[j] = cosp * p.m.v[j];
      }
#ifdef ROTATION
#ifdef ROTATIONAL_INERTIA
      sigmar = sqrt(temperature / p.p.rinertia[j]);
#endif
      if (sigmar > 0.0) {
        p.m.omega[j] =
          cosp * (p.m.omega[j]) + sinp * (sigmar * gaussian_random());
      } else {
        p.m.omega[j] = cosp * p.m.omega[j];
      }
#endif
    }
  }
}

/* momentum flip for ghmc */
void momentum_flip() {
  INTEG_TRACE(fprintf(stderr, "%d: ghmc_momentum_flip:\n", this_node));

  for (auto &p : local_cells.particles()) {
    for (int j = 0; j < 3; j++) {
      p.m.v[j] = -(p.m.v[j]);
#ifdef ROTATION
      p.m.omega[j] = -(p.m.omega[j]);
#endif
    }
  }
}
/*@}*/

#endif

/************************************************************/
/** \name Exported Functions */
/************************************************************/
/*@{*/

#ifdef GHMC

void thermo_init_ghmc() {
  ghmc_att = 0;
  ghmc_acc = 0;

  cosp = ghmc_phi;
  sinp = sin(acos(ghmc_phi));

  ghmcdata.hmlt_old = 0;
  ghmcdata.hmlt_new = 0;
  beta = 1.0 / temperature;

  THERMO_TRACE(fprintf(stderr,
                       "%d: thermo_init_ghmc: ghmc_csp=%f, ghmc_snp=%f \n",
                       this_node, cosp, sinp));
}

void ghmc_init() {
  INTEG_TRACE(fprintf(stderr, "%d: ghmc_init:\n", this_node));

  ghmcdata.att = 0;
  ghmcdata.acc = 0;

  save_last_state();
}

void ghmc_close() {

  INTEG_TRACE(fprintf(stderr, "%d: ghmc_close:\n", this_node));
  ghmc_att += ghmcdata.att;
  ghmc_acc += ghmcdata.acc;
}

/* monte carlo step of ghmc - evaluation stage */
void ghmc_mc() {
  INTEG_TRACE(fprintf(stderr, "%d: ghmc_mc:\n", this_node));

  int ekin_update_flag = 0;
  hamiltonian_calc(ekin_update_flag);

  // make MC decision only on the master
  if (this_node == 0) {

    ghmcdata.att++;

    // metropolis algorithm
    double boltzmann = ghmcdata.hmlt_new - ghmcdata.hmlt_old;
    if (boltzmann < 0)
      boltzmann = 1.0;
    else if (boltzmann > 30)
      boltzmann = 0.0;
    else
      boltzmann = exp(-beta * boltzmann);

    if (d_random() < boltzmann) {
      ghmcdata.acc++;
      ghmc_mc_res = GHMC_MOVE_ACCEPT;
    } else {
      ghmc_mc_res = GHMC_MOVE_REJECT;
    }
  }

  // let all nodes know about the MC decision result
  mpi_bcast_parameter(FIELD_GHMC_RES);

  if (ghmc_mc_res == GHMC_MOVE_ACCEPT) {
    save_last_state();
    // fprintf(stderr,"%d: mc move accepted\n",this_node);
  } else {
    load_last_state();
    // fprintf(stderr,"%d: mc move rejected\n",this_node);

    // if the move is rejected we might need to resort particles according to
    // the loaded configurations
    cells_resort_particles(CELL_GLOBAL_EXCHANGE);
    invalidate_obs();

    if (ghmc_mflip == GHMC_MFLIP_ON) {
      momentum_flip();
    } else if (ghmc_mflip == GHMC_MFLIP_RAND) {
      if (d_random() < 0.5)
        momentum_flip();
    }
  }

  // fprintf(stderr,"%d: temp after mc: %f\n",this_node,calc_local_temp());
}

/*@}*/

#endif
