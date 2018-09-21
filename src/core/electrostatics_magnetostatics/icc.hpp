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
//

/** \file

    ICCP3M is a method that allows to take into account the influence
    of arbitrarily shaped dielectric interfaces.  The dielectric
    properties of a dielectric medium in the bulk of the simulation
    box are taken into account by reproducing the jump in the electric
    field at the inface with charge surface segments. The charge
    density of the surface segments have to be determined
    self-consistently using an iterative scheme.  It can at presently
    - despite its name - be used with P3M, ELCP3M, MMM2D and MMM1D.
    For details see:<br> S. Tyagi, M. Suzen, M. Sega, C. Holm,
    M. Barbosa: A linear-scaling method for computing induced charges
    on arbitrary dielectric boundaries in large system simulations
    (Preprint)

    To set up ICCP3M first the dielectric boundary has to be modelled
    by espresso particles 0..n where n has to be passed as a parameter
    to ICCP3M. This is still a bit inconvenient, as it forces the user
    to reserve the first n particle ids to wall charges, but as the
    other parts of espresso do not suffer from a limitation like this,
    it can be tolerated.

    For the determination of the induced charges only the forces
    acting on the induced charges has to be determined. As P3M and the
    other coulomb solvers calculate all mutual forces, the force
    calculation was modified to avoid the calculation of the short
    range part of the source-source force calculation.  For different
    particle data organisation schemes this is performed differently.
    */

#ifndef _ICCP3M_H
#define _ICCP3M_H

#include "config.hpp"

#if defined(ELECTROSTATICS)

#include "cells.hpp"
#include "electrostatics_magnetostatics/mmm1d.hpp"
#include "electrostatics_magnetostatics/mmm2d.hpp"
#include "electrostatics_magnetostatics/p3m.hpp"
#include "ghosts.hpp"
#include "global.hpp"
#include "integrate.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "particle_data.hpp"
#include "topology.hpp"
#include "utils.hpp"
#include <ctime>

/* iccp3m data structures*/
typedef struct {
  int initialized;
  int n_ic;          /* Last induced id (can not be smaller then 2) */
  int num_iteration; /* Number of max iterations                    */
  double eout;       /* Dielectric constant of the bulk             */
  double *areas;     /* Array of area of the grid elements          */
  double *ein;       /* Array of dielectric constants at each surface element */
  double *sigma;     /* Surface Charge density */
  double convergence; /* Convergence criterion                       */
  double *nvectorx;
  double *nvectory;
  double *nvectorz; /* Surface normal vectors                      */
  double extx;
  double exty;
  double extz;    /* External field                              */
  double relax;   /* relaxation parameter for iterative                       */
  int citeration; /* current number of iterations*/
  int set_flag;   /* flag that indicates if ICCP3M has been initialized properly
                   */

  double *fx;
  double *fy;
  double *fz;
  int first_id;
} iccp3m_struct;
extern iccp3m_struct iccp3m_cfg; /* global variable with ICCP3M configuration */
extern int iccp3m_initialized;

int bcast_iccp3m_cfg(void);

/** The main iterative scheme, where the surface element charges are calculated
 * self-consistently.
 */
int iccp3m_iteration();

/** The initialisation of ICCP3M with zero values for all variables
 */
void iccp3m_init(void);

/** The allocation of ICCP3M lists for python interface
 */
void iccp3m_alloc_lists();

/** The set of the init flag from python interface
 */
void iccp3m_set_initialized();

/** check sanity of parameters for use with ICCP3M
 */
int iccp3m_sanity_check();

#endif /* ELECTROSTATICS */

#endif /* ICCP3M_H */
