/*
  Copyright (C) 2010,2011,2012,2013 The ESPResSo project
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

/** \file iccp3m.hpp 

    ICCP3M is a method that allows to take into account the influence
    of arbitrarliy shaped dielectric interfaces.  The dielectric
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
    acting on the induced charges has to be determined. As P3M an the
    other coulomb solvers calculate all mutual forces, the force
    calculation was modified to avoid the calculation of the short
    range part of the source-source force calculation.  For different
    particle data organisation schemes this is performed differently.
    */

#ifndef _ICCP3M_H 
#define _ICCP3M_H

#include <ctime>
#include "p3m.hpp"
#include "utils.hpp"
#include "mmm1d.hpp"
#include "mmm2d.hpp"
#include "domain_decomposition.hpp"
#include "cells.hpp"
#include "integrate.hpp"
#include "verlet.hpp"
#include "layered.hpp"
#include "global.hpp"
#include "ghosts.hpp"
#include "nsquare.hpp"
#include "interaction_data.hpp"
#include "particle_data.hpp"
#include "topology.hpp"
#include "ghosts.hpp"

#if defined(ELECTROSTATICS)

/* iccp3m data structures*/
typedef struct {
  int n_ic;                             /* Last induced id (can not be smaller then 2) */
  int num_iteration;                    /* Number of max iterations                    */
  double eout;                          /* Dielectric constant of the bulk             */
  double *areas;                        /* Array of area of the grid elements          */
  double *ein;                          /* Array of dielectric constants at each surface element */
  double *sigma;                        /* Surface Charge density */
  double convergence;                   /* Convergence criterion                       */
  double *nvectorx,*nvectory,*nvectorz; /* Surface normal vectors                      */
  double extx,exty,extz;             /* External field                              */
  double relax;                         /* relaxation parameter for iterative                       */
  int citeration ;                      /* current number of iterations*/
  int set_flag;                         /* flag that indicates if ICCP3M has been initialized properly */    
  double *fx;
  double *fy;
  double *fz;
  int first_id;
} iccp3m_struct;
extern iccp3m_struct iccp3m_cfg;        /* global variable with ICCP3M configuration */
extern int iccp3m_initialized;
#include "forces.hpp"

int bcast_iccp3m_cfg(void);

/** Calculation of the electrostatic forces between source charges (= real charges) and wall charges.
 *  For each electrostatic method the proper functions for short and long range parts are called.
 *  Long Range Parts are calculated directly, short range parts need helper functions according
 *  to the particle data organisation. A modified version of \ref force_calc in \ref forces.hpp.
 */
void force_calc_iccp3m();

/** Calculation of short range part of electrostatic interaction in layered systems 
 *  A modified version of \ref layered_calculate_ia in \ref layered.hpp
 */
void layered_calculate_ia_iccp3m();

/** Calculate the short range part of electrostatic interaction using verlet lists.
 * A modified version of \ref build_verlet_lists_and_calc_verlet_ia()	in \ref verlet.hpp
 */
void build_verlet_lists_and_calc_verlet_ia_iccp3m();

/** Calculate he short range part of electrostatic interaction using verlet lists, if verlet lists
 * have already been properly built. A modified version of \ref calculate_verlet_ia()	in \ref verlet.hpp
 */
void calculate_verlet_ia_iccp3m();

/** Calculate he short range part of electrostatic interaction using for linked cell
 * systems. A modified version of \ref calc_link_cell() in \ref domain_decomposition.hpp
 */
void calc_link_cell_iccp3m();

/** Calculate he short range part of electrostatic interaction using for n-squared cell
 * systems. A modified version of nsq_calculate_ia \ref nsquare.hpp.
 */
void nsq_calculate_ia_iccp3m();

/** The main iterative scheme, where the surface element charges are calculated self-consistently. 
 */
int iccp3m_iteration();

/** The initialisation of ICCP3M with zero values for all variables 
 */
void iccp3m_init(void);

/** check sanity of parameters for use with ICCP3M 
 */
int iccp3m_sanity_check();

/** Variant of add_non_bonded_pair_force where only coulomb 
 *  contributions are calculated   */
inline void add_non_bonded_pair_force_iccp3m(Particle *p1, Particle *p2, 
					double d[3], double dist, double dist2)
{
  /* IA_parameters *ia_params = get_ia_param(p1->p.type,p2->p.type);*/
  double force[3] = { 0, 0, 0 };
  int j;

  FORCE_TRACE(fprintf(stderr, "%d: interaction %d<->%d dist %f\n", this_node, p1->p.identity, p2->p.identity, dist));

  /***********************************************/
  /* long range electrostatics                   */
  /***********************************************/

  /* real space coulomb */
  double q1q2 = p1->p.q*p2->p.q;
  switch (coulomb.method) {
#ifdef P3M
  case COULOMB_ELC_P3M:
    if (q1q2) p3m_add_pair_force(q1q2,d,dist2,dist,force); 
    break;
  case COULOMB_P3M_GPU:
  case COULOMB_P3M:
    if (q1q2) p3m_add_pair_force(q1q2,d,dist2,dist,force);
    break;
#endif /* P3M */
  case COULOMB_MMM1D:
    if (q1q2) add_mmm1d_coulomb_pair_force(q1q2,d,dist2,dist,force);
    break;
  case COULOMB_MMM2D:
    if (q1q2) add_mmm2d_coulomb_pair_force(q1q2,d,dist2,dist,force);
    break;
  case COULOMB_NONE:
    break;
  }

  /***********************************************/
  /* add total nonbonded forces to particle      */
  /***********************************************/
   for (j = 0; j < 3; j++) { 
      p1->f.f[j] += force[j];
      p2->f.f[j] -= force[j];
   }
   /***********************************************/
}
#endif /* ELECTROSTATICS */

#endif /* ICCP3M_H */
