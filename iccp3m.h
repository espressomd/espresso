// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2006; all rights reserved unless otherwise stated.
/** \file iccp3m.h
    For more information on particle_data,
    see \ref iccp3m.c "iccp3m.c"

 */

#include <tcl.h>
#include <time.h>
#include "p3m.h"
#include "utils.h"
#include "mmm1d.h"
#include "mmm2d.h"
#include "domain_decomposition.h"
#include "cells.h"
#include "integrate.h"
#include "communication.h"
#include "verlet.h"
#include "layered.h"
#include "global.h"
#include "communication.h"
#include "ghosts.h"
#include "nsquare.h"
#include "interaction_data.h"
#include "particle_data.h"
#include "topology.h"
#include "forces.h"
#include "ghosts.h"

/* iccp3m data structures*/
typedef struct {
  int last_ind_id;   /* Last induced id can not be smaller then 2 */
  int num_iteration; /* Number of max iterations                  */
  double e1;         /* Dielectric constants                      */
  double e2; 
  double area;        /* Area of the grid element */
  double *areas;      /* Array of area of the grid elements */
  double convergence; /* Convergence criterion */
  double *nvectorx,*nvectory,*nvectorz; /* Surface normal vectors */
  int selection;  /* by default it is not selected*/
  double relax;   /* relaxation parameter for iterative */
  int update;         /* iccp3m update interval */
  double *fx,*fy,*fz; /* forces iccp3m will use*/ 
  int citeration ; /* current number of iterations*/
} iccp3m_struct;

iccp3m_struct iccp3m_cfg;
/* parse iccp3m options */
int iccp3m(ClientData data, Tcl_Interp *interp, int argc, char **argv);
/* iccp3m : get trial electrostatic forces from p3m and related functions*/
void force_calc_iccp3m();
void layered_calculate_ia_iccp3m();
void build_verlet_lists_and_calc_verlet_ia_iccp3m();
void calculate_verlet_ia_iccp3m();
void calc_link_cell_iccp3m();
void nsq_calculate_ia_iccp3m();
/* iterative sceme to obtain induced charges */
int iccp3m_iteration();

/* pass citeration */
int inter_parse_iccp3m(Tcl_Interp * interp, int argc, char ** argv);
int gettime(ClientData data, Tcl_Interp *interp, int argc, char **argv);

/* add non-bonded iccp3m */
MDINLINE void add_non_bonded_pair_force_iccp3m(Particle *p1, Particle *p2, 
					double d[3], double dist, double dist2)
{
  /* IA_parameters *ia_params = get_ia_param(p1->p.type,p2->p.type);*/
  double force[3] = { 0, 0, 0 };
  int j;

  FORCE_TRACE(fprintf(stderr, "%d: interaction %d<->%d dist %f\n", this_node, p1->p.identity, p2->p.identity, dist));

  /***********************************************/
  /* short range electrostatics                  */
  /***********************************************/

#ifdef ELECTROSTATICS
  if (coulomb.method == COULOMB_DH)
    add_dh_coulomb_pair_force(p1,p2,d,dist,force);
#endif

  /*********************************************************************/
  /* everything before this contributes to the virial pressure in NpT, */
  /* but nothing afterwards                                            */
  /*********************************************************************/
#ifdef NPT
  for (j = 0; j < 3; j++)
    if(integ_switch == INTEG_METHOD_NPT_ISO)
      nptiso.p_vir[j] += force[j] * d[j];
#endif

  /***********************************************/
  /* long range electrostatics                   */
  /***********************************************/

#ifdef ELECTROSTATICS

  /* real space coulomb */
  switch (coulomb.method) {
  case COULOMB_P3M: {
#ifdef NPT
    double eng = add_p3m_coulomb_pair_force(p1->p.q*p2->p.q,d,dist2,dist,force);
    if(integ_switch == INTEG_METHOD_NPT_ISO)
      nptiso.p_vir[0] += eng;
#else
    /* add_p3m_coulomb_pair_force(p1,p2,d,dist2,dist,force); 
     * this was the old form, simply changed the first two parameters 
     * multiplication of charges */
    add_p3m_coulomb_pair_force(p2->p.q* p1->p.q,d,dist2,dist,force);
#endif
    break;
  }
  case COULOMB_MMM1D:
    add_mmm1d_coulomb_pair_force(p1,p2,d,dist2,dist,force);
    break;
  case COULOMB_MMM2D:
    add_mmm2d_coulomb_pair_force(p1->p.q*p2->p.q,d,dist2,dist,force);
    break;
  case COULOMB_MAGGS:
    if(maggs.yukawa == 1)
      add_maggs_yukawa_pair_force(p1,p2,d,dist2,dist,force);
    break;
  case COULOMB_NONE:
    break;
  }

#endif

  /***********************************************/
  /* add total nonbonded forces to particle      */
  /***********************************************/
   for (j = 0; j < 3; j++) { 
      p1->f.f[j] += force[j];
      p2->f.f[j] -= force[j];
     }
    /***********************************************/
  /*if(!(p1->p.identity > iccp3m_cfg.last_ind_id && p2->p.identity >iccp3m_cfg.last_ind_id)) {
    * printf("peha identities (%d,%d) \n",p1->p.identity,p2->p.identity);*
      iccp3m_cfg.fx[p1->p.identity]=force[0];
      iccp3m_cfg.fx[p2->p.identity]=-force[0];
      iccp3m_cfg.fy[p1->p.identity]=force[1];
      iccp3m_cfg.fy[p2->p.identity]=-force[1];
      iccp3m_cfg.fz[p1->p.identity]=force[2];
      iccp3m_cfg.fz[p2->p.identity]=-force[2]; 
  }*/
}
