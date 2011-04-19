/*
  Copyright (C) 2010,2011 The ESPResSo project
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

/** \file iccp3m.h 

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

#ifndef ICCP3M_H 
#define ICCP3M_H

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

#if defined(ELECTROSTATICS)

/* iccp3m data structures*/
typedef struct {
  int last_ind_id;                      /* Last induced id (can not be smaller then 2) */
  int num_iteration;                    /* Number of max iterations                    */
  double eout;                          /* Dielectric constant of the bulk             */
  double *areas;                        /* Array of area of the grid elements          */
  double *ein;                          /* Array of dielectric constants at each surface element */
  double convergence;                   /* Convergence criterion                       */
  double *nvectorx,*nvectory,*nvectorz; /* Surface normal vectors                      */
  double *extx,*exty,*extz;             /* External field                              */
  int selection;                        /* by default it is not selected, WHO KNOWS WHAT THAT MEANS?*/
  double relax;                         /* relaxation parameter for iterative                       */
  int update;                           /* iccp3m update interval, currently not used */
  double *fx,*fy,*fz;                   /* forces iccp3m will use*/ 
  int citeration ;                      /* current number of iterations*/
  int set_flag;                         /* flag that indicates if ICCP3M has been initialized properly */    
} iccp3m_struct;

extern iccp3m_struct iccp3m_cfg;        /* global variable with ICCP3M configuration */
extern int iccp3m_initialized;
int bcast_iccp3m_cfg(void);
/** Implementation of the tcl-command <br>
    iccp3m  { \<last_ind_id\> \<e1\> \<num_iteration\> \<convergence\> \<relaxation\> \<area\> \<normal_components\> \<e_in/e_out\>  [\<ext_field\>] | iterate } 
    ICC sets up and calculates induced charges on dielectric surfaces. At the beginning of every simulation run particles on the surface boundary 
    have to be set up (before any real particle) together with the list of areas, normal vectors and dielectric constant associated with them. 
    After that the iterate flag can be used during the simulation to update the value of the induced charges.
    
    Parameters: <br>
                 \<last_ind_id\> ID of the last surface charge. Note that the IDs of the surface charges must range from 0 to \<last_ind_id\>
                 \<e1\>          = Dielectric Constant of the Bulk accessible to free particles 
                 \<num_iteration\> = Maximum number of ICCP3M iterations calculating the induced charges. 
                 \<relaxation\> = Relaxaxion parameter \f$omega\f$ for the successive over-relaxation scheme. 
                 \<area\>       = List of the areas of each surface element.
                 \<normal_components\> = List of normal vectors of each surface element. 3n list entries. Do not have to be normalized.
                 \<e_in/e_out\> = Ratio of dielectric co


                 iterate         = Indicates that a previous surface discretization shall be used. T
*/
int tclcommand_iccp3m(ClientData data, Tcl_Interp *interp, int argc, char **argv);

/** Calculation of the electrostatic forces between source charges (= real charges) and wall charges.
 *  For each electrostatic method the proper functions for short and long range parts are called.
 *  Long Range Parts are calculated directly, short range parts need helper functions according
 *  to the particle data organisation. A modified version of \ref force_calc in \ref forces.h.
 */
void force_calc_iccp3m();

/** Calculation of short range part of electrostatic interaction in layered systems 
 *  A modified version of \ref layered_calculate_ia in \ref layered.h
 */
void layered_calculate_ia_iccp3m();

/** Calculate the short range part of electrostatic interaction using verlet lists.
 * A modified version of \ref build_verlet_lists_and_calc_verlet_ia()	in \ref verlet.h
 */
void build_verlet_lists_and_calc_verlet_ia_iccp3m();

/** Calculate he short range part of electrostatic interaction using verlet lists, if verlet lists
 * have already been properly built. A modified version of \ref calculate_verlet_ia()	in \ref verlet.h
 */
void calculate_verlet_ia_iccp3m();

/** Calculate he short range part of electrostatic interaction using for linked cell
 * systems. A modified version of \ref calc_link_cell() in \ref domain_decomposition.h
 */
void calc_link_cell_iccp3m();

/** Calculate he short range part of electrostatic interaction using for n-squared cell
 * systems. A modified version of nsq_calculate_ia \ref nsquare.h.
 */
void nsq_calculate_ia_iccp3m();

/** The main iterative scheme, where the surface element charges are calculated self-consistently. 
 */
int iccp3m_iteration();

/** The initialisation of ICCP3M with zero values for all variables 
 */
void iccp3m_init(void);


/** The short range part of the electrostatic interation between two particles.  
 *  The appropriate function from the underlying electrostatic method is called. */
MDINLINE void add_non_bonded_pair_force_iccp3m(Particle *p1, Particle *p2, 
					double d[3], double dist, double dist2)
{
  /* IA_parameters *ia_params = get_ia_param(p1->p.type,p2->p.type);*/
  double force[3] = { 0, 0, 0 };
  int j;
  char* errtxt;

  FORCE_TRACE(fprintf(stderr, "%d: interaction %d<->%d dist %f\n", this_node, p1->p.identity, p2->p.identity, dist));

  /***********************************************/
  /* short range electrostatics                  */
  /***********************************************/

  if (coulomb.method == COULOMB_DH) {
    errtxt = runtime_error(128);
	  ERROR_SPRINTF(errtxt, "ICCP3M does not work with Debye-Hueckel iccp3m.h");
  }
  if (coulomb.method == COULOMB_RF) {
    errtxt = runtime_error(128);
	  ERROR_SPRINTF(errtxt, "ICCP3M does not work with COULOMB_RF iccp3m.h");
  }

  /*********************************************************************/
  /* everything before this contributes to the virial pressure in NpT, */
  /* but nothing afterwards                                            */
  /*********************************************************************/
#ifdef NPT
  for (j = 0; j < 3; j++)
    if(integ_switch == INTEG_METHOD_NPT_ISO)
  if (coulomb.method == COULOMB_RF) {
    errtxt = runtime_error(128);
	  ERROR_SPRINTF(errtxt, "ICCP3M does not work in the NPT ensemble");
  }
#endif

  /***********************************************/
  /* long range electrostatics                   */
  /***********************************************/

  /* real space coulomb */
  switch (coulomb.method) {

    /*  if ELCP3M */
    #ifdef ELC_P3M
    case COULOMB_ELC_P3M: {
      add_p3m_coulomb_pair_force(p1->p.q*p2->p.q,d,dist2,dist,force); 
      if (elc_params.dielectric_contrast_on) {
                  errtxt = runtime_error(128);
                  ERROR_SPRINTF(errtxt, "{ICCP3M conflicts with ELC dielectric constrast");
      }
      break;

    }
    #endif /* ELCP3M */

    /*  if P3M */
    #ifdef ELP3M
      case COULOMB_P3M: {
  
      #ifdef NPT /* ICCP3M does not work in NPT ensemble */
          if(integ_switch == INTEG_METHOD_NPT_ISO){
              errtxt = runtime_error(128);
              ERROR_SPRINTF(errtxt, "{ICCP3M cannot be used with pressure coupling} ");
          }
      #endif
  
      #ifdef DIPOLES /* ICCP3M still does not work with dipoles, so abort if compiled in */
            errtxt = runtime_error(128);
            ERROR_SPRINTF(errtxt, "{ICCP3M and dipoles not implemented yet} ");
      
      #else        /* If it is just electrostatic P3M */
          add_p3m_coulomb_pair_force(p2->p.q* p1->p.q,d,dist2,dist,force);
      #endif  /* DIPOLES */

    break;
    }
    #endif /* P3M */

    case COULOMB_MMM1D:
      add_mmm1d_coulomb_pair_force(p1,p2,d,dist2,dist,force);
      break;
    case COULOMB_MMM2D:
      add_mmm2d_coulomb_pair_force(p1->p.q*p2->p.q,d,dist2,dist,force);
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
