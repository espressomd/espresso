/*
  Copyright (C) 2010,2011,2012 The ESPResSo project
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
#ifndef _MORSE_H
#define _MORSE_H

/** \file morse.h
 *  Routines to calculate the lennard jones energy and/or  force 
 *  for a particle pair.
 *  \ref forces.c
*/

#include "utils.h"

#ifdef MORSE

#include "morse.h"
#include "interaction_data.h"
#include "communication.h"

/** set the force cap for the Morse interaction.
    @param morseforcecap the maximal force, 0 to disable, -1 for individual cutoff
    for each of the interactions.
*/

MDINLINE int morseforcecap_set_params(double morseforcecap)
{
  if (morse_force_cap != -1.0)
    mpi_morse_cap_forces(morse_force_cap);
  
  return TCL_OK;
}


MDINLINE int morse_set_params(int part_type_a, int part_type_b,
				      double eps, double alpha,
              double rmin, double cut, double cap_radius)
{
  double add1, add2;
  IA_parameters *data, *data_sym;

  make_particle_type_exist(part_type_a);
  make_particle_type_exist(part_type_b);
    
  data     = get_ia_param(part_type_a, part_type_b);
  data_sym = get_ia_param(part_type_b, part_type_a);

  if (!data || !data_sym) {
    return TCL_ERROR;
  }

  /* MORSE should be symmetrically */
  data->MORSE_eps       = data_sym->MORSE_eps       = eps;
  data->MORSE_alpha    = data_sym->MORSE_alpha    = alpha;
  data->MORSE_rmin    = data_sym->MORSE_rmin    = rmin;
  data->MORSE_cut    = data_sym->MORSE_cut    = cut;

  /* calculate dependent parameter */
  add1 = exp(-2.0*data->MORSE_alpha*(data->MORSE_cut - data->MORSE_rmin));
  add2 = 2.0*exp(-data->MORSE_alpha*(data->MORSE_cut - data->MORSE_rmin));
  data->MORSE_rest    = data_sym->MORSE_rest    = data->MORSE_eps * (add1 - add2); 
 
  if (cap_radius > 0) {
    data->MORSE_capradius = cap_radius;
    data_sym->MORSE_capradius = cap_radius;
  }

  /* broadcast interaction parameters */
  mpi_bcast_ia_params(part_type_a, part_type_b);
  mpi_bcast_ia_params(part_type_b, part_type_a);

  if (morse_force_cap != -1.0)
    mpi_morse_cap_forces(morse_force_cap);

  return TCL_OK;
}

/** Calculate Morse force between particle p1 and p2 */
MDINLINE void add_morse_pair_force(Particle *p1, Particle *p2, IA_parameters *ia_params,
				double d[3], double dist, double force[3])
{
  double add1, add2, fac=0.0;
  int j;
  if(dist < ia_params->MORSE_cut) { 
    /* normal case: resulting force/energy smaller than capping. */

    if(dist > ia_params->MORSE_capradius ) {
      add1 = exp(-2.0 * ia_params->MORSE_alpha * (dist - ia_params->MORSE_rmin));
      add2 = exp( -ia_params->MORSE_alpha * (dist - ia_params->MORSE_rmin));
      fac   = -ia_params->MORSE_eps * 2.0 * ia_params->MORSE_alpha * (add2 - add1) / dist;

      for(j=0;j<3;j++)
	force[j] += fac * d[j];
#ifdef MORSE_WARN_WHEN_CLOSE
      if(fac*dist > 1000) fprintf(stderr,"%d: Morse-Warning: Pair (%d-%d) force=%f dist=%f\n",
				  this_node,p1->p.identity,p2->p.identity,fac*dist,dist);
#endif

    }
    /* capped part of morse potential. */

    else if(dist > 0.0) {

        add1 = exp(-2.0 * ia_params->MORSE_alpha * (ia_params->MORSE_capradius - ia_params->MORSE_rmin));
        add2 = exp( -ia_params->MORSE_alpha * (ia_params->MORSE_capradius - ia_params->MORSE_rmin));
        fac   = -ia_params->MORSE_eps * 2.0 * ia_params->MORSE_alpha * (add2 - add1) / ia_params->MORSE_capradius;

	/* vector d is rescaled to length MORSE_capradius */
        for(j=0;j<3;j++)
	  force[j] += fac * d[j];
    }
    /* this should not happen! */

    else {
      MORSE_TRACE(fprintf(stderr, "%d: Morse warning: Particles id1=%d id2=%d exactly on top of each other\n",this_node,p1->p.identity,p2->p.identity));

      add1 = exp(-2.0 * ia_params->MORSE_alpha * (ia_params->MORSE_capradius - ia_params->MORSE_rmin));
      add2 = exp( -ia_params->MORSE_alpha * (ia_params->MORSE_capradius - ia_params->MORSE_rmin));
      fac   = -ia_params->MORSE_eps * 2.0 * ia_params->MORSE_alpha * (add2 - add1) / ia_params->MORSE_capradius;
      
      force[0] += fac * ia_params->MORSE_capradius;
    }

    ONEPART_TRACE(if(p1->p.identity==check_id) fprintf(stderr,"%d: OPT: MORSE   f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p1->f.f[0],p1->f.f[1],p1->f.f[2],p2->p.identity,dist,fac));
    ONEPART_TRACE(if(p2->p.identity==check_id) fprintf(stderr,"%d: OPT: MORSE   f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p2->f.f[0],p2->f.f[1],p2->f.f[2],p1->p.identity,dist,fac));

    MORSE_TRACE(fprintf(stderr,"%d: LJ: Pair (%d-%d) dist=%.3f: force+-: (%.3e,%.3e,%.3e)\n",
		     this_node,p1->p.identity,p2->p.identity,dist,fac*d[0],fac*d[1],fac*d[2]));

  }
}

/** calculate Morse energy between particle p1 and p2. */
MDINLINE double morse_pair_energy(Particle *p1, Particle *p2, IA_parameters *ia_params,
				double d[3], double dist)
{

  double add1, add2, fac;

  if(dist < ia_params->MORSE_cut) {

    /* normal case: resulting force/energy smaller than capping. */
    if(dist > ia_params->MORSE_capradius) { 
        add1 = exp(-2.0 * ia_params->MORSE_alpha * (dist - ia_params->MORSE_rmin));
        add2 = 2.0 * exp( -ia_params->MORSE_alpha * (dist - ia_params->MORSE_rmin));
        fac   = ia_params->MORSE_eps * (add1 - add2) - ia_params->MORSE_rest;
      return fac;
    }

    /* capped part of morse potential. */
    else if(dist > 0.0) {
        add1 = exp(-2.0 * ia_params->MORSE_alpha * (ia_params->MORSE_capradius - ia_params->MORSE_rmin));
        add2 = 2.0 * exp( -ia_params->MORSE_alpha * (ia_params->MORSE_capradius - ia_params->MORSE_rmin));
        fac   = ia_params->MORSE_eps * (add1 - add2) - ia_params->MORSE_rest;
      return fac;
    }
    /* this should not happen! */

    else {
        add1 = exp(-2.0 * ia_params->MORSE_alpha * (ia_params->MORSE_capradius - ia_params->MORSE_rmin));
        add2 = 2.0 * exp( -ia_params->MORSE_alpha * (ia_params->MORSE_capradius - ia_params->MORSE_rmin));
        fac   = ia_params->MORSE_eps * (add1 - add2) - ia_params->MORSE_rest;
        return fac;
    }

  }

  return 0.0;
}

/** calculate morse_capradius from morse_force_cap */

MDINLINE void calc_morse_cap_radii(double force_cap)
{
  int i,j,cnt=0;
  IA_parameters *params;
  double force=0.0, rad=0.0, step, add1, add2;

  for(i=0; i<n_particle_types; i++) {
    for(j=0; j<n_particle_types; j++) {
      params = get_ia_param(i,j);
      if(force_cap > 0.0 && params->MORSE_eps > 0.0) {

	/* I think we have to solve this numerically... and very crude as well */

	cnt=0;
	rad = params->MORSE_rmin - 0.69314719 / params->MORSE_alpha;
	step = -0.1 * rad;
	force=0.0;
	
	while(step != 0) {

          add1 = exp(-2.0 * params->MORSE_alpha * (rad - params->MORSE_rmin));
          add2 = exp( -params->MORSE_alpha * (rad - params->MORSE_rmin));
          force   = -params->MORSE_eps * 2.0 * params->MORSE_alpha * (add2 - add1) / rad;

	  if((step < 0 && force_cap < force) || (step > 0 && force_cap > force)) {
	    step = - (step/2.0); 
	  }
	  if(fabs(force-force_cap) < 1.0e-6) step=0;
	  rad += step; cnt++;
	} 
      	params->MORSE_capradius = rad;
      }
      else {
	params->MORSE_capradius = 0.0; 
      }
      FORCE_TRACE(fprintf(stderr,"%d: Ptypes %d-%d have cap_radius %f and cap_force %f (iterations: %d)\n",
			  this_node,i,j,rad,force,cnt));
    }
  }
}

#endif /* ifdef MORSE */
#endif /* ifdef _MORSE_H */
