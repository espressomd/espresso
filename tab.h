// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2003; all rights reserved unless otherwise stated.
#ifndef TAB_H
#define TAB_H

/** \file tab.h
 *  Routines to calculate the  energy and/or  force 
 *  for a particle pair via interpolating from lookup tables.
 *  \ref forces.c
 *
 *  <b>Responsible:</b>
 *  <a href="mailto:cooke@mpip-mainz.mpg.de">Ira</a>
*/

#include "utils.h"

/** Add a pair force by linear interpolation from a table */
MDINLINE void add_tabulated_pair_force(Particle *p1, Particle *p2, IA_parameters *ia_params,
				double d[3], double dist)
{
  double phi, dindex, fac;
  int tablepos, table_start,j;
  double rescaled_force_cap = tab_force_cap/dist;
  double maxval = ia_params->TAB_maxval;
  double minval = ia_params->TAB_minval;

  fac = 0.0;

  if ( maxval > 0 ) {
    if ( dist < maxval){ 
      table_start = ia_params->TAB_startindex;
      dindex = (dist-minval)/ia_params->TAB_stepsize;
      tablepos = (int)(floor(dindex));  

      if ( dist > minval ) {
       phi = dindex - tablepos;	  
       fac = tabulated_forces.e[table_start + tablepos]*(1-phi) + tabulated_forces.e[table_start + tablepos+1]*phi;
    }
    else {
      /* Use an extrapolation beyond the table */
      if ( dist > 0 ) {
	    tablepos = 0;
	    phi = dindex - tablepos;	  
	    fac = (tabulated_forces.e[table_start]*minval)*(1-phi) + 
	        (tabulated_forces.e[table_start+1]*(minval+ia_params->TAB_stepsize))*phi;
	    fac = fac/dist;
      }
      else { /* Particles on top of each other .. leave fac as 0.0 */
      }
    }
    
    if ( rescaled_force_cap < fac && tab_force_cap > 0.0) {
      fac = rescaled_force_cap;
    }
  }
  /* Now update the forces */  
    for(j=0;j<3;j++) {
      p1->f.f[j] += fac * d[j];
      p2->f.f[j] -= fac * d[j];
#ifdef NPT
      if (piston != 0.0) 
	p_vir += fac*d[j] * d[j];
#endif
    }
  }
}

/** Add a pair energy by linear interpolation from a table */
MDINLINE double tabulated_pair_energy(Particle *p1, Particle *p2, IA_parameters *ia_params,
				      double d[3], double dist) {
  double phi, dindex;
  int tablepos, table_start;


  if ( dist < ia_params->TAB_maxval){ 
    dindex = (dist-ia_params->TAB_minval)/ia_params->TAB_stepsize;
    tablepos = (int)(floor(dindex)); 

    if ( tablepos < 0 ){
      tablepos = 0; /* Just gives a linear extrapolation in the region outside the table */     
    }

    phi = (dindex - tablepos);
    table_start = ia_params->TAB_startindex;

    return  tabulated_energies.e[table_start + tablepos]*(1-phi) 
      + tabulated_energies.e[table_start + tablepos+1]*phi;
  }
  return 0.0;
}


/** check the tabulated forcecap to see that it is sensible \warning
    This routine will probably give strange results if forcecap is
    applied before the table is loaded */
MDINLINE void check_tab_forcecap(double force_cap)
{
  int i,j,startindex;
  IA_parameters *params;

  for(i=0; i<n_particle_types; i++) {
    for(j=0; j<n_particle_types; j++) {
      params = get_ia_param(i,j);
      startindex = params->TAB_startindex;
      if ( tabulated_forces.max < (params->TAB_npoints + startindex )) { /* Make sure forces are initialized */
	if(force_cap > 0.0 && params->TAB_maxval > 0.0 && tabulated_forces.e[startindex] > force_cap) {
	  for ( i = 0 ; i < params->TAB_npoints ; i++) {
	    if ( tabulated_forces.e[startindex + i] < force_cap ) {
	      return; /* Everything is OK nothing to say :) */
	    }	  
	  }
	  if ( i == params->TAB_npoints - 1) {
	    tab_force_cap = -1.0;
	    /* Force cap is below all known forces .. turn force capping off */
	  }
	}    
	if ( force_cap > tabulated_forces.e[startindex] ) {
	  fprintf(stderr,"calc_tab_cap_radii: Capped region is outside the force table");
	  
	  if ( tabulated_forces.e[startindex] < tabulated_forces.e[startindex+1] ) {
	    fprintf(stderr,"Extrapolation does not make sense outside this table \n");
	    errexit();
	  }
	  
	  fprintf(stderr,", will extrapolate the force outside table \n");
	  printf("fc: %f \n",tab_force_cap);	
	}
      }
    }
  }
}

#endif
