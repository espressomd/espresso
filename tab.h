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
 *  for a particle pair or bonds via interpolating from lookup tables.
 *  \ref forces.c
 *  Needs feature TABULATED compiled in (see \ref config.h).

 *  <b>Responsible:</b>
 *  <a href="mailto:cooke@mpip-mainz.mpg.de">Ira</a>
*/


#include "dihedral.h"

/** Add a non-bonded pair force by linear interpolation from a table.
    Needs feature TABULATED compiled in (see \ref config.h). */
MDINLINE void add_tabulated_pair_force(Particle *p1, Particle *p2, IA_parameters *ia_params,
				       double d[3], double dist)
{
#ifdef TABULATED
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
      if(integ_switch == INTEG_METHOD_NPT_ISO)
	nptiso.p_vir[j] += fac*d[j] * d[j];
#endif
    }
  }
#endif
}

/** Add a non-bonded pair energy by linear interpolation from a table.
    Needs feature TABULATED compiled in (see \ref config.h). */
MDINLINE double tabulated_pair_energy(Particle *p1, Particle *p2, IA_parameters *ia_params,
				      double d[3], double dist) {
#ifdef TABULATED
  double phi, dindex;
  int tablepos, table_start;


  if ( dist < ia_params->TAB_maxval){ 
    dindex = (dist-ia_params->TAB_minval)/ia_params->TAB_stepsize;
    tablepos = (int)(floor(dindex)); 
    table_start = ia_params->TAB_startindex;

    if ( tablepos < 0 ){
      /* For distances smaller than the tabulated minimim take quadratic
	 extrapolation from first two values of the force table! 
	 This corresponds to the linear extrapolation of the force at that point.
	 This sould not occur too often, since it is quite expensive!
      */
      tablepos = 0;   
      double x0, b;
      b = (tabulated_forces.e[table_start + tablepos + 1]-tabulated_forces.e[table_start + tablepos])/ia_params->TAB_stepsize;
      x0 = ia_params->TAB_minval-tabulated_forces.e[table_start + tablepos]/b;
      return ( (tabulated_energies.e[table_start + tablepos]
		+ 0.5*b*SQR(ia_params->TAB_minval-x0))
	       - 0.5*b*SQR(dist-x0) );
    }

    phi = (dindex - tablepos);
 
    return  tabulated_energies.e[table_start + tablepos]*(1-phi) 
      + tabulated_energies.e[table_start + tablepos+1]*phi;
  }
#endif
  return 0.0;
}


/** check the tabulated forcecap to see that it is sensible \warning
    This routine will probably give strange results if forcecap is
    applied before the table is loaded.
    Needs feature TABULATED compiled in (see \ref config.h). */
MDINLINE void check_tab_forcecap(double force_cap)
{
#ifdef TABULATED
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
#endif
}

/* BONDED INTERACTIONS */


/** Force factor lookup in a force table for bonded interactions (see
    \ref Bonded_ia_parameters). The force is calculated by linear
    interpolation between the closest tabulated values. There is no
    check for the upper bound! 
    Needs feature TABULATED compiled in (see \ref config.h).*/
MDINLINE double bonded_tab_force_lookup(double val, int type_num)
{
#ifdef TABULATED
  int    ind;
  double dind;
  
  dind = (val - bonded_ia_params[type_num].p.tab.minval) 
    * bonded_ia_params[type_num].p.tab.invstepsize;

  if( dind < 0.0 ) ind = 0;
  else ind = (int)dind;

  dind = dind - ind;
  /* linear interpolation between data points */
  return  bonded_ia_params[type_num].p.tab.f[ind]*(1.0-dind) 
    + bonded_ia_params[type_num].p.tab.f[ind+1]*dind;
#else
  return 0.0;
#endif
}

/** Energy lookup in a energy table for bonded interactions (see \ref
    Bonded_ia_parameters). The force is calculated by linear
    interpolation between the closest tabulated values. There is no
    check for the upper bound! 
    Needs feature TABULATED compiled in (see \ref config.h). */
MDINLINE double bonded_tab_energy_lookup(double val, int type_num)
{
#ifdef TABULATED
  int ind;
  double dind;
  
  dind = (val - bonded_ia_params[type_num].p.tab.minval) 
    * bonded_ia_params[type_num].p.tab.invstepsize;

  if( dind < 0.0 ) ind = 0;
  else ind = (int)dind;

  dind = dind - ind;
  /* linear interpolation between data points */
  return  bonded_ia_params[type_num].p.tab.e[ind]*(1.0-dind) 
    + bonded_ia_params[type_num].p.tab.e[ind+1]*dind;
#else
  return 0.0;
#endif
}

/** Calculate a tabulated bond length force with number type_num (see
    \ref Bonded_ia_parameters) between particles p1 and p2 and add it
    to the particle forces. The force acts in the direction of the
    connecting vector between the particles. For distances smaller
    than the tabulated range it uses a linear extrapolation based on
    the first two tabulated force values.
    Needs feature TABULATED compiled in (see \ref config.h). */
MDINLINE void add_tab_bond_force(Particle *p1, Particle *p2, int type_num) 
{
#ifdef TABULATED
  int i;
  double dx[3], dist, fac;


  get_mi_vector(dx, p1->r.p, p2->r.p);
  dist=sqrt(sqrlen(dx));

  if(dist > bonded_ia_params[type_num].p.tab.maxval) {
    fprintf(stderr,"%d: add_tab_bond_force: ERROR: Bond  between Pair (%d,%d) broken: dist=%f\n"
	    ,this_node, p1->p.identity, p2->p.identity, dist); 
    errexit();
  }

  fac = bonded_tab_force_lookup(dist, type_num)/dist;
  
  for(i=0;i<3;i++) {
    p1->f.f[i] -= fac*dx[i];
    p2->f.f[i] += fac*dx[i];
#ifdef NPT
    if(integ_switch == INTEG_METHOD_NPT_ISO)
      nptiso.p_vir[i] -= fac*dx[i] * dx[i];
#endif
  }

  ONEPART_TRACE(if(p1->p.identity==check_id) fprintf(stderr,"%d: OPT: TAB BOND f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p1->f.f[0],p1->f.f[1],p1->f.f[2],p2->p.identity,dist,fac));
  ONEPART_TRACE(if(p2->p.identity==check_id) fprintf(stderr,"%d: OPT: TAB BOND f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p2->f.f[0],p2->f.f[1],p2->f.f[2],p1->p.identity,dist,fac));

#endif
}

/** Calculate and return a tabulated bond length energy with number
    type_num (see \ref Bonded_ia_parameters) between particles p1 and
    p2. For distances smaller than the tabulated range it uses a
    quadratic extrapolation based on the first two tabulated force
    values and the first tabulated energy value. 
    Needs feature TABULATED compiled in (see \ref config.h). */
MDINLINE double tab_bond_energy(Particle *p1, Particle *p2, int type_num) 
{
#ifdef TABULATED
  double dx[3], dist;

  get_mi_vector(dx, p1->r.p, p2->r.p);
  dist=sqrt(sqrlen(dx));

  /* For distances smaller than the tabulated minimim take quadratic
     extrapolation from first two values of the force table! 
     This corresponds to the linear extrapolation of the force at that point.
     This sould not occur too often, since it is quite expensive!
  */
  if( dist <  bonded_ia_params[type_num].p.tab.minval) {
    double x0, b;
    b = (bonded_ia_params[type_num].p.tab.f[1]-bonded_ia_params[type_num].p.tab.f[0])
      *bonded_ia_params[type_num].p.tab.invstepsize;
    x0 = bonded_ia_params[type_num].p.tab.minval-bonded_ia_params[type_num].p.tab.f[0]/b;
    return ( (bonded_ia_params[type_num].p.tab.e[0]
	      +0.5*b*SQR(bonded_ia_params[type_num].p.tab.minval-x0))
	     -0.5*b*SQR(dist-x0) );
  }

  return bonded_tab_energy_lookup(dist, type_num);
#else
  return 0.0;
#endif
}


/** Calculate a tabulated bond angle force with number type_num (see
    \ref Bonded_ia_parameters) between particles p_left, p_mid and
    p_right and add it to the particle forces. The force on p_left and
    p_right acts perpendicular to the connecting vector between the
    particle and p_mid and in the plane defined by the three
    particles. The force on the middle particle balances the other two
    forces. The forces are scaled with the invers length of the
    connecting vectors. It is assumed that the potential is tabulated
    for all angles between 0 and Pi. 
    Needs feature TABULATED compiled in (see \ref config.h). */
MDINLINE void add_tab_angle_force(Particle *p_mid, Particle *p_left, 
				  Particle *p_right, int type_num)
{
#ifdef TABULATED
    double cosine, phi, invsinphi, vec1[3], vec2[3], d1i, d2i, dist2,  fac, f1=0.0, f2=0.0;
  int j;

  /* vector from p_left to p_mid */
  get_mi_vector(vec1, p_mid->r.p, p_left->r.p);
  dist2 = sqrlen(vec1);
  d1i = 1.0 / sqrt(dist2);
  for(j=0;j<3;j++) vec1[j] *= d1i;
  /* vector from p_mid to p_right */
  get_mi_vector(vec2, p_right->r.p, p_mid->r.p);
  dist2 = sqrlen(vec2);
  d2i = 1.0 / sqrt(dist2);
  for(j=0;j<3;j++) vec2[j] *= d2i;
  /* scalar produvt of vec1 and vec2 */
  cosine = scalar(vec1, vec2);
  phi = acos(cosine);
  invsinphi = sin(phi);
  if(invsinphi < TINY_SIN_VALUE) invsinphi = TINY_SIN_VALUE;
  invsinphi = 1.0/invsinphi;
  /* look up force factor */
  fac = bonded_tab_force_lookup(phi, type_num);
  /* apply bend forces */
  for(j=0;j<3;j++) {
    f1               = fac * (cosine * vec1[j] - vec2[j])*invsinphi * d1i;
    f2               = fac * (cosine * vec2[j] - vec1[j])*invsinphi * d2i;
    p_left->f.f[j]  -= f1;
    p_mid->f.f[j]   += (f1-f2);
    p_right->f.f[j] += f2;
  }
#endif

}

/** Calculate and return tabulated bond angle energy with number
    type_num (see \ref Bonded_ia_parameters) between particles p_left,
    p_mid and p_right. It is assumed that the potential is tabulated
    for all angles between 0 and Pi. 
    Needs feature TABULATED compiled in (see \ref config.h). */
MDINLINE double tab_angle_energy(Particle *p_mid, Particle *p_left, 
				 Particle *p_right, int type_num)
{
#ifdef TABULATED
  double phi, vec1[3], vec2[3], vl1, vl2; 

  /* vector from p_mid to p_left */
  get_mi_vector(vec1, p_mid->r.p, p_left->r.p);
  vl1 = sqrt(sqrlen(vec1));
  /* vector from p_right to p_mid */
  get_mi_vector(vec2, p_right->r.p, p_mid->r.p);
  vl2 = sqrt(sqrlen(vec2));
  /* calculate phi */
  phi = acos( scalar(vec1, vec2) / (vl1*vl2) );

  return bonded_tab_energy_lookup(phi, type_num);
#else
  return 0.0;
#endif
}

/** Calculate a tabulated dihedral force with number type_num (see
    \ref Bonded_ia_parameters) between particles p1. p2, p3 and p4 and
    add it to the particle forces. This function is not tested yet.
    Needs feature TABULATED compiled in (see \ref config.h). */
MDINLINE void add_tab_dihedral_force(Particle *p1, Particle *p2, 
				     Particle *p3, Particle *p4, int type_num)
{
#ifdef TABULATED
  int i;
  /* vectors in plane of particle triples (p1,p2,p3) and (p2,p3,p4) and thier length. */
  double vec1[3], vecl1, vec2[3], vecl2;
  /* dihedral angle, cosine of the dihedral angle, cosine of the bond angles */
  double phi, cosphi, cos1223, cos2334;
  /* force factors */
  double fac, f1, f2, f4;

  phi = calc_dihedral_angle(p1, p2, p3, p4, vec1, &vecl1, vec2, &vecl2, 
			    &cos1223, &cos2334, &cosphi);

  fac = bonded_tab_force_lookup(phi, type_num);
  /* apply dihedral forces */
  vecl1 = 1.0/vecl1; vecl2 = 1.0/vecl2;
  for(i=0;i<3;i++) {
    f1 = fac * (vec2[i] - vec1[i]*cosphi) * vecl1;
    f4 = fac * (vec1[i] - vec2[i]*cosphi) * vecl2;
    f2 = (cos1223 - 1.0)*f1 - cos2334*f4;
    p1->f.f[i] += f1;
    p2->f.f[i] += f2;
    p3->f.f[i] -= (f1 + f2 + f4);
    p4->f.f[i] += f4;
  }
#endif
}

/** Calculate and return a tabulated dihedral energy with number
    type_num (see \ref Bonded_ia_parameters) between particles p1. p2,
    p3 and p4. This function is not tested yet. 
    Needs feature TABULATED compiled in (see \ref config.h). */
MDINLINE double tab_dihedral_energy(Particle *p1, Particle *p2, 
				  Particle *p3, Particle *p4, int type_num)
{
#ifdef TABULATED
  /* vectors in plane of particle triples (p1,p2,p3) and (p2,p3,p4) and thier length. */
  double vec1[3], vecl1, vec2[3], vecl2;
  /* dihedral angle, cosine of the dihedral angle, cosine of the bond angles */
  double phi, cosphi, cos1223, cos2334;

  phi = calc_dihedral_angle(p1, p2, p3, p4, vec1, &vecl1, vec2, &vecl2, 
			    &cos1223, &cos2334, &cosphi);

  return bonded_tab_energy_lookup(phi, type_num);
#else
  return 0.0;
#endif
}

#endif
