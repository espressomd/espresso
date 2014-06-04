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
#ifndef _TAB_H
#define _TAB_H

/** \file tab.hpp
 *  Routines to calculate the  energy and/or  force 
 *  for a particle pair or bonds via interpolating from lookup tables.
 *  \ref forces.cpp
 *  Needs feature TABULATED compiled in (see \ref config.hpp).
*/

#include "utils.hpp"
#include "interaction_data.hpp"
#include "particle_data.hpp"
#include "mol_cut.hpp"
#include "dihedral.hpp"

/* should be changed to file containing force caps */
#include "forcecap.hpp"


#ifdef TABULATED

/** Non-Bonded tabulated potentials:
    Reads tabulated parameters and force and energy tables from a file.
    ia_params and force/energy tables are then communicated to each node

    @param part_type_a particle type for which the interaction is defined
    @param part_type_b particle type for which the interaction is defined
    @param filename from which file to fetch the data

    @return <ul>
    <li> 0 on success
    <li> 1 on particle type mismatches
    <li> 2 file name too long
    <li> 3 cannot open the file
    <li> 4 file too short
    <li> 5 file broken, cannot parse numbers
    <li> 6 number of points of existing potential changed
    </ul>
*/
int tabulated_set_params(int part_type_a, int part_type_b, char* filename);

/** Bonded tabulated potentials: Reads tabulated parameters and force
    and energy tables from a file.  ia_params and force/energy tables
    are then communicated to each node.

    @param bond_type bond type for which the interaction is defined
    @param tab_type table type, TAB_BOND_LENGTH, TAB_BOND_ANGLE, TAB_BOND_DIHEDRAL
    @param filename from which file to fetch the data

    @return <ul>
    <li> 0 on success
    <li> 1 if wrong bond type
    <li> 2 currently unused
    <li> 3 cannot open the file
    <li> 4 file too short
    <li> 5 file broken, cannot parse numbers
    <li> 6 parameter out of bounds
    </ul>
*/
int tabulated_bonded_set_params(int bond_type, int tab_type, char * filename);

/** Add a non-bonded pair force by linear interpolation from a table.
    Needs feature TABULATED compiled in (see \ref config.hpp). */
inline void add_tabulated_pair_force(Particle *p1, Particle *p2, IA_parameters *ia_params,
				       double d[3], double dist, double force[3])
{
  if (CUTOFF_CHECK(dist < ia_params->TAB_maxval)){ 
    double phi, dindex, fac;
    int tablepos, table_start,j;
    double rescaled_force_cap = force_cap/dist;
    
    fac = 0.0;

    table_start = ia_params->TAB_startindex;
    dindex = (dist-ia_params->TAB_minval)/ia_params->TAB_stepsize;
    tablepos = (int)(floor(dindex));  

    if ( dist > ia_params->TAB_minval ) {
      phi = dindex - tablepos;	  
      fac = tabulated_forces.e[table_start + tablepos]*(1-phi) + tabulated_forces.e[table_start + tablepos+1]*phi;
    }
    else {
      /* Use an extrapolation beyond the table */
      if ( dist > 0 ) {
	tablepos = 0;
	phi = dindex - tablepos;	  
	fac = (tabulated_forces.e[table_start]*ia_params->TAB_minval)*(1-phi) + 
	  (tabulated_forces.e[table_start+1]*(ia_params->TAB_minval+ia_params->TAB_stepsize))*phi;
	fac = fac/dist;
      }
      else { /* Particles on top of each other .. leave fac as 0.0 */
      }
    }
    
    if ( rescaled_force_cap < fac && force_cap > 0.0) {
      fac = rescaled_force_cap;
    }

    for(j=0;j<3;j++)
      force[j] += fac * d[j];
  }
}

/** Add a non-bonded pair energy by linear interpolation from a table.
    Needs feature TABULATED compiled in (see \ref config.hpp). */
inline double tabulated_pair_energy(Particle *p1, Particle *p2, IA_parameters *ia_params,
				      double d[3], double dist) {
  double phi, dindex;
  int tablepos, table_start;
  double x0, b;
  
  if (CUTOFF_CHECK(dist < ia_params->TAB_maxval)) { 
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
  return 0.0;
}

/** check the tabulated forcecap to see that it is sensible \warning
    This routine will probably give strange results if forcecap is
    applied before the table is loaded.
    Needs feature TABULATED compiled in (see \ref config.hpp). */
void check_tab_forcecap(double force_cap);

/* BONDED INTERACTIONS */

/** Force factor lookup in a force table for bonded interactions (see
    \ref Bonded_ia_parameters). The force is calculated by linear
    interpolation between the closest tabulated values. There is no
    check for the upper bound! 
    Needs feature TABULATED compiled in (see \ref config.hpp).*/
inline double bonded_tab_force_lookup(double val, Bonded_ia_parameters *iaparams)
{
  int    ind;
  double dind;
  
  dind = (val - iaparams->p.tab.minval)*iaparams->p.tab.invstepsize;

  if( dind < 0.0 ) ind = 0;
  else ind = (int)dind;

  dind = dind - ind;
  /* linear interpolation between data points */
  return  iaparams->p.tab.f[ind]*(1.0-dind) + iaparams->p.tab.f[ind+1]*dind;
}

/** Energy lookup in a energy table for bonded interactions (see \ref
    Bonded_ia_parameters). The force is calculated by linear
    interpolation between the closest tabulated values. There is no
    check for the upper bound! 
    Needs feature TABULATED compiled in (see \ref config.hpp). */
inline double bonded_tab_energy_lookup(double val, Bonded_ia_parameters *iaparams)
{
  int ind;
  double dind;
  
  dind = (val - iaparams->p.tab.minval)*iaparams->p.tab.invstepsize;

  if( dind < 0.0 ) ind = 0;
  else ind = (int)dind;

  dind = dind - ind;
  /* linear interpolation between data points */
  return iaparams->p.tab.e[ind]*(1.0-dind) + iaparams->p.tab.e[ind+1]*dind;
}

/** Calculate a tabulated bond length force with number type_num (see
    \ref Bonded_ia_parameters) between particles p1 and p2 and add it
    to the particle forces. The force acts in the direction of the
    connecting vector between the particles. For distances smaller
    than the tabulated range it uses a linear extrapolation based on
    the first two tabulated force values.
    Needs feature TABULATED compiled in (see \ref config.hpp). */
inline int calc_tab_bond_force(Particle *p1, Particle *p2, Bonded_ia_parameters *iaparams, double dx[3], double force[3]) 
{
  int i;
  double fac, dist = sqrt(sqrlen(dx));

  if(dist > iaparams->p.tab.maxval)
    return 1;

  fac = bonded_tab_force_lookup(dist, iaparams);
  
  for(i=0;i<3;i++)
    force[i] = fac*dx[i];

  ONEPART_TRACE(if(p1->p.identity==check_id) fprintf(stderr,"%d: OPT: TAB BOND f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p1->f.f[0],p1->f.f[1],p1->f.f[2],p2->p.identity,dist,fac));
  ONEPART_TRACE(if(p2->p.identity==check_id) fprintf(stderr,"%d: OPT: TAB BOND f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p2->f.f[0],p2->f.f[1],p2->f.f[2],p1->p.identity,dist,fac));

  return 0;
}

/** Calculate and return a tabulated bond length energy with number
    type_num (see \ref Bonded_ia_parameters) between particles p1 and
    p2. For distances smaller than the tabulated range it uses a
    quadratic extrapolation based on the first two tabulated force
    values and the first tabulated energy value. 
    Needs feature TABULATED compiled in (see \ref config.hpp). */
inline int tab_bond_energy(Particle *p1, Particle *p2, Bonded_ia_parameters *iaparams, double dx[3], double *_energy) 
{
  double dist = sqrt(sqrlen(dx));

  if(dist > iaparams->p.tab.maxval)
    return 1;

  /* For distances smaller than the tabulated minimim take quadratic
     extrapolation from first two values of the force table! 
     This corresponds to the linear extrapolation of the force at that point.
     This sould not occur too often, since it is quite expensive!
  */
  if( dist <  iaparams->p.tab.minval) {
    double x0, b;
    b = (iaparams->p.tab.f[1]-iaparams->p.tab.f[0])*iaparams->p.tab.invstepsize;
    x0 = iaparams->p.tab.minval - iaparams->p.tab.f[0]/b;
    *_energy = ( (iaparams->p.tab.e[0] + 0.5*b*SQR(iaparams->p.tab.minval-x0)) -
		 0.5*b*SQR(dist-x0) );
  }
  else
    *_energy = bonded_tab_energy_lookup(dist, iaparams);

  return 0;
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
    Needs feature TABULATED compiled in (see \ref config.hpp). */
inline int calc_tab_angle_force(Particle *p_mid, Particle *p_left, 
				  Particle *p_right, Bonded_ia_parameters *iaparams,
				  double force1[3], double force2[3])
{
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
#ifdef TABANGLEMINUS
  phi = acos(-cosine);
#else
  phi = acos(cosine);
#endif
  invsinphi = sin(phi);
  if(invsinphi < TINY_SIN_VALUE) invsinphi = TINY_SIN_VALUE;
  invsinphi = 1.0/invsinphi;
  /* look up force factor */
  fac = bonded_tab_force_lookup(phi, iaparams);
  /* apply bend forces */
  for(j=0;j<3;j++) {
    f1               = fac * (cosine * vec1[j] - vec2[j])*invsinphi * d1i;
    f2               = fac * (cosine * vec2[j] - vec1[j])*invsinphi * d2i;
    force1[j] = (f1-f2);
    force2[j] = -f1;
  }

  return 0;
}

/* The force on each particle due to a three-body bonded tabulated
   potential is computed. */
inline void calc_angle_3body_tabulated_forces(Particle *p_mid, Particle *p_left,
              Particle *p_right, Bonded_ia_parameters *iaparams,
              double force1[3], double force2[3], double force3[3]) {

  int j;
  double pot_dep;
  double cos_phi;
  double sin_phi;
  double vec31[3];
  double vec21[3];
  double vec12[3]; // espresso convention
  double vec21_sqr;
  double vec31_sqr;
  double vec21_magn;
  double vec31_magn;
  double fj[3];
  double fk[3];
  double phi, dU; // bond angle and d/dphi of U(phi)

  get_mi_vector(vec12, p_mid->r.p, p_left->r.p);
  for(j = 0; j < 3; j++)
    vec21[j] = -vec12[j];

  get_mi_vector(vec31, p_right->r.p, p_mid->r.p);
  vec21_sqr = sqrlen(vec21);
  vec21_magn = sqrt(vec21_sqr);
  vec31_sqr = sqrlen(vec31);
  vec31_magn = sqrt(vec31_sqr);
  cos_phi = scalar(vec21, vec31) / (vec21_magn * vec31_magn);
  sin_phi = sqrt(1.0 - SQR(cos_phi));

  if(cos_phi < -1.0) cos_phi = -TINY_COS_VALUE;
  if(cos_phi >  1.0) cos_phi =  TINY_COS_VALUE;
  phi = acos(cos_phi);

  dU = bonded_tab_force_lookup(phi, iaparams);

  // potential dependent term (dU/dphi * 1 / sin(phi))
  pot_dep = dU / sin_phi;

  for(j = 0; j < 3; j++) {
    fj[j] = vec31[j] / (vec21_magn * vec31_magn) - cos_phi * vec21[j] / vec21_sqr;
    fk[j] = vec21[j] / (vec21_magn * vec31_magn) - cos_phi * vec31[j] / vec31_sqr;
  }

  // note that F1 = -(F2 + F3) in analytical case
  for(j = 0; j < 3; j++) {
    force1[j] = force1[j] - pot_dep * (fj[j] + fk[j]);
    force2[j] = force2[j] + pot_dep * fj[j];
    force3[j] = force3[j] + pot_dep * fk[j];
  }
}


/** Calculate and return tabulated bond angle energy with number
    type_num (see \ref Bonded_ia_parameters) between particles p_left,
    p_mid and p_right. It is assumed that the potential is tabulated
    for all angles between 0 and Pi. 
    Needs feature TABULATED compiled in (see \ref config.hpp). */
inline int tab_angle_energy(Particle *p_mid, Particle *p_left, 
			      Particle *p_right, Bonded_ia_parameters *iaparams,
			      double *_energy)
{
  double phi, vec1[3], vec2[3], vl1, vl2; 

  /* vector from p_mid to p_left */
  get_mi_vector(vec1, p_mid->r.p, p_left->r.p);
  vl1 = sqrt(sqrlen(vec1));
  /* vector from p_right to p_mid */
  get_mi_vector(vec2, p_right->r.p, p_mid->r.p);
  vl2 = sqrt(sqrlen(vec2));
  /* calculate phi */
  phi = acos( scalar(vec1, vec2) / (vl1*vl2) );

  *_energy = bonded_tab_energy_lookup(phi, iaparams);

  return 0;
}

/** Calculate a tabulated dihedral force with number type_num (see
    \ref Bonded_ia_parameters) between particles p1. p2, p3 and p4 and
    add it to the particle forces. This function is not tested yet.
    Needs feature TABULATED compiled in (see \ref config.hpp). */
inline int calc_tab_dihedral_force(Particle *p2, Particle *p1,
				     Particle *p3, Particle *p4, Bonded_ia_parameters *iaparams,
				     double force2[3], double force1[3], double force3[3])
{
  int i;
  /* vectors for dihedral angle calculation */
  double v12[3], v23[3], v34[3], v12Xv23[3], v23Xv34[3], l_v12Xv23, l_v23Xv34;
  double v23Xf1[3], v23Xf4[3], v34Xf4[3], v12Xf1[3];
  /* dihedral angle, cosine of the dihedral angle, cosine of the bond angles */
  double phi, cosphi;
  /* force factors */
  double fac, f1[3], f4[3];

  /* dihedral angle */
  calc_dihedral_angle(p1, p2, p3, p4, v12, v23, v34, v12Xv23, &l_v12Xv23, v23Xv34, &l_v23Xv34, &cosphi, &phi);
  /* dihedral angle not defined - force zero */
  if ( phi == -1.0 ) { 
    for(i=0;i<3;i++) { force1[i] = 0.0; force2[i] = 0.0; force3[i] = 0.0; }
    return 0;
  }

  /* calculate force components (directions) */
  for(i=0;i<3;i++)  {
    f1[i] = (v23Xv34[i] - cosphi*v12Xv23[i])/l_v12Xv23;;
    f4[i] = (v12Xv23[i] - cosphi*v23Xv34[i])/l_v23Xv34;
  }
  vector_product(v23, f1, v23Xf1);
  vector_product(v23, f4, v23Xf4);
  vector_product(v34, f4, v34Xf4);
  vector_product(v12, f1, v12Xf1);

  /* table lookup */
  fac = bonded_tab_force_lookup(phi, iaparams);

  /* store dihedral forces */
  for(i=0;i<3;i++) {
      force1[i] = fac*v23Xf1[i];
      force2[i] = fac*(v34Xf4[i] - v12Xf1[i] - v23Xf1[i]);
      force3[i] = fac*(v12Xf1[i] - v23Xf4[i] - v34Xf4[i]);
  }

  return 0;
}

/** Calculate and return a tabulated dihedral energy with number
    type_num (see \ref Bonded_ia_parameters) between particles p1. p2,
    p3 and p4. This function is not tested yet. 
    Needs feature TABULATED compiled in (see \ref config.hpp). */
inline int tab_dihedral_energy(Particle *p2, Particle *p1, 
				 Particle *p3, Particle *p4, Bonded_ia_parameters *iaparams,
				 double *_energy)
{
  /* vectors for dihedral calculations. */
  double v12[3], v23[3], v34[3], v12Xv23[3], v23Xv34[3], l_v12Xv23, l_v23Xv34;
  /* dihedral angle, cosine of the dihedral angle */
  double phi, cosphi;

  calc_dihedral_angle(p1, p2, p3, p4, v12, v23, v34, v12Xv23, &l_v12Xv23, v23Xv34, &l_v23Xv34, &cosphi, &phi);

  *_energy = bonded_tab_energy_lookup(phi, iaparams);

  return 0;
}
#endif

#endif
