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
#ifndef DIHEDRAL_H
#define DIHEDRAL_H
/** \file dihedral.hpp Routines to calculate the dihedral energy or/and
 *  and force for a particle quadruple.  Note that usage of dihedrals
 *  increases the interaction range of bonded interactions to 2 times
 *  the maximal bond length!  \ref forces.cpp
*/

#include "utils.hpp"
#include "interaction_data.hpp"
#include "grid.hpp"

#define ANGLE_NOT_DEFINED -100

/// set dihedral parameters
int dihedral_set_params(int bond_type, int mult, double bend, double phase);

/** Calculates the dihedral angle between particle quadruple p1, p2,
p3 and p4. The dihedral angle is the angle between the planes
specified by the particle triples (p1,p2,p3) and (p2,p3,p4). 
Vectors a, b and c are the bond vectors between consequtive particles.
If the a,b or b,c are parallel the dihedral angle is not defined in which
case the routine returns phi=-1. Calling functions should check for that
(Written by: Arijit Maitra) */
inline void calc_dihedral_angle(Particle *p1, Particle *p2, Particle *p3, Particle *p4, 
				  double a[3], double b[3], double c[3], 
				  double aXb[3], double *l_aXb, double bXc[3], double *l_bXc, 
				  double *cosphi, double *phi)
{
  int i;

  get_mi_vector(a, p2->r.p, p1->r.p);
  get_mi_vector(b, p3->r.p, p2->r.p);
  get_mi_vector(c, p4->r.p, p3->r.p);

  /* calculate vector product a X b and b X c */
  vector_product(a, b, aXb);
  vector_product(b, c, bXc);

  /* calculate the unit vectors */
  *l_aXb = sqrt(sqrlen(aXb));
  *l_bXc = sqrt(sqrlen(bXc));

  /* catch case of undefined dihedral angle */
   if ( *l_aXb <= TINY_LENGTH_VALUE || *l_bXc <= TINY_LENGTH_VALUE ) { *phi = -1.0; *cosphi = 0; return;}

  for (i=0;i<3;i++) {
    aXb[i] /= *l_aXb;
    bXc[i] /= *l_bXc;
  }

  *cosphi = scalar(aXb, bXc);

  if ( fabs(fabs(*cosphi)-1)  < TINY_SIN_VALUE  ) *cosphi = dround(*cosphi);

  /* Calculate dihedral angle */
  *phi = acos(*cosphi);
  if( scalar(aXb, c) < 0.0 ) *phi = (2.0*PI) - *phi;

}

/** calculate dihedral force between particles p1, p2 p3 and p4 
    Written by Arijit Maitra, adapted to new force interface by Hanjo,
    more general new dihedral form by Ana.
*/
inline int calc_dihedral_force(Particle *p2, Particle *p1, Particle *p3, Particle *p4,
				 Bonded_ia_parameters *iaparams, double force2[3],
				 double force1[2], double force3[2])
{
  int i;
  /* vectors for dihedral angle calculation */
  double v12[3], v23[3], v34[3], v12Xv23[3], v23Xv34[3], l_v12Xv23, l_v23Xv34;
  double v23Xf1[3], v23Xf4[3], v34Xf4[3], v12Xf1[3];
  /* dihedral angle, cosine of the dihedral angle */
  double phi, cosphi, sinmphi_sinphi;
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

  /* calculate force magnitude */
#ifdef OLD_DIHEDRAL
  fac = iaparams->p.dihedral.bend * iaparams->p.dihedral.phase * iaparams->p.dihedral.mult;
#else
  fac = -iaparams->p.dihedral.bend * iaparams->p.dihedral.mult;
#endif

  if(fabs(sin(phi)) < TINY_SIN_VALUE) {
#ifdef OLD_DIHEDRAL
    sinmphi_sinphi = iaparams->p.dihedral.mult * cos(2.0*PI - iaparams->p.dihedral.mult*phi)/cos(phi);
#else
    /*(comes from taking the first term of the MacLaurin expansion of
      sin(n*phi - phi0) and sin(phi) and then making the division).
      The original code had a 2PI term in the cosine (cos(2PI - nPhi))
      but I removed it because it wasn't doing anything. AnaVV*/
    sinmphi_sinphi = iaparams->p.dihedral.mult*
      cos(iaparams->p.dihedral.mult*phi - iaparams->p.dihedral.phase)/cosphi;
#endif
  }
  else {
#ifdef OLD_DIHEDRAL
    sinmphi_sinphi = sin(iaparams->p.dihedral.mult*phi)/sin(phi);
#else
    sinmphi_sinphi = sin(iaparams->p.dihedral.mult*phi - iaparams->p.dihedral.phase)/sin(phi);
#endif
  }

  fac *= sinmphi_sinphi;


  /* store dihedral forces */
  for(i=0;i<3;i++) {
      force1[i] = fac*v23Xf1[i];
      force2[i] = fac*(v34Xf4[i] - v12Xf1[i] - v23Xf1[i]);
      force3[i] = fac*(v12Xf1[i] - v23Xf4[i] - v34Xf4[i]);
  }
  return 0;
}

/** calculate dihedral energy between particles p1, p2 p3 and p4 
    Written by Arijit Maitra, adapted to new force interface by Hanjo */
inline int dihedral_energy(Particle *p1, Particle *p2, Particle *p3, Particle *p4,
			     Bonded_ia_parameters *iaparams, double *_energy) 
{
  /* vectors for dihedral calculations. */
  double v12[3], v23[3], v34[3], v12Xv23[3], v23Xv34[3], l_v12Xv23, l_v23Xv34;
  /* dihedral angle, cosine of the dihedral angle */
  double phi, cosphi;
  /* energy factors */
  double fac;

  calc_dihedral_angle(p1, p2, p3, p4, v12, v23, v34, v12Xv23, &l_v12Xv23, v23Xv34, &l_v23Xv34, &cosphi, &phi);
#ifdef OLD_DIHEDRAL
  fac = iaparams->p.dihedral.phase * cos(iaparams->p.dihedral.mult*phi);
#else
  fac = -cos(iaparams->p.dihedral.mult*phi -iaparams->p.dihedral.phase);
#endif
  fac += 1.0;
  fac *=  iaparams->p.dihedral.bend;

  *_energy = fac;

  return 0;
}

#endif
