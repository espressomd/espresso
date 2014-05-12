/*
  Copyright (C) 2010,2012,2013 The ESPResSo project
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
#include "dihedral.hpp"

#ifndef _OVERLAP_H
#define _OVERLAP_H

/** \file overlap.hpp
 *  Routines to calculate the energy and/or force 
 *  for bonds, angles and dihedrals as a sum of N functions in the forms:
 *  bonds --- parameter: [N, a_i, b_i, c_i], function: U(bond)  = sum_(i=1,N) {a_i*exp[-(bond-b_i)^2 /(c_i^2)]}.
 *  angles--- parameter: [N, a_i, b_i, c_i], function: U(cos(angl)) = sum_(i=1,N) {a_i*exp[-(cos(angl)-b_i)^2 / (2 * c_i^2)]}.
 *  dihedrals---parameter: [N, a_i, b_i, c_i], function: U(dihe)= sum_(i=1,N) {a_i*(1+Math.cos(c_i*dihe+b_i))}.
 *  Require feature OVERLAPPED compiled in myconfig.hpp (for more info of FEATURES, see \ref config.hpp ).
*/

#ifdef OVERLAPPED 

/** Bonded overlapped potentials: Reads overlapped parameters from a file.  
    ia_params are then communicated to each node \warning No checking is
    performed for the file read!! */
int overlapped_bonded_set_params(int bond_type, int overlap_type,
				 char * filename);

/** Computes the two body overlapped bonded force.
    Adds this force to the particle forces in forces.hpp (see \ref tclcommand_inter). 
    @param p1        Pointer to first particle.
    @param p2        Pointer to second/middle particle.
    @param iaparams  bond type number of the angle interaction (see \ref tclcommand_inter).
    @param dx        particle distance vector
    @param force     returns force of particle 1
    @return 0.
    Needs feature OVERLAPPED compiled in (see \ref config.hpp).
*/
inline int calc_overlap_bond_force(Particle *p1, Particle *p2, Bonded_ia_parameters *iaparams, double dx[3], double force[3]) 
{
  int i;
  double fac;

  int ig;
  double ev;
  double prob=0.; 
  double Apart=0.;

  double dist2 = sqrlen(dx);
  double dist = sqrt(dist2);

  /* compute fac = - 1/r * dU/dr */
  /* prob = sum_(i=1,N) [ a_i * exp(-(r -b_i)^2/c_i^2))] */
  /* Apart = sum_(i=1,N) [ a_i * exp(-(r -b_i)^2/c_i^2)) * (x-b_i) / c_i^2] */
  /* fac = - 2.0 * ( 1/r + Apart / prob ) */
  for(ig=0; ig<iaparams->p.overlap.noverlaps; ig++) {
        ev = (dist - iaparams->p.overlap.para_b[ig]) / iaparams->p.overlap.para_c[ig];
        prob = prob + iaparams->p.overlap.para_a[ig] * exp (-1.0*ev*ev); 
        Apart = Apart + iaparams->p.overlap.para_a[ig] * exp (-1.0*ev*ev) * (dist-iaparams->p.overlap.para_b[ig]) / (iaparams->p.overlap.para_c[ig] * iaparams->p.overlap.para_c[ig]);
  }
  fac = -2. * ( 1 / dist + Apart / prob);
  fac /= dist;

  /* compute force */
  for(i=0;i<3;i++)
    force[i] = fac*dx[i];

  ONEPART_TRACE(if(p1->p.identity==check_id) fprintf(stderr,"%d: OPT: HARMONIC f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p1->f.f[0],p1->f.f[1],p1->f.f[2],p2->p.identity,dist2,fac));
  ONEPART_TRACE(if(p2->p.identity==check_id) fprintf(stderr,"%d: OPT: HARMONIC f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p2->f.f[0],p2->f.f[1],p2->f.f[2],p1->p.identity,dist2,fac));

  return 0;
}

/** Computes the two body overlapped angle interaction energy (see \ref tclcommand_inter, \ref tclcommand_analyze). 
    @param p1        Pointer to first particle.
    @param p2        Pointer to second/middle particle.
    @param iaparams  bond type number of the angle interaction (see \ref tclcommand_inter).
    @param dx        particle distance vector
    @param _energy   returns energy of this interaction
    @return 0.
    Needs feature OVERLAPPED compiled in (see \ref config.hpp).
*/
inline int overlap_bond_energy(Particle *p1, Particle *p2, Bonded_ia_parameters *iaparams, double dx[3], double *_energy) 
{

  int ig;
  double ev;
  double prob=0.;
  
  double dist2 = sqrlen(dx);
  double dist = sqrt(dist2);

  /*compute overlapped bond energy */
  /*prob = sum_(i=1,N) [ a_i * exp(-(r -b_i)^2/c_i^2))] */
  /*pot  = -1.0 * Log { prob / r^2 } */
  for(ig=0; ig<iaparams->p.overlap.noverlaps; ig++) {
        ev = (dist - iaparams->p.overlap.para_b[ig]) / iaparams->p.overlap.para_c[ig];
        prob = prob + iaparams->p.overlap.para_a[ig] * exp (-1.0*ev*ev);
   }
  *_energy = - log (prob/dist2);

  return 0;
}

/** Computes the three body overlapped angle interaction force.
    Adds this force to the particle forces in forces.hpp (see \ref tclcommand_inter). 
    @param p_mid     Pointer to second/middle particle.
    @param p_left    Pointer to first/left particle.
    @param p_right   Pointer to third/right particle.
    @param iaparams  bond type number of the angle interaction (see \ref tclcommand_inter).
    @param force1 returns force of particle 1
    @param force2 returns force of particle 2
    @return 0
    Needs feature OVERLAPPED compiled in (see \ref config.hpp). 
*/
inline int calc_overlap_angle_force(Particle *p_mid, Particle *p_left, 
				  Particle *p_right, Bonded_ia_parameters *iaparams,
				  double force1[3], double force2[3])
{
  double cosine, vec1[3], vec2[3], d1i, d2i, dist2,  fac, f1=0.0, f2=0.0;
  int j; 
  
  int ig;
  double ev;
  double prob=0.; 
  double Apart=0.;

  cosine=0.0;
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
  /* Notice: cosine = - costheta */
  cosine = scalar(vec1, vec2);

  if ( cosine >  TINY_COS_VALUE ) cosine = TINY_COS_VALUE;
  if ( cosine < -TINY_COS_VALUE)  cosine = -TINY_COS_VALUE;

  /* compute fac = - dU/d(costheta) */
  /* prob = sum_(i=1,N) [ a_i * exp(-(costheta -b_i)^2/c_i^2))] */
  /* Apart = sum_(i=1,N) [ a_i * exp(-(costheta -b_i)^2/c_i^2)) * (costheta-b_i) / c_i^2] */
  /* fac = -2.0 * ( Apart / prob ) */

  for(ig=0; ig<iaparams->p.overlap.noverlaps; ig++) {
        ev = (-cosine - iaparams->p.overlap.para_b[ig]) / iaparams->p.overlap.para_c[ig];
        prob = prob + iaparams->p.overlap.para_a[ig] * exp (-1.0*ev*ev);
        Apart = Apart + iaparams->p.overlap.para_a[ig] * exp (-1.0*ev*ev) * (-cosine-iaparams->p.overlap.para_b[ig]) / (iaparams->p.overlap.para_c[ig] * iaparams->p.overlap.para_c[ig]);
  }
  fac = -2. * ( Apart / prob);

  /* compute force */
  for(j=0;j<3;j++) {
    f1               = fac * (cosine * vec1[j] - vec2[j]) * d1i;
    f2               = fac * (cosine * vec2[j] - vec1[j]) * d2i;

    force1[j] = (f1-f2);
    force2[j] = -f1;
  }
  return 0;
}

/** Computes the three body overlapped angle interaction energy (see \ref tclcommand_inter, \ref tclcommand_analyze). 
    @param p_mid        Pointer to first particle.
    @param p_left        Pointer to second/middle particle.
    @param p_right        Pointer to third particle.
    @param iaparams  bond type number of the angle interaction (see \ref tclcommand_inter).
    @param _energy   return energy pointer.
    @return 0.
    Needs feature OVERLAPPED compiled in (see \ref config.hpp). 
*/
inline int overlap_angle_energy(Particle *p_mid, Particle *p_left, 
			      Particle *p_right, Bonded_ia_parameters *iaparams,
			      double *_energy)
{
  double cosine, vec1[3], vec2[3],  d1i, d2i, dist2;
  int j;

  int ig;
  double ev;
  double prob=0.; 

  cosine=0.0;
  /* vector from p_mid to p_left */
  get_mi_vector(vec1, p_mid->r.p, p_left->r.p);
  dist2 = sqrlen(vec1);
  d1i = 1.0 / sqrt(dist2);
  for(j=0;j<3;j++) vec1[j] *= d1i;
  /* vector from p_right to p_mid */
  get_mi_vector(vec2, p_right->r.p, p_mid->r.p);
  dist2 = sqrlen(vec2);
  d2i = 1.0 / sqrt(dist2);
  for(j=0;j<3;j++) vec2[j] *= d2i;
  /* scalar produvt of vec1 and vec2 */
  /* Notice: cosine = - costheta */
  cosine = scalar(vec1, vec2);

  if ( cosine >  TINY_COS_VALUE)  cosine = TINY_COS_VALUE;
  if ( cosine < -TINY_COS_VALUE)  cosine = -TINY_COS_VALUE;

  /* compute overlapped cosangl energy */
  /*prob = sum_(i=1,N) [ a_i * exp(-(costheta -b_i)^2/c_i^2))] */
  /*pot  = -1.0 * Log { prob } */
  for(ig=0; ig<iaparams->p.overlap.noverlaps; ig++) {
        ev = (-cosine - iaparams->p.overlap.para_b[ig]) / iaparams->p.overlap.para_c[ig];
        prob = prob + iaparams->p.overlap.para_a[ig] * exp (-1.0*ev*ev);
   }
  *_energy = - log (prob);

  return 0;

}

/** Computes the four body overlapped dihedral interaction force.
    Adds this force to the particle forces in forces.hpp (see \ref tclcommand_inter). 
    @param p1, p2, p3, p4 define the angle between the planes p1,p2,p3 and p2,p3,p4
    @param iaparams  bond type number of the angle interaction (see \ref tclcommand_inter).
    @param force1 returns force of particle 1
    @param force2 returns force of particle 2
    @param force3 returns force of particle 3
    @return 0
    Needs feature OVERLAPPED compiled in (see \ref config.hpp). 
*/
inline int calc_overlap_dihedral_force(Particle *p2, Particle *p1,
					 Particle *p3, Particle *p4, Bonded_ia_parameters *iaparams,
					 double force2[3], double force1[3], double force3[3])
{
 int i;
  /* vectors for dihedral angle calculation */
  double v12[3], v23[3], v34[3], v12Xv23[3], v23Xv34[3], l_v12Xv23, l_v23Xv34;
  double v23Xf1[3], v23Xf4[3], v34Xf4[3], v12Xf1[3];
  /* dihedral angle, cosine of the dihedral angle */
  double phi, cosphi, sinmphi_sinphi;
  /* force factors */
  double fac, f1[3], f4[3];
  
  int ig;
  double f_dihe = 0.;

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
                //fac = sum_(i=1,N) { a_i * c_i * c_i * cos(c_i*x + b_i)/cos(phi) }
                //fac = sum_(i=1,N) { a_i * c_i * sin(c_i*phi + b_i) /sin(phi)}
  for(ig=0; ig<iaparams->p.overlap.noverlaps; ig++) {
  	fac = iaparams->p.overlap.para_a[ig] * iaparams->p.overlap.para_c[ig];

  	if(fabs(sin(phi)) < TINY_SIN_VALUE) {
    		/*(comes from taking the first term of the MacLaurin expansion of
      		sin(n*phi - phi0) and sin(phi) and then making the division).
      		The original code had a 2PI term in the cosine (cos(2PI - nPhi))
      		but I removed it because it wasn't doing anything. AnaVV*/
    		sinmphi_sinphi = iaparams->p.overlap.para_c[ig]*
      			cos(iaparams->p.overlap.para_c[ig] * phi + iaparams->p.overlap.para_b[ig])/cosphi;
  	}
  	else {
    		sinmphi_sinphi = sin(iaparams->p.overlap.para_c[ig] * phi - iaparams->p.overlap.para_b[ig])/sin(phi);
 	}
  	fac *= sinmphi_sinphi;
  	f_dihe += fac;
  }

  /* store dihedral forces */
  for(i=0;i<3;i++) {
      force1[i] = f_dihe*v23Xf1[i];
      force2[i] = f_dihe*(v34Xf4[i] - v12Xf1[i] - v23Xf1[i]);
      force3[i] = f_dihe*(v12Xf1[i] - v23Xf4[i] - v34Xf4[i]);
  }
  return 0;
}

/** Computes the four body overlapped dihedral interaction energy (see \ref tclcommand_inter, \ref tclcommand_analyze). 
    @param p1, p2, p3, p4 define the angle between the planes p1,p2,p3 and p2,p3,p4
    @param iaparams  bond type number of the angle interaction (see \ref tclcommand_inter).
    @param _energy   return energy pointer.
    @return 0.
    Needs feature OVERLAPPED compiled in (see \ref config.hpp). 
*/
inline int overlap_dihedral_energy(Particle *p2, Particle *p1, 
				 Particle *p3, Particle *p4, Bonded_ia_parameters *iaparams,
				 double *_energy)
{
 /* vectors for dihedral calculations. */
  double v12[3], v23[3], v34[3], v12Xv23[3], v23Xv34[3], l_v12Xv23, l_v23Xv34;
  /* dihedral angle, cosine of the dihedral angle */
  double phi, cosphi;
  /* energy factors */
  double fac ;

  int ig;
  double pot = 0.;

  calc_dihedral_angle(p1, p2, p3, p4, v12, v23, v34, v12Xv23, &l_v12Xv23, v23Xv34, &l_v23Xv34, &cosphi, &phi);

  /* compute overlapped dihedral energy */
  /* pot (phi) = sum_(i=1,N) { a_i* [1 + cos(int(c_i)*phi + b_i)]} */
  for(ig=0; ig<iaparams->p.overlap.noverlaps; ig++) {
  	fac = cos(iaparams->p.overlap.para_c[ig] * phi + iaparams->p.overlap.para_b[ig]);
  	fac += 1.0;
  	fac *=  iaparams->p.overlap.para_a[ig];
     	pot = pot + fac; 
  }
  *_energy = pot;

  return 0;

}
#endif /* ifdef OVERLAPPED */

#endif /* ifndef OVERLAP_H */
