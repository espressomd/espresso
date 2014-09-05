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

#ifndef BENDING_FORCE_IBM_H
#define BENDING_FORCE_IBM_H

#include "utils.hpp"
#include "interaction_data.hpp"
#include "particle_data.hpp"
#include "grid.hpp"

#define TOANGLE 360.0/(2.0*PI)

#ifdef BENDING_FORCE_IMMERSED_BOUNDARY

/*@}*/
  /** Some general remarks:
   * This file implements membrane bending stiffness as described
   * in this paper "Computer simulation study of collective phenomena in dense suspensions   * of red blood cells under shear* by Kruger 2011, pg 149
   */

/** Setting the reference configuration of bending interaction between 4 particles
 *
 *  @param bond_type type of bond refer to interaction_data.hpp
 *  @param ind1 particle id of first particle in the bond
 *  @param ind2 particle id of second particle in the bond
 *  @param ind3 particle id of third particle in the bond
 *  @param ind4 particle id of fourth particle in the bond
 *  @param max maximum stretch distance
 *  @param boo direction of cross product between surface normals
 *  @param kb bending modulus as described in literature refer kruger 2011
 */
int bending_force_ibm_set_params(int bond_type, int ind1, int ind2, int ind3, int ind4, int boo, double max, double kb);

/** Resetting the reference configuration of bending interaction between 4 particles
 *
 *  @param bond_type type of bond refer to interaction_data.hpp
 *  @param theta0 angle between surface normals of the triangles composed of the 4 particles aligned to each other
 *  @param max maximum stretch distance
 *  @param bood direction of cross product between surface normals
 *  @param kb bending modulus as described in literature refer kruger 2011
 */
int bending_force_ibm_reset_params(int bond_type, double bood, double theta0, double kb, double max);

/** Main method for calculation of bending forces described in literature
 *  @param p_ind1 particle object for particle 1
 *  @param p_ind2 particle object for particle 2
 *  @param p_ind3 particle object for particle 3
 *  @param p_ind3 particle object for particle 4
 *  @param iaparams parameters of the interaction such as max stretch, kb and theta0
 *  @param force forces in xyz on particle 1
 *  @param force2 forces in xyz on particle 2
 *  @param force3 forces in xyz on particle 3
 *  @param force4 forces in xyz on particle 4
 */
inline int calc_bending_force_ibm(Particle *p1, Particle *p2, Particle *p3, Particle *p4, Bonded_ia_parameters *iaparams, double force[3], double force2[3], double force3[3], double force4[3]) { 
  double theta, Ai, Aj;
  double dx1[3], dx2[3], dx3[3], n1[3], n2[3];
  double Pre, sc, len;
  double v1l[3], v2l[3], v1[3], v2[3], tmp[3], tmp2[3], term1[3], term2[3];
  double direc[3];
  double desc, DTh;
  int i;
  
  
  //Get vectors making up the two triangles
  get_mi_vector(dx1, p1->r.p, p3->r.p);
  get_mi_vector(dx2, p2->r.p, p3->r.p);
  get_mi_vector(dx3, p4->r.p, p3->r.p);
    
  //printf("dx1: %lf %lf %lf\n", dx1[0], dx1[1], dx1[2]);
  //printf("dx2: %lf %lf %lf\n", dx2[0], dx2[1], dx2[2]);
  //printf("dx3: %lf %lf %lf\n", dx3[0], dx3[1], dx3[2]);
  
  //Get normals on triangle; pointing outwards by definition of indices sequence
  vector_product(dx1, dx2, n1);
  vector_product(dx1, dx3, n2);
    
  if(iaparams->p.bending_force_ibm.boo == 0) {
      n2[0]=-1*n2[0]; n2[1]=-1*n2[1]; n2[2]=-1*n2[2];
  }
   
  //Get 2*area of triangles out of the magnitude of the resulting normals and make the latter unity 
  Ai = sqrt(n1[0]*n1[0] + n1[1]*n1[1] + n1[2]*n1[2]);
  n1[0] = n1[0]/Ai; n1[1]=n1[1]/Ai; n1[2]=n1[2]/Ai; 
  
  Aj = sqrt(n2[0]*n2[0] + n2[1]*n2[1] + n2[2]*n2[2]);
  n2[0] = n2[0]/Aj; n2[1]=n2[1]/Aj; n2[2]=n2[2]/Aj; 
  
  //printf("n1 = %lf %lf %lf\n", n1[0], n1[1], n1[2]);
  //printf("n2 = %lf %lf %lf\n", n2[0], n2[1], n2[2]);
  //printf("Ai = %lf   Aj = %lf\n", Ai, Aj);
  
  //Get the prefactor for the force term
  sc = scalar(n1,n2);
  
  if(sc>1.0) {
      sc = 1.0;
  }
  
  //Get theta as angle between normals
  theta = acos(sc);
  
  //printf("Theta_Pure = %lf\n", theta);
  //printf("n1: %lf %lf %lf\n", n1[0], n1[1], n1[2]);
  //printf("n2: %lf %lf %lf\n", n2[0], n2[1], n2[2]);
  
  vector_product(n1,n2,direc);
  
  //printf("direc = %lf %lf %lf\n", direc[0], direc[1], direc[2]);
  
  desc = scalar(dx1,direc);
  
  //printf("dec = %lf\n", desc);
  
  if(desc<0) {
      theta = -1.0*theta;
  }
  
  //printf("theta = %lf  theta0 = %lf\n", theta, iaparams->p.bending_force_ibm.theta0);
  
  
  DTh = theta-iaparams->p.bending_force_ibm.theta0;  
  
 
  //printf("DTh = %lf\n", DTh);
  
  if(isnan(theta)) {
      printf("Triangle-Pair: %d %d %d %d\n", p1->p.identity, p2->p.identity, p3->p.identity, p4->p.identity);
      printf("n1: %lf %lf %lf\n", n1[0], n1[1], n1[2]);
      printf("n2: %lf %lf %lf\n", n2[0], n2[1], n2[2]);
      printf("scalar: %lf\n", scalar(n1,n2));
      printf("p1: %lf %lf %lf\n", p1->r.p[0], p1->r.p[1], p1->r.p[2]);
      printf("p2: %lf %lf %lf\n", p2->r.p[0], p2->r.p[1], p2->r.p[2]);
      printf("p3: %lf %lf %lf\n", p3->r.p[0], p3->r.p[1], p3->r.p[2]);
      printf("p4: %lf %lf %lf\n", p4->r.p[0], p4->r.p[1], p4->r.p[2]);
  }
  
  
  //printf("%lf %lf\n", theta*TOANGLE, iaparams->p.bending_force_ibm.theta0*TOANGLE);
  
  //printf("sc = %lf\n", sc);
  
  
  if(theta>0) {
    Pre = 1.0*iaparams->p.bending_force_ibm.kb * sin(DTh);  
  } else {
    Pre = -1.0*iaparams->p.bending_force_ibm.kb * sin(DTh); 
  }
  
  //printf("Pre = %lf  kb=%lf\n", Pre, iaparams->p.bending_force_ibm.kb); 
  
  for(i=0; i<3; i++) {
      v1l[i] = n2[i]-sc*n1[i];
      v2l[i] = n1[i]-sc*n2[i];
  }
  
  //printf("v1l: %lf %lf %lf\n", v1l[0], v1l[1], v1l[2]);
  //printf("v2l: %lf %lf %lf\n", v2l[0], v2l[1], v2l[2]);
  
  len = sqrt(sqrlen(v1l));
  
  //printf("len1 = %lf\n", len);
  
  if(len>0) {
      for(i=0;i <3; i++)
	v1[i]=v1l[i]/len;
  }
  
  len = sqrt(sqrlen(v2l));
  
   //printf("len2 = %lf\n", len);
  
  if(len>0) {
      for(i=0;i <3; i++)
	v2[i]=v2l[i]/len;
  }
  
  //printf("v1: %lf %lf %lf\n", v1[0], v1[1], v1[2]);
  //printf("v2: %lf %lf %lf\n", v2[0], v2[1], v2[2]);
  
  //Force for particle 1:
  get_mi_vector(tmp,p2->r.p,p3->r.p); get_mi_vector(tmp2, p3->r.p, p4->r.p);
  vector_product(tmp,v1, term1); vector_product(tmp2,v2, term2);
  
  //printf("tmp: %lf %lf %lf\n", tmp[0], tmp[1], tmp[2]);
  //printf("tmp2: %lf %lf %lf\n", tmp2[0], tmp2[1], tmp2[2]);
  
  //printf("p1f: ");
  
  for(i=0;i<3;i++) {
      force[i] = Pre*(term1[i]/Ai + term2[i]/Aj);
      //printf("%lf ", Pre*(term1[i]/Ai + term2[i]/Aj));
  }
  
  //Force for particle 2:
  get_mi_vector(tmp,p3->r.p,p1->r.p);
  vector_product(tmp,v1, term1);
  
  //printf("\np2f: ");
  
  for(i=0;i<3;i++) {
      force2[i] = Pre*(term1[i]/Ai);
      //printf("%lf ", term1[i]/Ai);
  }
  
  //Force for Particle 3:
  get_mi_vector(tmp,p1->r.p,p2->r.p); get_mi_vector(tmp2, p4->r.p, p1->r.p);
  vector_product(tmp,v1, term1); vector_product(tmp2,v2, term2);
  
  //printf("\np3f: ");
  
  for(i=0;i<3;i++) {
      force3[i] = Pre*(term1[i]/Ai + term2[i]/Aj);
      //printf("%lf ", Pre*(term1[i]/Ai + term2[i]/Aj));
  }
  
  //Force for Particle 4:
  get_mi_vector(tmp,p1->r.p,p3->r.p);
  vector_product(tmp,v2, term1);
  
  //printf("\np4f: ");
  
  for(i=0;i<3;i++) {
      force4[i] = Pre*(term1[i]/Aj);
      //printf("%lf ", Pre*(term1[i]/Aj));
  }
  
  //printf("\n");
  
  return 0;
  
}

#endif
#endif
