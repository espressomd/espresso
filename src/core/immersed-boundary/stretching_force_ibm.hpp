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

#ifndef STRETCHING_FORCE_IBM_H
#define STRETCHING_FORCE_IBM_H

#include "utils.hpp"
#include "interaction_data.hpp"
#include "particle_data.hpp"
#include "grid.hpp"
#include "lb.hpp"
#include "integrate.hpp"

extern int stretching_force_law_ibm;

#ifdef STRETCHING_FORCE_IMMERSED_BOUNDARY

const int NEO_HOOK = 0;
const int SKALAK = 1;

/*@}*/
  /** Some general remarks:
   * This file implements membrane strain and area dilation force calculation as described
   * in this paper "Computer simulation study of collective phenomena in dense suspensions   * of red blood cells under shear* by Kruger 2011, pg 145
   */

/** Setting the reference position of a triangular interaction between 3 particles
 *
 *  @param bond_type type of bond refer to interaction_data.hpp
 *  @param ind1 particle id of first particle in the bond
 *  @param ind2 particle id of second particle in the bond
 *  @param ind3 particle id of third particle in the bond
 *  @param max maximum stretch distance between any 2 particles in the bond
 *  @param ks elastic shear modulus as described in literature refer kruger 2011
 *  @param ka area dilation modulus as described in literature refer kruger 2011
 */
int stretching_force_ibm_set_params(int bond_type, int ind1, int ind2, int ind3, double max, double ks, double ka);

/** Resetting the reference position of a triangular interaction between 3 particles
 *
 *  @param bond_type type of bond refer to interaction_data.hpp
 *  @param lo distance between particle 3 and particle 1 in the interaction
 *  @param lp distance between particle 2 and particle 1 in the interaction
 *  @param cospo cosine angle function between lo and lp 
 *  @param sinpo sine angle function between lo and lp 
 *  @param Area0 initial area between the 3 particles in the interaction 
 *  @param max maximum stretch distance between any 2 particles in the bond
 *  @param ks elastic shear modulus as described in literature refer kruger 2011
 *  @param ka area dilation modulus as described in literature refer kruger 2011
 */
int stretching_force_ibm_reset_params(int bond_type, double lo, double lpo, double cospo, double sinpo, double Area0, double max, double ks, double ka);



/** Rotating calculated trielastic forces in the 2d plane back to the 3d plane
 *Use knowledge that the x-axis in rotates system is parallel to r(p1->p2) in original system;
 *To find the corresponding unit vector to y in the rotated system, construct vector perpendicular to r(p1->p2);
 * note that f3 is not calculated here but is implicitly calculated in forces.hpp by f3 = -(f1+f2) which is consistent with the literature
 *  @param f1_rot forces in xy on particle 1
 *  @param f2_rot forces in xy on particle 2
 *  @param f1 forces in xyz on particle 1
 *  @param f2 forces in xyz on particle 2
 *  @param r1 distance between particle 2 and particle 1 in the interaction
 *  @param r2 distance between particle 3 and particle 1 in the interaction
 */
inline void RotateForces(double f1_rot[2], double f2_rot[2], double f1[3], double f2[3], double r1[3], double r2[3]) {
  double xu[3] = {0., 0., 0.};
  double y[3] = {0., 0., 0.};
  double yu[3] = {0., 0., 0.};
  double sca;
  int i;
  unit_vector(r1,xu);
  sca = scalar(r2,xu);
  for(i=0; i<3; i++) {
    y[i] = r2[i] - (sca * xu[i]);
  }
	
  unit_vector(y,yu);
   
  f1[0] = (f1_rot[0] * xu[0]) + (f1_rot[1] * yu[0]);
  f1[1] = (f1_rot[0] * xu[1]) + (f1_rot[1] * yu[1]); 
  f1[2] = (f1_rot[0] * xu[2]) + (f1_rot[1] * yu[2]); 
	
  f2[0] = (f2_rot[0] * xu[0]) + (f2_rot[1] * yu[0]); 
  f2[1] = (f2_rot[0] * xu[1]) + (f2_rot[1] * yu[1]);
  f2[2] = (f2_rot[0] * xu[2]) + (f2_rot[1] * yu[2]);
    
}

/** Main method for calculation of area dilation forces described in literature
 *  @param p_ind1 particle object for particle 1
 *  @param p_ind2 particle object for particle 2
 *  @param p_ind3 particle object for particle 3
 *  @param iaparams parameters of the interaction such as max stretch and ks,ka
 *  @param force1 forces in xy on particle 1
 *  @param force2 forces in xy on particle 2
 */
inline int calc_stretching_force_ibm(Particle *p_ind1, Particle *p_ind2, Particle *p_ind3,
			      Bonded_ia_parameters *iaparams, double force1[3], double force2[3]) 
{
  double dxy, dxx, dyy, dyx; //Displacment gradient tensor matrix elements
  double gxy, gyx, gxx, gyy; //matrix calculations to get I1 and I2
  double e1, e2, i1, i2, i11, i12, i21, i22, i23, i24; // eigen values 1, 2 and D matrix elements  
    double A0, a1, a2, a3, b1, b2, b3;
    double l, lp, sinp, cosp;  // l, l' and cross and dot product between them 
    double vec1[3] = {0., 0., 0.};
    double vec2[3] = {0., 0., 0.}; 
    double vecpro[3] = {0., 0., 0.};
    double f1_rot[2] = {0., 0.};
    double f2_rot[2] = {0., 0.};
    
    //Calculate the current shape of the triangle (l,lp,cos(phi),sin(phi));
    //l = length between 1 and 3

    get_mi_vector(vec2, p_ind3->r.p, p_ind1->r.p);
    // vecsub(p_ind3->r.p,p_ind1->r.p,vec2);

    l = sqrt (sqrlen(vec2));
    // lp = lenght between 1 and 2
    get_mi_vector(vec1, p_ind2->r.p, p_ind1->r.p);
    // vecsub(p_ind2->r.p,p_ind1->r.p,vec1);

    lp = sqrt (sqrlen(vec1));
    //cosp / sinp angle functions between these vectors; calculated directly via the producs
    cosp = scalar(vec1,vec2)/(lp*l);
    vector_product(vec1, vec2, vecpro);
    sinp = sqrt(sqrlen(vecpro))/(l*lp);
    
    if( (lp-iaparams->p.stretching_force_ibm.lpo > iaparams->p.stretching_force_ibm.maxdist) ||  (l-iaparams->p.stretching_force_ibm.lo > iaparams->p.stretching_force_ibm.maxdist)) {
      // triel_params[0] = 1;
      return 1;
    }
    
    
    //Calculate forces in common plane (after assumed triangle rotation/translation in xy-plane);
    //Note that certain geometries and parameters (e.g. a3=0) can be used to speed the code up.
    //For now it will be played safe and done in detail. 
	dxx = lp/iaparams->p.stretching_force_ibm.lpo;
	dxy = ((l*cosp/iaparams->p.stretching_force_ibm.lo) - (lp*iaparams->p.stretching_force_ibm.cospo/iaparams->p.stretching_force_ibm.lpo)) / iaparams->p.stretching_force_ibm.sinpo;
	dyx = 0.0;
	dyy = (l*sinp)/(iaparams->p.stretching_force_ibm.lo * iaparams->p.stretching_force_ibm.sinpo);
	gxx = SQR(dxx)+SQR(dyx);
	gxy = dxx*dxy + dyx*dyy;
	gyx = dxx*dxy + dyy*dyx;
	gyy = SQR(dxy) + SQR(dyy);
	i1 = (gxx + gyy) - 2;
	i2 = ((gxx * gyy) - (gxy * gyx)) - 1;
	i11 = 1.0; i12 = 1.0;
	i21 = gyy; i22 = -gyx; i23 = i22; i24 = gxx;
	
	if(stretching_force_law_ibm == NEO_HOOK) { 
	  e1 = iaparams->p.stretching_force_ibm.ks/6.0;
	  e2 = (-1)*iaparams->p.stretching_force_ibm.ks/(6.0*(i2+1.0)*(i2+1.0));
	}
	  else if(stretching_force_law_ibm == SKALAK) { 
	  //Skalak-Law = Standard
	  e1 = iaparams->p.stretching_force_ibm.ks*(i1+1)/6.0;
	  e2 = (-1)*iaparams->p.stretching_force_ibm.ks/6.0 + iaparams->p.stretching_force_ibm.ka*i2/6.0;
	}
    
    //For sake of better readability shorten the call for the triangle's constants:
	A0 = iaparams->p.stretching_force_ibm.Area0; 
	a1 = iaparams->p.stretching_force_ibm.a1; a2 = iaparams->p.stretching_force_ibm.a2; a3 = iaparams->p.stretching_force_ibm.a3;
	b1 = iaparams->p.stretching_force_ibm.b1; b2 = iaparams->p.stretching_force_ibm.b2; b3 = iaparams->p.stretching_force_ibm.b3;
	
    f1_rot[0] = A0*((-1)*e1*((i11*2*a1*dxx)+(i12*2*b1*dxy))+ (-1)*e2*((i21*2*a1*dxx)+(i22*(a1*dxy+b1*dxx))+(i23*(a1*dxy+b1*dxx))+(i24*2*b1*dxy)));
    f1_rot[1] = A0*((-1)*e1*((i11*0.0)+(i12*2*b1*dyy))+ (-1)*e2*((i21*0.0)+(i22*a1*dyy)+(i23*a1*dyy)+(i24*2*b1*dyy)));

    f2_rot[0] = A0*((-1)*e1*((i11*2*a2*dxx)+(i12*2*b2*dxy))+ (-1)*e2*((i21*2*a2*dxx)+(i22*(a2*dxy+b2*dxx))+(i23*(a2*dxy+b2*dxx))+(
i24*2*b2*dxy)));
    f2_rot[1] = A0*((-1)*e1*((i11*0.0)+(i12*2*b2*dyy))+ (-1)*e2*((i21*0.0)+(i22*a2*dyy)+(i23*a2*dyy)+(i24*2*b2*dyy)));

    //Rotate forces back into original position of triangle
    RotateForces(f1_rot,f2_rot, force1,force2, vec1, vec2);

    return 0;
}

#endif
#endif
