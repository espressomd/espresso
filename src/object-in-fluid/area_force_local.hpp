/*
  Copyright (C) 2012,2013 The ESPResSo project
  
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
#ifndef _OBJECT_IN_FLUID_AREA_FORCE_LOCAL_H
#define _OBJECT_IN_FLUID_AREA_FORCE_LOCAL_H
/** \file area_force_local.hpp
 *  Routines to calculate the AREA_FORCE_LOCAL energy or/and and force 
 *  for a particle triple (triangle from mesh). (Dupin2007)
 *  \ref forces.cpp
*/

#include "utils.hpp"
#include "interaction_data.hpp"
#include "particle_data.hpp"
#include "grid.hpp"

/************************************************************/

/** set parameters for the AREA_FORCE_LOCAL potential. 
*/
int area_force_local_set_params(int bond_type, double A0_l, double ka_l);

/** Computes the local area force (Dupin2007 eqn. 15) and adds this
    force to the particle forces (see \ref tclcommand_inter). 
    @param p1,p2,p3     Pointers to triangle particles.
    @param iaparams  elastic area modulus ka, initial area A0 (see \ref tclcommand_inter).
    @param force1 returns force of particle 1
    @param force2 returns force of particle 2
    @param force3 returns force of particle 3
    @return 0
*/
inline int calc_area_force_local(Particle *p1, Particle *p2, Particle *p3,
			      Bonded_ia_parameters *iaparams, double force1[3], double force2[3], double force3[3]) //first-fold-then-the-same approach
{
	int k;	
	double A, aa, h[3], rh[3], hn;
	double p11[3],p22[3],p33[3];
	int img[3];
	
	memcpy(p11, p1->r.p, 3*sizeof(double));
	memcpy(img, p1->l.i, 3*sizeof(int));
	fold_position(p11, img);
					
	memcpy(p22, p2->r.p, 3*sizeof(double));
	memcpy(img, p2->l.i, 3*sizeof(int));
	fold_position(p22, img);

	memcpy(p33, p3->r.p, 3*sizeof(double));
	memcpy(img, p3->l.i, 3*sizeof(int));
	fold_position(p33, img);

	for(k=0;k<3;k++) h[k]=1.0/3.0 *(p11[k]+p22[k]+p33[k]);
	//volume+=A * -n[2]/dn * h[2];
	A=area_triangle(p11,p22,p33);
	aa=(A - iaparams->p.area_force_local.A0_l)/iaparams->p.area_force_local.A0_l;
	
	//aminusb(3,h,p11,rh);				// area_forces for each triangle node
	vecsub(h,p11,rh);				// area_forces for each triangle node
	hn=normr(rh);
	for(k=0;k<3;k++) {
		force1[k] =  iaparams->p.area_force_local.ka_l * aa * rh[k]/hn;
		//(&part1)->f.f[k]+=force[k];
	}
	//aminusb(3,h,p22,rh);				// area_forces for each triangle node
	vecsub(h,p22,rh);				// area_forces for each triangle node
	hn=normr(rh);
	for(k=0;k<3;k++) {
		force2[k] =  iaparams->p.area_force_local.ka_l * aa * rh[k]/hn;
		//(&part2)->f.f[k]+=force[k];
	}
	//aminusb(3,h,p33,rh);				// area_forces for each triangle node
	vecsub(h,p33,rh);				// area_forces for each triangle node
	hn=normr(rh);
	for(k=0;k<3;k++) {
		force3[k] =  iaparams->p.area_force_local.ka_l * aa * rh[k]/hn;
		//(&part3)->f.f[k]+=force[k];
	}
  return 0;
}


inline int calc_area_force_local_complicated(Particle *p1, Particle *p2, Particle *p3,
         Bonded_ia_parameters *iaparams, double force1[3], double force2[3], double force3[3]) // more complicated approach
{
 int k; 
 double A, aa, uh[3], vh[3], rh[3], hn;
 double uu[3],vv[3];

 get_mi_vector(uu, p2->r.p, p1->r.p);
 get_mi_vector(vv, p3->r.p, p1->r.p);

 for(k=0;k<3;k++) rh[k]=1.0/3.0 *(uu[k]+vv[k]);
 for(k=0;k<3;k++) uh[k]=rh[k]-uu[k];
 for(k=0;k<3;k++) vh[k]=rh[k]-vv[k];

 A=area_triangle_new(uu,vv);
 //printf("area %e uu %e %e %e vv %e %e %e\n", A, uu[0], uu[1], uu[2], vv[0], vv[1], vv[2]);
 aa=(A - iaparams->p.area_force_local.A0_l)/iaparams->p.area_force_local.A0_l;

 hn=normr(rh);
 for(k=0;k<3;k++) {
  force1[k] =  iaparams->p.area_force_local.ka_l * aa * rh[k]/hn;
 }
 hn=normr(uh);
 for(k=0;k<3;k++) {
  force2[k] =  iaparams->p.area_force_local.ka_l * aa * uh[k]/hn;
 }
 hn=normr(vh);
 for(k=0;k<3;k++) {
  force3[k] =  iaparams->p.area_force_local.ka_l * aa * vh[k]/hn;
 }
  return 0;
}


#endif

