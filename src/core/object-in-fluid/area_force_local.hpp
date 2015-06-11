/*
  Copyright (C) 2012,2013,2014 The ESPResSo project
  
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
	double AA[3],BB[3];
	#ifdef GHOST_FLAG
	// first find out which particle out of p1, p2 (possibly p3, p4) is not a ghost particle. In almost all cases it is p1, however, it might be other one. we call this particle reference particle.
	if (p1->l.ghost != 1) {
		//unfold non-ghost particle using image, because for physical particles, the structure p->l.i is correctly set
		memmove(p11, p1->r.p, 3*sizeof(double));
		memmove(img, p1->l.i, 3*sizeof(int));
		unfold_position(p11,img);
		// other coordinates are obtained from its relative positions to the reference particle
		get_mi_vector(AA, p2->r.p, p11);
		get_mi_vector(BB, p3->r.p, p11);
		for (int i=0; i < 3; i++) { p22[i] = p11[i] + AA[i]; p33[i] = p11[i] + BB[i]; }
	} else {
		// in case the first particle is a ghost particle
		if (p2->l.ghost != 1) {
			memmove(p22, p2->r.p, 3*sizeof(double));
			memmove(img, p2->l.i, 3*sizeof(int));
			unfold_position(p22,img);
			get_mi_vector(AA, p1->r.p, p22);
			get_mi_vector(BB, p3->r.p, p22);
			for (int i=0; i < 3; i++) { p11[i] = p22[i] + AA[i]; p33[i] = p22[i] + BB[i]; }
		} else {
			// in case the first and the second particle are ghost particles
			if (p3->l.ghost != 1) {
				memmove(p33, p3->r.p, 3*sizeof(double));
				memmove(img, p3->l.i, 3*sizeof(int));
				unfold_position(p33,img);
				get_mi_vector(AA, p1->r.p, p33);
				get_mi_vector(BB, p2->r.p, p33);
				for (int i=0; i < 3; i++) { p11[i] = p33[i] + AA[i]; p22[i] = p33[i] + BB[i]; }
			} else {
				printf("Something wrong in area_force_local.hpp: All particles in a bond are ghost particles, impossible to unfold the positions...");
				return 0;
			}
		}
	}
	#endif
	#ifndef GHOST_FLAG
		// if ghost flag was not defined we have no other option than to assume the first particle is a physical one.
		memmove(p11, p1->r.p, 3*sizeof(double));
		memmove(img, p1->l.i, 3*sizeof(int));
		unfold_position(p11,img);
		// other coordinates are obtained from its relative positions to the reference particle
		get_mi_vector(AA, p2->r.p, p11);
		get_mi_vector(BB, p3->r.p, p11);
		for (int i=0; i < 3; i++) { p22[i] = p11[i] + AA[i]; p33[i] = p11[i] + BB[i]; }
	#endif

	for(k=0;k<3;k++) h[k]=1.0/3.0 *(p11[k]+p22[k]+p33[k]);
	//volume+=A * -n[2]/dn * h[2];
	A=area_triangle(p11,p22,p33);

	//aa=(A - iaparams->p.area_force_local.A0_l)/iaparams->p.area_force_local.A0_l;
	//corrected with square root of the triangle's area, see "R.Tothova: Comparison of Different Formulas for Local Area Conservation Modulus in Spring Network Models"
	aa=(A - iaparams->p.area_force_local.A0_l)/sqrt(iaparams->p.area_force_local.A0_l);
		
	//aminusb(3,h,p11,rh);				// area_forces for each triangle node
	vecsub(h,p11,rh);				// area_forces for each triangle node
	hn=normr(rh);
	for(k=0;k<3;k++) {
		force1[k] =  iaparams->p.area_force_local.ka_l * aa * rh[k]/hn;
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
#endif

