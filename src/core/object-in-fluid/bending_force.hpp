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
#ifndef _OBJECT_IN_FLUID_BENDING_FORCE_H
#define _OBJECT_IN_FLUID_BENDING_FORCE_H
/** \file bending_force.hpp Routines to calculate the bending_force energy or/and
 *  and force for a particle quadruple (two triangles that have 2 particles in common)
*/

#include "utils.hpp"
#include "interaction_data.hpp"
#include "particle_data.hpp"
#include "grid.hpp"

/// set bending_force parameters
int bending_force_set_params(int bond_type, double phi0, double kb);

/** Computes the bending force (Dupin2007 eqn. 20 and 21) and adds this
    force to the particle forces (see \ref tclcommand_inter). 
    @param p1,p2,p3     Pointers to particles of triangle 1.
    @param p2,p3,p4     Pointers to particles of triangle 2.
    (triangles have particles p2 and p3 in common)
    @param iaparams  bending stiffness kb, initial rest angle phi0 (see \ref tclcommand_inter).
    @param force1 returns force on particles of triangle 1
    @param force2 returns force on particles of triangle 2
    (p1 += force1; p2 += 0.5*force1+0.5*force2; p3 += 0.5*force1+0.5*force2; p4 += force2;
    @return 0
*/
inline int calc_bending_force(Particle *p2, Particle *p1, Particle *p3, Particle *p4,
				 Bonded_ia_parameters *iaparams, double force1[3],
				 double force2[2])// first-fold-then-the-same approach
{		
	double n1[3],n2[3],dn1,dn2,phi,aa,fac,penal;
	int k;
	double fp1[3],fp2[3],fp3[3],fp4[3];
	int img[3];

	memcpy(fp1, p1->r.p, 3*sizeof(double));
	memcpy(img, p1->l.i, 3*sizeof(int));
	fold_position(fp1, img);
					
	memcpy(fp2, p2->r.p, 3*sizeof(double));
	memcpy(img, p2->l.i, 3*sizeof(int));
	fold_position(fp2, img);

	memcpy(fp3, p3->r.p, 3*sizeof(double));
	memcpy(img, p3->l.i, 3*sizeof(int));
	fold_position(fp3, img);

	memcpy(fp4, p4->r.p, 3*sizeof(double));
	memcpy(img, p4->l.i, 3*sizeof(int));
	fold_position(fp4, img);
	
	get_n_triangle(fp2,fp1,fp3,n1);
	dn1=normr(n1);
	get_n_triangle(fp2,fp3,fp4,n2);
	dn2=normr(n2);
	phi = angle_btw_triangles(fp1,fp2,fp3,fp4);		
	
	if (iaparams->p.bending_force.phi0 < 0.001 || iaparams->p.bending_force.phi0 > 2*M_PI - 0.001) 
		printf("bending_force.h, calc_bending_force: Resting angle is close to zero!!!\n");

	aa = (phi-iaparams->p.bending_force.phi0)/iaparams->p.bending_force.phi0;
	fac = iaparams->p.bending_force.kb * aa;

    penal = (1+1/pow(10*(2*M_PI-phi),2) + 1/pow(10*(phi),2));
	if (penal > 5.) penal = 5.;

//	fac = fac*penal; // This is to penalize the angles smaller than some threshold tr and also it penalizes angles greater than 2*Pi - tr. It prevents the objects to have negative angles.
	if (phi < 0.001 || phi > 2*M_PI - 0.001) printf("bending_force.h, calc_bending_force: Angle approaches 0 or 2*Pi\n");
	
	for(k=0;k<3;k++) {
		force1[k]=fac * n1[k]/dn1;
		force2[k]=fac * n2[k]/dn2;
	}	
  return 0;
}

#endif

