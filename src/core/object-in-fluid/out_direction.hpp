/*
  Copyright (C) 2012,2013,2016 The ESPResSo project
  
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
#ifndef _OBJECT_IN_FLUID_OUT_DIRECTION_H
#define _OBJECT_IN_FLUID_OUT_DIRECTION_H
/** \file out_direction.hpp Routines to calculate the outward direction of the membrane
 *  using a particle quadruple (one particle and its 3 strategically placed neighbors)
*/

#include "utils.hpp"
#include "interaction_data.hpp"
#include "particle_data.hpp"
#include "grid.hpp"

#ifdef MEMBRANE_COLLISION

// set out_direction parameters
int out_direction_set_params(int bond_type);

/** Computes the outward direction of the membrane from one particle and its three neighbors 
 
    (see \ref tclcommand_inter).
    @param p1           Pointer to the central particle.
    @param p2,p3,p4     Pointers to the neighboring particles.

    computes the normal of triangle p2p3p4
    this triangle was initially oriented in such a way that its normal already points out of the object
    normalizes and stores the result as out_direction in p1 data
    @return 0
*/
inline int calc_out_direction(Particle *p1, Particle *p2, Particle *p3, Particle *p4,
				 Bonded_ia_parameters *iaparams)// first-fold-then-the-same approach
{		
    double n[3],dn;
	int j;
	double fp1[3],fp2[3],fp3[3],fp4[3];
	int img[3];
    double AA[3],BB[3],CC[3];

#ifdef GHOST_FLAG
    // first find out which particle out of p1, p2 (possibly p3, p4) is not a ghost particle. In almost all cases it is p2, however, it might be other one. we call this particle reference particle.
    if (p2->l.ghost != 1) {
        //unfold non-ghost particle using image, because for physical particles, the structure p->l.i is correctly set
        memmove(fp2, p2->r.p, 3*sizeof(double));
        memmove(img, p2->l.i, 3*sizeof(int));
        unfold_position(fp2,img);
        // other coordinates are obtained from its relative positions to the reference particle
        get_mi_vector(AA, p1->r.p, fp2);
        get_mi_vector(BB, p3->r.p, fp2);
        get_mi_vector(CC, p4->r.p, fp2);
        for (int i=0; i < 3; i++) { fp1[i] = fp2[i] + AA[i]; fp3[i] = fp2[i] + BB[i]; fp4[i] = fp2[i] + CC[i]; }
    } else {
        // in case  particle p2 is a ghost particle
        if (p1->l.ghost != 1) {
            memmove(fp1, p1->r.p, 3*sizeof(double));
            memmove(img, p1->l.i, 3*sizeof(int));
            unfold_position(fp1,img);
            get_mi_vector(AA, p2->r.p, fp1);
            get_mi_vector(BB, p3->r.p, fp1);
            get_mi_vector(CC, p4->r.p, fp1);
            for (int i=0; i < 3; i++) { fp2[i] = fp1[i] + AA[i]; fp3[i] = fp1[i] + BB[i];  fp4[i] = fp1[i] + CC[i];}
        } else {
            // in case the first and the second particle are ghost particles
            if (p3->l.ghost != 1) {
                memmove(fp3, p3->r.p, 3*sizeof(double));
                memmove(img, p3->l.i, 3*sizeof(int));
                unfold_position(fp3,img);
                get_mi_vector(AA, p1->r.p, fp3);
                get_mi_vector(BB, p2->r.p, fp3);
                get_mi_vector(CC, p4->r.p, fp3);
                for (int i=0; i < 3; i++) { fp1[i] = fp3[i] + AA[i]; fp2[i] = fp3[i] + BB[i]; fp4[i] = fp3[i] + CC[i]; }
            } else {
                // in case the first and the second particle are ghost particles
                if (p4->l.ghost != 1) {
                    memmove(fp4, p4->r.p, 3*sizeof(double));
                    memmove(img, p4->l.i, 3*sizeof(int));
                    unfold_position(fp4,img);
                    get_mi_vector(AA, p1->r.p, fp4);
                    get_mi_vector(BB, p2->r.p, fp4);
                    get_mi_vector(CC, p3->r.p, fp4);
                    for (int i=0; i < 3; i++) { fp1[i] = fp4[i] + AA[i]; fp2[i] = fp4[i] + BB[i]; fp3[i] = fp4[i] + CC[i]; }
                } else {
                    printf("Something wrong in out_direction.hpp: All particles in a bond are ghost particles, impossible to unfold the positions...");
                    return 0;
                }
            }
        }
    }
#endif
#ifndef GHOST_FLAG
    // if ghost flag was not defined we have no other option than to assume the first particle (p1) is a physical one.
    memmove(fp1, p1->r.p, 3*sizeof(double));
    memmove(img, p1->l.i, 3*sizeof(int));
    unfold_position(fp1,img);
    // other coordinates are obtained from its relative positions to the reference particle
    get_mi_vector(AA, p2->r.p, fp1);
    get_mi_vector(BB, p3->r.p, fp1);
    get_mi_vector(CC, p4->r.p, fp1);
    for (int i=0; i < 3; i++) {
        fp2[i] = fp1[i] + AA[i];
        fp3[i] = fp1[i] + BB[i];
        fp4[i] = fp1[i] + CC[i];
    }
#endif
	
	get_n_triangle(fp2,fp3,fp4,n);
	dn=normr(n);
    if ( fabs(dn) < 0.001 )
		printf("out_direction.hpp, calc_out_direction: Length of outward vector is close to zero!\n");
	for(j=0;j<3;j++){
            p1->p.out_direction[j] = n[j]/dn;
    }
    return 0;
}

#endif /* ifdef MEMBRANE_COLLISION */
#endif
