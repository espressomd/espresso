/*
  Copyright (C) 2012,2013,2014,2015,2016 The ESPResSo project
  
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
#ifndef _OBJECT_IN_FLUID_OIF_GLOBAL_FORCES_H
#define _OBJECT_IN_FLUID_OIF_GLOBAL_FORCES_H
/** \file oif_global_forces.hpp
 *  Routines to calculate the OIF_GLOBAL_FORCES energy or/and and force 
 *  for a particle triple (triangle from mesh). (Dupin2007)
 *  \ref forces.cpp
*/

#include "utils.hpp"
#include "interaction_data.hpp"
#include "particle_data.hpp"
#include "cells.hpp"
#include "lb.hpp"
#include "grid.hpp"
#include "errorhandling.hpp"

/** set parameters for the OIF_GLOBAL_FORCES potential. 
*/
int oif_global_forces_set_params(int bond_type, double A0_g, double ka_g, double V0, double kv);

/************************************************************/

/** called in force_calc() from within forces.cpp
 *  calculates the global area and global volume for a cell before the forces are handled
 *  sums up parts for area with mpi_reduce from local triangles
 *  synchronization with allreduce
 *
 *  !!! loop over particles from domain_decomposition !!!
 */  

inline void calc_oif_global(double *area_volume, int molType){ //first-fold-then-the-same approach
	double partArea=0.0;
    double part_area_volume[2]; //added

// z volume
	double VOL_partVol=0.,VOL_A,VOL_norm[3],VOL_dn,VOL_hz;


	/** loop over particles */
	int c, np, i ,j;
	Cell *cell;
	Particle *p, *p1, *p2, *p3;
	double p11[3],p22[3],p33[3];
	int img[3];
	double AA[3],BB[3];
	Bonded_ia_parameters *iaparams;
    int type_num, n_partners,id;
    BondedInteraction type;

	int test=0;

	/* Loop local cells */
	for (c = 0; c < local_cells.n; c++) {
		cell = local_cells.cell[c];
		p   = cell->part;
		np  = cell->n;
		/* Loop cell particles */
		for(i=0; i < np; i++) {				
			j = 0;
			p1 = &p[i];
			while(j<p1->bl.n){
				/* bond type */
				type_num = p1->bl.e[j++];
				iaparams = &bonded_ia_params[type_num];			
				type = iaparams->type;
				n_partners = iaparams->num;
				id=p1->p.mol_id;
				if(type == BONDED_IA_OIF_GLOBAL_FORCES && id == molType){ // BONDED_IA_OIF_GLOBAL_FORCES with correct molType  
					test++;
					/* fetch particle 2 */
					p2 = local_particles[p1->bl.e[j++]];
                    if (!p2) {
                      runtimeErrorMsg() <<"oif global calc: bond broken between particles " << p1->p.identity << " and " << p1->bl.e[j-1] << " (particles not stored on the same node - oif_global_forces1); n " << p1->bl.n << " max " << p1->bl.max ;
                      return;
					}
					/* fetch particle 3 */
					//if(n_partners>2){
					p3 = local_particles[p1->bl.e[j++]];
                    if (!p3) {
                      runtimeErrorMsg() <<"oif global calc: bond broken between particles " << p1->p.identity << ", " << p1->bl.e[j-2] << " and " << p1->bl.e[j-1] << " (particles not stored on the same node - oif_global_forces1); n " << p1->bl.n << " max " << p1->bl.max ;
                      return;
					}
					// remaining neighbors fetched
					
					// getting unfolded positions of all particles
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
								printf("Something wrong in oif_global_forces.hpp: All particles in a bond are ghost particles, impossible to unfold the positions...");
								return;
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
					// unfolded positions correct
					VOL_A=area_triangle(p11,p22,p33);
					partArea += VOL_A;

					get_n_triangle(p11,p22,p33,VOL_norm);
					VOL_dn=normr(VOL_norm);
					VOL_hz=1.0/3.0 *(p11[2]+p22[2]+p33[2]);
					VOL_partVol += VOL_A * -1*VOL_norm[2]/VOL_dn * VOL_hz;	
				}
				else{
					j+=n_partners;
				}	
			}
		}
    }
	part_area_volume[0] = partArea;
	part_area_volume[1] = VOL_partVol;
	
	MPI_Allreduce(part_area_volume, area_volume, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}


inline void add_oif_global_forces(double *area_volume, int molType){  //first-fold-then-the-same approach
	double area = area_volume[0];
	double VOL_volume = area_volume[1];
	double VOL_force[3];
	double VOL_A, VOL_norm[3], VOL_dn, VOL_vv; 

	double aa, force1[3], force2[3], force3[3], rh[3], hn, h[3];
	int k;
	
	/** loop over particles */
	int c, np, i ,j;
	Cell *cell;
	Particle *p, *p1, *p2, *p3;
	double p11[3],p22[3],p33[3];
	double AA[3],BB[3];
	int img[3];

	Bonded_ia_parameters *iaparams;
    int type_num, n_partners,id;
    BondedInteraction type;

	int test=0;
	
	/* Loop local cells */
	for (c = 0; c < local_cells.n; c++) {
		cell = local_cells.cell[c];
		p   = cell->part;
		np  = cell->n;
		
		/* Loop cell particles */
		for(i=0; i < np; i++) {				
			j = 0;
			p1=&p[i];
			//printf("i=%d neigh=%d\n", i, p1->bl.n);
			while(j<p1->bl.n){
				/* bond type */
				type_num = p1->bl.e[j++];
				iaparams = &bonded_ia_params[type_num];
				type = iaparams->type;
				n_partners = iaparams->num;
				id=p1->p.mol_id;
				//printf("neigh=%d, type=%d type_num=%d\n", p1->bl.n-1, type, type_num);
				//printf("id %d molType %d\n", id, molType); 
				if(type == BONDED_IA_OIF_GLOBAL_FORCES && id == molType){ // BONDED_IA_OIF_GLOBAL_FORCES with correct molType
					test++;
					/* fetch particle 2 */
					p2 = local_particles[p1->bl.e[j++]];
                                        if (!p2) {
                                          runtimeErrorMsg() <<"add area: bond broken between particles " << p1->p.identity << " and " << p1->bl.e[j-1] << " (particles not stored on the same node - oif_globalforce2); n " << p1->bl.n << " max " << p1->bl.max ;
                                          return;
					}
					/* fetch particle 3 */
					//if(n_partners>2){
					p3 = local_particles[p1->bl.e[j++]];
                                        if (!p3) {
                                          runtimeErrorMsg() <<"add area: bond broken between particles " << p1->p.identity << ", " << p1->bl.e[j-2] << " and " << p1->bl.e[j-1] << " (particles not stored on the same node); n " << p1->bl.n << " max " << p1->bl.max;
                                          return;
					}
					
					// getting unfolded positions of all particles
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
								printf("Something wrong in oif_global_forces.hpp: All particles in a bond are ghost particles, impossible to unfold the positions...");
								return;
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
					// unfolded positions correct
					/// starting code from volume force
					get_n_triangle(p11,p22,p33,VOL_norm);
					VOL_dn=normr(VOL_norm);
					VOL_A=area_triangle(p11,p22,p33);
					VOL_vv=(VOL_volume - iaparams->p.oif_global_forces.V0)/iaparams->p.oif_global_forces.V0;					
					for(k=0;k<3;k++) {
						VOL_force[k]=iaparams->p.oif_global_forces.kv * VOL_vv * VOL_A * VOL_norm[k]/VOL_dn * 1.0 / 3.0;
						//printf("%e ",force[k]);
						p1->f.f[k] += VOL_force[k]; 
						p2->f.f[k] += VOL_force[k];
						p3->f.f[k] += VOL_force[k];
					}
					///  ending code from volume force

					for(k=0;k<3;k++){
						h[k]=1.0/3.0 *(p11[k]+p22[k]+p33[k]);
					}
					
					aa=( area - iaparams->p.oif_global_forces.A0_g) / iaparams->p.oif_global_forces.A0_g;

					//aminusb(3,h,p11,rh);				// area_forces for each triangle node
					vecsub(h,p11,rh);				// area_forces for each triangle node
					hn=normr(rh);
					for(k=0;k<3;k++) {
						force1[k] =  iaparams->p.oif_global_forces.ka_g * aa * rh[k]/hn;
						//(&part1)->f.f[k]+=force[k];
					}
					//aminusb(3,h,p22,rh);				// area_forces for each triangle node
					vecsub(h,p22,rh);				// area_forces for each triangle node
					hn=normr(rh);
					for(k=0;k<3;k++) {
						force2[k] =  iaparams->p.oif_global_forces.ka_g * aa * rh[k]/hn;
						//(&part2)->f.f[k]+=force[k];
					}
					//aminusb(3,h,p33,rh);				// area_forces for each triangle node
					vecsub(h,p33,rh);				// area_forces for each triangle node
					hn=normr(rh);
					for(k=0;k<3;k++) {
						force3[k] =  iaparams->p.oif_global_forces.ka_g * aa * rh[k]/hn;
						//(&part3)->f.f[k]+=force[k];
					}
	
					for(k=0;k<3;k++) {
						p1->f.f[k] += force1[k]; 
						p2->f.f[k] += force2[k];
						p3->f.f[k] += force3[k];
					}

				}
				else{
					j+=n_partners;
				}
			}
		}
    }
	
}

#endif 
