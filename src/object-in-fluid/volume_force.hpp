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
#ifndef _OBJECT_IN_FLUID_VOLUME_FORCE_H
#define _OBJECT_IN_FLUID_VOLUME_FORCE_H
/** \file volume_force.hpp
 *  Routines to calculate the VOLUME_FORCE energy or/and and force 
 *  for a particle triple (triangle from mesh). (Dupin2007)
 *  \ref forces.cpp
*/

#include "utils.hpp"
#include "interaction_data.hpp"
#include "particle_data.hpp"
#include "cells.hpp"
#include "lb.hpp"
#include "grid.hpp"


/** set parameters for the VOLUME_FORCE potential. 
*/
int volume_force_set_params(int bond_type, double V0, double kv);


/************************************************************/

/** called in force_calc() from within forces.cpp
 *  calculates the volume for a cell before the forces are handled
 *  sums up parts for volume with mpi_reduce from Dupin2007 (Equ. 13)
 *  sends volume via mpi_send to nodes
 *
 *  !!! loop over particles from domain_decomposition !!!
 */  

inline void calc_volume(double *volume, int molType){ //first-fold-then-the-same approach
	double partVol=0.,A,norm[3],dn,hz;
	
	/** loop over particles */
	int c, np, i ,j;
	Cell *cell;
	Particle *p, *p1, *p2, *p3;
	double p11[3],p22[3],p33[3];
	int img[3];
	
	Bonded_ia_parameters *iaparams;
	int type_num, type, n_partners, id;
	char *errtxt;

	//int test=0;
	//printf("rank%d, molType2: %d\n", rank,molType);
	/* Loop local cells */
	for (c = 0; c < local_cells.n; c++) {
		cell = local_cells.cell[c];
		p   = cell->part;
		np  = cell->n;
		/* Loop cell particles */
		for(i=0; i < np; i++) {				
			j = 0;
			p1 = &p[i];
			//printf("rank%d, i=%d neigh=%d\n", rank, i, p1->bl.n);
			while(j<p1->bl.n){
				/* bond type */
				type_num = p1->bl.e[j++];				//bond_number
				iaparams = &bonded_ia_params[type_num];
				type = iaparams->type;		  					//type of interaction 14...volume_force
				n_partners = iaparams->num;  					//bonded_neigbours
				id=p1->p.mol_id;						//mol_id of blood cell
				if(type == BONDED_IA_VOLUME_FORCE && id == molType){ // BONDED_IA_VOLUME_FORCE with correct molType   !!!!!!!!!!!!! needs area force local !!!!!!!!!!!!!!!!!!
					p2 = local_particles[p1->bl.e[j++]];
					if (!p2) {
						printf("broken: particles sum %d, id %d, partn %d, bond %d\n", np,id,n_partners,type_num); 
						errtxt = runtime_error(128 + 2*ES_INTEGER_SPACE);
						ERROR_SPRINTF(errtxt,"{volume calc 078 bond broken between particles %d and %d (particles not stored on the same node - volume_force1)} ",
						  p1->p.identity, p1->bl.e[j-1]);
						return;
					}
					/* fetch particle 3 */
					p3 = local_particles[p1->bl.e[j++]];
					if (!p3) {
						errtxt = runtime_error(128 + 3*ES_INTEGER_SPACE);
						ERROR_SPRINTF(errtxt,"{volume calc 079 bond broken between particles %d, %d and %d (particles not stored on the same node); n %d max %d} ",
							p1->p.identity, p1->bl.e[j-2], p1->bl.e[j-1],p1->bl.n,p1->bl.max);
						return;
					}
					memcpy(p11, p1->r.p, 3*sizeof(double));
					memcpy(img, p1->l.i, 3*sizeof(int));
					fold_position(p11, img);
									
					memcpy(p22, p2->r.p, 3*sizeof(double));
					memcpy(img, p2->l.i, 3*sizeof(int));
					fold_position(p22, img);
				
					memcpy(p33, p3->r.p, 3*sizeof(double));
					memcpy(img, p3->l.i, 3*sizeof(int));
					fold_position(p33, img);
				
					get_n_triangle(p11,p22,p33,norm);
					dn=normr(norm);
					A=area_triangle(p11,p22,p33);
					hz=1.0/3.0 *(p11[2]+p22[2]+p33[2]);
					partVol += A * -1*norm[2]/dn * hz;	
				}
				else{
					j+=n_partners;
				}	
			}
		}
    }
	MPI_Allreduce(&partVol, volume, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}

inline void add_volume_force(double volume, int molType){  //first-fold-then-the-same approach
	double A,norm[3],dn; //partVol=0.

	double vv, force[3];
	int k;
	int img[3];
	
	/** loop over particles */
	int c, np, i ,j;
	Cell *cell;
	Particle *p, *p1, *p2, *p3;
	double p11[3],p22[3],p33[3];
	Bonded_ia_parameters *iaparams;
	int type_num, type, n_partners, id;
	char *errtxt;

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
				if(type == BONDED_IA_VOLUME_FORCE && id == molType){ // BONDED_IA_VOLUME_FORCE with correct molType
					test++;
					/* fetch particle 2 */
					p2 = local_particles[p1->bl.e[j++]];
					if (!p2) {
						errtxt = runtime_error(128 + 2*ES_INTEGER_SPACE);
						ERROR_SPRINTF(errtxt,"{volume add 078 bond broken between particles %d and %d (particles not stored on the same node - volume_force2)} ",
						  p1->p.identity, p1->bl.e[j-1]);
						return;
					}
					/* fetch particle 3 */
					p3 = local_particles[p1->bl.e[j++]];
					if (!p3) {
						errtxt = runtime_error(128 + 3*ES_INTEGER_SPACE);
						ERROR_SPRINTF(errtxt,"{volume add 079 bond broken between particles %d, %d and %d (particles not stored on the same node); n %d max %d} ",
							p1->p.identity, p1->bl.e[j-2], p1->bl.e[j-1],p1->bl.n,p1->bl.max);
						return;
					}
					memcpy(p11, p1->r.p, 3*sizeof(double));
					memcpy(img, p1->l.i, 3*sizeof(int));
					fold_position(p11, img);
									
					memcpy(p22, p2->r.p, 3*sizeof(double));
					memcpy(img, p2->l.i, 3*sizeof(int));
					fold_position(p22, img);
				
					memcpy(p33, p3->r.p, 3*sizeof(double));
					memcpy(img, p3->l.i, 3*sizeof(int));
					fold_position(p33, img);
				

					
					get_n_triangle(p11,p22,p33,norm);
					dn=normr(norm);
					A=area_triangle(p11,p22,p33);
					{
				}
					vv=(volume - iaparams->p.volume_force.V0)/iaparams->p.volume_force.V0;
					
					for(k=0;k<3;k++) {
						force[k]=iaparams->p.volume_force.kv * vv * A * norm[k]/dn * 1.0 / 3.0;
						//printf("%e ",force[k]);
						p1->f.f[k] += force[k]; 
						p2->f.f[k] +=  force[k];
						p3->f.f[k] +=  force[k];
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
