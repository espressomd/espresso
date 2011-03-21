/** This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
	It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
	and by which you are legally bound while utilizing this file in any form or way.
	There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
	You should have received a copy of that license along with this program;
	if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
	write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
	Copyright (c) 2002-2006; all rights reserved unless otherwise stated.
	\file maggs.c  
	Local Maggs algorithm for long range coulomb interaction.
 */
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"
#include "global.h"
#include "grid.h"
#include "integrate.h"
#include "initialize.h"
#include "interaction_data.h"
#include "particle_data.h"
#include "communication.h"
#include "maggs.h"
#include "thermostat.h"
#include "cells.h"
#include "domain_decomposition.h"

/* MPI tags for the maggs communications: */
/** Tag for communication in Maggs_init() -> calc_glue_patch(). */
#define REQ_MAGGS_SPREAD 300
#define REQ_MAGGS_EQUIL  301
#define ELECTROSTATICS /* just so XCode indexes the code! Later remove! */

#ifdef ELECTROSTATICS /* later to remove!!!! */

/************************************************
 * data types
 ************************************************/

typedef int      t_ivector [SPACE_DIM];
typedef double   t_dvector [SPACE_DIM];
typedef int      t_dirs    [NDIRS];

/** Structure of local lattice parameters. */
typedef struct {
	/* local mesh characterization. */
	t_dvector ld_pos;        /** spacial position of left down grid point */  
	t_dvector ur_pos;        /** spacial positon of upper right grid point */ 
	t_dvector low_bound;     /** low bound for the coord of particles of ghost region */ 
	t_dvector up_bound;      /** upper bound for the coord of particles of ghost region */
	int globin_ld[3];        /** inner left down global grid point    */
	int in_ld[3];            /** inner left down grid point    */
	int in_ur[3];            /** inner up right grid point + (1,1,1) */
	int halo_ld[3];          /** halo-region left down grid point  */
	int halo_ur[3];          /** halo-region up right global grid point  */
	int margin[SPACE_DIM*2]; /** number of margin mesh points (even index - left, odd - right). */
	int r_margin[6];         /** number of margin mesh points from neighbour nodes */
	t_ivector dim;           /** grid dimension (size + glue_patch region) of local mesh.  */
	t_ivector size;          /** dimension of mesh inside node domain.          */
	int volume;
	int inner_vol;
} lattice_param;

/** surface_patch structure */
typedef struct {
	int offset;        /** source offset for the site index */
	int doffset;       /** desitnation offset for the site index */
	int stride;        /** minimal contiguous block */  
	int skip;          /** gap between two strides (from the first element of one stride 
	 ** to the first elem. of next stride
	 */
	int nblocks;       /** number of strides */
	int coord[2];      /** coordinates of the vector fields which has to be exchanged */
	int volume;
} t_surf_patch;

typedef struct {
	double    charge;
	double    permittivity[SPACE_DIM];            /* dielectric properties of the site. */
	short     r[SPACE_DIM];
} t_site;



/************************************************
 * variables
 ************************************************/
int IND_ORDER[SPACE_DIM] = {2,1,0};  /* loop order for dimensions (for ESPRESSO z is shortest) */

//static double maggs_pref1;
//static double maggs_pref2;
MAGGS_struct maggs = { 0. ,0. ,0. ,0. ,0. ,0 , 0., 0.};

/** local mesh. */
static lattice_param lparam;
/** local lattice */
static t_site* lattice;
/** local E field */
static double* Efield;
/** local B field */
static double* Bfield;
/** site neighbors */
static t_dirs* neighbor;


/** Array to store glue_patch data to send. */
double *send_databuf = NULL; 
/** Array to store glue_patch data to recv */
double *recv_databuf = NULL;

/** Private Functions */
/************************************************************/

void maggs_interpolate_charge(int *first, double *rel, double q);
void maggs_accumulate_charge_from_ghosts();
void maggs_calc_self_influence(Particle* P);
void maggs_calc_charge_gradients(double *rel, double q, double *grad);
short maggs_check_intersect_1D(double delta, double r_new, int dir, int first, double *t_step, int identity);
void maggs_calc_e_field_on_link_1D(int index, double *flux, double v, int dir);
void maggs_setup_neighbors();


int maggs_get_linear_index(int x, int y, int z, int latticedim[SPACE_DIM])
/** Turns the 3D index into a linear index. */
/** adim contains the dimensions of the lattice in */
/** the 3 directions: adim[0] = x_max, ..., adim[2] = z_max */
/** z is the first index!!! */
{
	return (z + latticedim[ZPLUS]*(y + latticedim[YPLUS]*x));
}

int maggs_count_charged_particles()
{  
	Cell *cell;
	Particle *part;
	int i,c,np;
	double node_sum, tot_sum;
	
	node_sum=0.0; 
	tot_sum =0.0;
	
	for (c = 0; c < local_cells.n; c++) {
		cell = local_cells.cell[c];
		part = cell->part;
		np   = cell->n;
		for(i=0;i<np;i++) 
			if( part[i].p.q != 0.0 ) node_sum += 1.0;
	}
	
	MPI_Reduce(&node_sum, &tot_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	
	return tot_sum;
}


int maggs_get_offset(int index_shift, int index_base, int axes, int adim[3])
/** Index shift is calculated for moving in various */
/** directions with the linear index. */
{
	int dif;
	dif = index_shift - index_base;
	if(axes <= 1) dif *= adim[2];
	if(axes == 0) dif *= adim[1]; 
	
	return (dif);
}

double maggs_interpol1D(double x)
/** Interpolation function in one dimension. */
/** Currently only linear interpolation conserves charge. */
{	
	return x;
	/* Cosine interpolation would be: */
	/* return sqr(sin(M_PI_2*x)); */
}

void maggs_calc_directions(int j, int* dir1, int*dir2)
/** For any direction j, write the two other directions */
/** into the pointers in circular permutation. */
{
	*dir1 = *dir2 = -1;
	switch(j) {
		case 0 :
			*dir1 = 2;
			*dir2 = 1;
			break;
		case 1 :
			*dir1 = 2;
			*dir2 = 0;
			break;  
		case 2 :
			*dir1 = 1;
			*dir2 = 0;
			break;
	}  
}

double maggs_calc_curl(int mue, int nue, double* field, int* Neighbor, int index)
/** Calculates the finite differences rotation in real space in mue-nue plane: */
/** \dot{D} = \nabla \times B (and prefactors plus current)*/
/** The given "double* field" should be a B-Field!! */
{
	double result;
	
	result = field[index+mue] + field[3*Neighbor[OPP_DIR(mue)]+nue] -
    field[3*Neighbor[OPP_DIR(nue)]+mue] - field[index+nue];
	
	return result;
}

double maggs_calc_dual_curl(int mue, int nue, double* field, int* Neighbor, int index)
/** Calculates the finite differences rotation in dual space in mue-nue plane: */
/** \dot{B} = - \nabla \times D / \epsilon */
/** The given "double* field" should be a D-Field!! */
{
	double res;
	
	res = field[index+mue] + field[3*Neighbor[mue]+nue] -
	field[3*Neighbor[nue]+mue] - field[index+nue];
	
	return res;
}


void maggs_update_plaquette(int mue, int nue, int* Neighbor, int index, double delta)
/** updates all E-fields on links of the plaquette */
/** and the surroundings. delta was calculated before */
/** in function "maggs_perform_rot_move_inplane". */
{
	int i = 3*index;
	Efield[i+mue]             += delta;
	Efield[3*Neighbor[mue]+nue] += delta;
	Efield[3*Neighbor[nue]+mue] -= delta;
	Efield[i+nue]             -= delta;  
}


void maggs_setup_local_lattice()
/** Set up lattice, calculate dimensions and lattice parameters */
/** Allocate memory for lattice sites and fields */
{
	int i;
	int ix = 0;
	int iy = 0;
	int iz = 0;
	int linearindex = 0;
	int xyzcube;	
	
	xyzcube = 1;
	FOR3D(i) {
		/* inner left down grid point (global index) */
		lparam.in_ld[i] = (int)ceil(my_left[i]*maggs.inva); 
		/* inner up right grid point (global index) */
		lparam.in_ur[i] = (int)floor(my_right[i]*maggs.inva); 
		/* correct roundof errors at boundary */
		if(my_right[i]*maggs.inva-lparam.in_ur[i]<ROUND_ERROR_PREC) lparam.in_ur[i]--;
		if(1.0+my_left[i]*maggs.inva-lparam.in_ld[i]<ROUND_ERROR_PREC) lparam.in_ld[i]--;
		lparam.globin_ld[i] = lparam.in_ld[i];
		/* inner grid dimensions */
		lparam.size[i] = lparam.in_ur[i] - lparam.in_ld[i] + 1;
		/* spacial position of left down grid point */
		lparam.ld_pos[i] = my_left[i] - maggs.a;  
		/* spacial position of upper right grid point */
		lparam.ur_pos[i] = my_right[i] + maggs.a;  
		/* left down margin */
		lparam.margin[i*2] = 1;
		/* up right margin */
		lparam.margin[(i*2)+1] = 1;
		
		lparam.dim[i] = lparam.size[i] + lparam.margin[i*2] + lparam.margin[i*2+1];
		xyzcube *= lparam.dim[i];
		/* reduce inner grid indices from global to local */
		lparam.in_ld[i] = lparam.margin[i*2];
		lparam.in_ur[i] = lparam.margin[i*2]+lparam.size[i];
		lparam.halo_ld[i] = 0;
		lparam.halo_ur[i] = lparam.in_ur[i];      
	}
	
	
	lparam.volume    = xyzcube;
	lparam.inner_vol = lparam.size[0]*lparam.size[1]*lparam.size[2];
	/* allocate memory for sites and neighbors */
	lattice  = (t_site*) malloc(xyzcube*sizeof(t_site));
	neighbor = (t_dirs*) malloc(xyzcube*sizeof(t_dirs));
	
	Bfield   = (double*) malloc(3*xyzcube*sizeof(double));
	Efield   = (double*) malloc(3*xyzcube*sizeof(double));
	
	/* set up lattice sites */
	FORALL_SITES(ix, iy, iz) {
		linearindex = maggs_get_linear_index(ix, iy, iz, lparam.dim);
		
		lattice[linearindex].r[0] = ix;
		lattice[linearindex].r[1] = iy;
		lattice[linearindex].r[2] = iz;

		FOR3D(i) {
			Bfield[3*linearindex+i]  = 0.;
			Efield[3*linearindex+i]  = 0.;
		}

		lattice[linearindex].charge = 0.;
		
		/* Here, we need a function to set PERMITTIVITY!!!!! */
		FOR3D(i) lattice[linearindex].permittivity[i] = 1.;
	}
	
	maggs_setup_neighbors();
}

void maggs_setup_neighbors()
/** sets up nearest neighbors for each site in linear index */
{
	int ix = 0;
	int iy = 0;
	int iz = 0;
	
	int xsize = lparam.dim[0];
	int ysize = lparam.dim[1];
	int zsize = lparam.dim[2];
	
	int ixplus = 0;
	int ixminus = 0;
	int iyplus = 0;
	int iyminus = 0;
	int izplus = 0;
	int izminus = 0;
	
	int kount = 0;
	
	int kountxplus = 0;
	int kountxminus = 0;
	int kountyplus = 0;
	int kountyminus = 0;
	int kountzplus = 0;
	int kountzminus = 0;
	
	for (ix = 0; ix < xsize; ix++) 
    {
		ixplus  = ix + 1;
		ixminus = ix - 1;
		for(iy = 0; iy < ysize; iy ++)
		{
			iyplus  = iy + 1;
			iyminus = iy - 1;
			for(iz = 0; iz < zsize; iz ++)
			{
				izplus  = iz + 1;
				izminus = iz - 1;
				
				kount         = maggs_get_linear_index(ix,      iy,      iz,      lparam.dim);
				kountzplus    = maggs_get_linear_index(ix,      iy,      izplus,  lparam.dim);
				kountzminus   = maggs_get_linear_index(ix,      iy,      izminus, lparam.dim);
				kountyplus    = maggs_get_linear_index(ix,      iyplus,  iz,      lparam.dim);
				kountyminus   = maggs_get_linear_index(ix,      iyminus, iz,      lparam.dim);
				kountxplus    = maggs_get_linear_index(ixplus,  iy,      iz,      lparam.dim);
				kountxminus   = maggs_get_linear_index(ixminus, iy,      iz,      lparam.dim);
				
				if(ixminus < 0)     neighbor[kount][XMINUS] = -1;
				else                neighbor[kount][XMINUS] = kountxminus;
				if(ixplus >= xsize) neighbor[kount][XPLUS]  = -1;
				else                neighbor[kount][XPLUS]  = kountxplus;
				
				if(iyminus < 0)     neighbor[kount][YMINUS] = -1;
				else                neighbor[kount][YMINUS] = kountyminus;
				if(iyplus >= ysize) neighbor[kount][YPLUS]  = -1;
				else                neighbor[kount][YPLUS]  = kountyplus;
				
				if(izminus < 0)     neighbor[kount][ZMINUS] = -1;
				else                neighbor[kount][ZMINUS] = kountzminus;
				if(izplus >= zsize) neighbor[kount][ZPLUS]  = -1;
				else                neighbor[kount][ZPLUS]  = kountzplus;
			}
		}
    }
	return;
}

void maggs_calc_surface_patches(t_surf_patch* surface_patch)
/** sets up surface patches for all domains. */
{
	/* x=lparam.size[0] plane */
	surface_patch[0].offset   = lparam.dim[2]*lparam.dim[1]*lparam.size[0];    /*(size[0],0,0) point */
	surface_patch[0].doffset  = 0;                                             /*(0,0,0) point */
	surface_patch[0].stride   = lparam.dim[2]*lparam.dim[1];
	surface_patch[0].skip     = 0;
	surface_patch[0].nblocks  = 1;
	surface_patch[0].coord[0] = 2;
	surface_patch[0].coord[1] = 1;
	surface_patch[0].volume   = lparam.dim[2]*lparam.dim[1];
	
	/* x=1 plane */
	surface_patch[1].offset   = lparam.dim[2]*lparam.dim[1];                    /*(1,0,0) point */
	surface_patch[1].doffset  = lparam.dim[2]*lparam.dim[1]*lparam.in_ur[0];    /*(halo[0],0,0) point */
	surface_patch[1].stride   = lparam.dim[2]*lparam.dim[1];
	surface_patch[1].skip     = 0;
	surface_patch[1].nblocks  = 1;
	surface_patch[1].coord[0] = 2;
	surface_patch[1].coord[1] = 1;
	surface_patch[1].volume   = lparam.dim[2]*lparam.dim[1]; 
	
	/* y=lparam.size[1] plane */
	surface_patch[2].offset   = lparam.dim[2]*lparam.size[1];               /*(0,size[1],0) point */
	surface_patch[2].doffset  = 0;                                          /*(0,0,0) point */
	surface_patch[2].stride   = lparam.dim[2];
	surface_patch[2].skip     = lparam.dim[2]*lparam.dim[1];
	surface_patch[2].nblocks  = lparam.dim[0];  
	surface_patch[2].coord[0] = 2;
	surface_patch[2].coord[1] = 0;
	surface_patch[2].volume   = lparam.dim[2]*lparam.dim[0];
	
	/* y=1 plane */
	surface_patch[3].offset   = lparam.dim[2];                             /*(0,1,0) point */
	surface_patch[3].doffset  = lparam.dim[2]*lparam.in_ur[1];             /*(0,in_ur[1],0) point */
	surface_patch[3].stride   = lparam.dim[2];
	surface_patch[3].skip     = lparam.dim[2]*lparam.dim[1];
	surface_patch[3].nblocks  = lparam.dim[0];
	surface_patch[3].coord[0] = 2;
	surface_patch[3].coord[1] = 0;
	surface_patch[3].volume   = lparam.dim[2]*lparam.dim[0];
	
	/* z=lparam.size[2] plane */
	surface_patch[4].offset   = lparam.size[2];    /*(0,0,size[2]) point */
	surface_patch[4].doffset  = 0;                 /*(0,0,0) point */
	surface_patch[4].stride   = 1;
	surface_patch[4].skip     = lparam.dim[2];
	surface_patch[4].nblocks  = lparam.dim[0]*lparam.dim[1];
	surface_patch[4].coord[0] = 1;
	surface_patch[4].coord[1] = 0;
	surface_patch[4].volume   = lparam.dim[0]*lparam.dim[1];
	
	/* z=1 plane for z it must be higher*/
	surface_patch[5].offset   = 1;                   /*(0,0,1) point */
	surface_patch[5].doffset  = lparam.in_ur[2];     /*(0,0,in_ur[2]) point */
	surface_patch[5].stride   = 1;
	surface_patch[5].skip     = lparam.dim[2];
	surface_patch[5].nblocks  = lparam.dim[0]*lparam.dim[1];
	surface_patch[5].coord[0] = 1;
	surface_patch[5].coord[1] = 0;
	surface_patch[5].volume   = lparam.dim[0]*lparam.dim[1];
}


void maggs_prepare_surface_planes(int dim, MPI_Datatype *xy, MPI_Datatype *xz, MPI_Datatype *yz, 
							t_surf_patch *surface_patch)
/** sets up MPI communications for domain surfaces */
{
	MPI_Type_contiguous(dim*surface_patch[0].stride*sizeof(double),MPI_BYTE,yz);  
	MPI_Type_commit(yz);
	MPI_Type_vector(surface_patch[4].nblocks, dim*surface_patch[4].stride,
					dim*surface_patch[4].skip, MPI_DOUBLE,xy);
	MPI_Type_commit(xy); 
	MPI_Type_vector(surface_patch[2].nblocks, dim*surface_patch[2].stride,
					dim*surface_patch[2].skip, MPI_DOUBLE,xz);
	MPI_Type_commit(xz);
}


void maggs_exchange_surface_patch(double *field, int dim, int e_equil)
/** MPI communication of surface region. */
/** works for D- and B-fields. */
{
	static int init = 1;
	static int flag_free = 1;
	static MPI_Datatype xyPlane,xzPlane,yzPlane; 
	static MPI_Datatype xzPlane2D, xyPlane2D, yzPlane2D;
	/*  int coord[2]; */
	int l, s_dir, r_dir;
	/*  int pos=0; */
	MPI_Status status[2];
	MPI_Request request[]={MPI_REQUEST_NULL, MPI_REQUEST_NULL};
	int offset, doffset, skip, stride, nblocks;
	/** surface_patch */
	static t_surf_patch  surface_patch[6];
	
	if(init) {
		MPI_Datatype xz_plaq, oneslice;
		
		maggs_calc_surface_patches(surface_patch);
		maggs_prepare_surface_planes(dim, &xyPlane, &xzPlane, &yzPlane, surface_patch);
		
		MPI_Type_vector(surface_patch[0].stride, 2, 3, MPI_DOUBLE,&yzPlane2D);    
		MPI_Type_commit(&yzPlane2D);
		
		/* create data type for xz plaquette */
		MPI_Type_hvector(2,1*sizeof(double),2*sizeof(double), MPI_BYTE, &xz_plaq);
		/* create data type for a 1D section */
		MPI_Type_contiguous(surface_patch[2].stride, xz_plaq, &oneslice); 
		/* create data type for a 2D xz plane */
		MPI_Type_hvector(surface_patch[2].nblocks, 1, dim*surface_patch[2].skip*sizeof(double), oneslice, &xzPlane2D);
		MPI_Type_commit(&xzPlane2D);    
		/* create data type for a 2D xy plane */
		MPI_Type_vector(surface_patch[4].nblocks, 2, dim*surface_patch[4].skip, MPI_DOUBLE, &xyPlane2D);
		MPI_Type_commit(&xyPlane2D); 
		
		init = 0;
	}
	
		
	/* direction loop */
	for(s_dir=0; s_dir < 6; s_dir++) { 
		offset = dim * surface_patch[s_dir].offset;
		doffset= dim * surface_patch[s_dir].doffset;
		
		if(s_dir%2==0) r_dir = s_dir+1;
		else           r_dir = s_dir-1;
		/* pack send halo-plane data */
		if(node_neighbors[s_dir] != this_node) {
			/* communication */
			switch(s_dir) {
				case 0 :
				case 1 :
					if(e_equil || dim == 1) {
						MPI_Irecv (&field[doffset],1,yzPlane,node_neighbors[s_dir],REQ_MAGGS_SPREAD,MPI_COMM_WORLD,&request[0]);
						MPI_Isend(&field[offset],1,yzPlane,node_neighbors[r_dir],REQ_MAGGS_SPREAD,MPI_COMM_WORLD,&request[1]);
					}
					else {
						MPI_Irecv (&field[doffset+1],1,yzPlane2D,node_neighbors[s_dir],REQ_MAGGS_SPREAD,MPI_COMM_WORLD,&request[0]);
						MPI_Isend(&field[offset+1],1,yzPlane2D,node_neighbors[r_dir],REQ_MAGGS_SPREAD,MPI_COMM_WORLD,&request[1]);
					}	  
					
					MPI_Waitall(2,request,status);
					break;
				case 2 :
				case 3 :
					if(e_equil || dim == 1) {
						MPI_Irecv (&field[doffset],1,xzPlane,node_neighbors[s_dir],REQ_MAGGS_SPREAD,MPI_COMM_WORLD,&request[0]);
						MPI_Isend(&field[offset],1,xzPlane,node_neighbors[r_dir],REQ_MAGGS_SPREAD,MPI_COMM_WORLD,&request[1]);
					}
					else {
						MPI_Irecv (&field[doffset],1,xzPlane2D,node_neighbors[s_dir],REQ_MAGGS_SPREAD,MPI_COMM_WORLD,&request[0]);
						MPI_Isend(&field[offset],1,xzPlane2D,node_neighbors[r_dir],REQ_MAGGS_SPREAD,MPI_COMM_WORLD,&request[1]);
					}	  
					MPI_Waitall(2,request,status);
					break;
				case 4 :
				case 5 : 
					if(e_equil || dim == 1) {
						MPI_Irecv (&field[doffset],1,xyPlane,node_neighbors[s_dir],REQ_MAGGS_SPREAD,MPI_COMM_WORLD,&request[0]);
						MPI_Isend(&field[offset],1,xyPlane,node_neighbors[r_dir],REQ_MAGGS_SPREAD,MPI_COMM_WORLD,&request[1]);
					}
					else {
						MPI_Irecv (&field[doffset],1,xyPlane2D,node_neighbors[s_dir],REQ_MAGGS_SPREAD,MPI_COMM_WORLD,&request[0]);
						MPI_Isend(&field[offset],1,xyPlane2D,node_neighbors[r_dir],REQ_MAGGS_SPREAD,MPI_COMM_WORLD,&request[1]);
					}
					MPI_Waitall(2,request,status);
					break;
			}
		}

		else {
			/** copy locally */
			skip    = dim * surface_patch[s_dir].skip;
			stride  = dim * surface_patch[s_dir].stride * sizeof(double);
			nblocks = surface_patch[s_dir].nblocks;
			
			for(l=0; l<nblocks; l++){
				memcpy(&(field[doffset]), &(field[offset]), stride);
				offset  += skip;
				doffset += skip;
			}
			
		}
	}
}

void  maggs_distribute_particle_charges()
/** finds current lattice site of each particle. */
/** calculates charge interpolation on cube. */
{	
	Cell *cell;
	Particle* p;
	int i, c, d;
	int np;
	int first[SPACE_DIM];
	double q;
	double pos[SPACE_DIM], rel[SPACE_DIM];
		
	for(i=0;i<lparam.volume;i++) lattice[i].charge = 0.;
	/* === charge assignment === */ 
	/* loop over inner cells */
	for (c = 0; c < local_cells.n; c++) {
		cell = local_cells.cell[c];
		p  = cell->part;
		np = cell->n;
		for(i = 0; i < np; i++) {
			if( (q=p[i].p.q) != 0.0 ) {
				FOR3D(d) {
					pos[d]        = (p[i].r.p[d] - lparam.ld_pos[d])* maggs.inva;
					first[d]      = (int) pos[d];
					rel[d]        = pos[d] - first[d];
				}
				maggs_interpolate_charge(first, rel, q);
			}
		}      
	}
	maggs_accumulate_charge_from_ghosts();	
}
	

void maggs_accumulate_charge_from_ghosts()
/** add charges from ghost cells to lattice sites. */
{
	Cell *cell;
	Particle* p;
	int i, c, d;
	int np;
	int flag_inner=0;
	int first[SPACE_DIM];
	double q;
	double pos[SPACE_DIM], rel[SPACE_DIM];
		
	/* loop over ghost cells */
	for (c = 0; c < ghost_cells.n; c++) {
		cell = ghost_cells.cell[c];
		p  = cell->part;
		np = cell->n;
		for(i = 0; i < np; i++) {
			if( (q=p[i].p.q) != 0.0 ) {
				flag_inner=1;
				FOR3D(d) {
					if(p[i].r.p[d]<lparam.ld_pos[d]||p[i].r.p[d]>=lparam.ur_pos[d])
					{flag_inner=0; break;}
				}
			}
			if(flag_inner) {
				FOR3D(d) {
					pos[d]        = (p[i].r.p[d] - lparam.ld_pos[d])* maggs.inva;
					first[d]      = (int) pos[d];
					rel[d]        = pos[d] - first[d];
				}
				maggs_interpolate_charge(first, rel, q);
			}
		}      
	}  

}


void  maggs_calc_grad_and_point_forces(double *grad)
/** finds correct lattice cube for particles. */
/** calculates charge gradients on this cube. */
{
	Cell *cell;
	Particle* p;
	int i, c, d, ip;
	int np;
	int index = -1;
	int first[SPACE_DIM];
	double q;
	double pos[SPACE_DIM], rel[SPACE_DIM];
		
	/* === grad assignment for real particles and self-force === */ 
	ip = 0;
	for (c = 0; c < local_cells.n; c++) {
		cell = local_cells.cell[c];
		p  = cell->part;
		np = cell->n;
		for(i = 0; i < np; i++) {
			if( (q=p[i].p.q) != 0.0 ) {
				FOR3D(d) {
					pos[d]        = (p[i].r.p[d] - lparam.ld_pos[d])* maggs.inva;
					first[d]      = (int) pos[d];
					rel[d]        = pos[d] - first[d];
				}
				maggs_calc_charge_gradients(rel, q, &grad[ip]);
				ip += 12;
			}
		}      
	}  
	
} 

double maggs_check_curl_E()
/** For energy minimization: */
/** returns maximum of curl(E) from all sites. */
/** May also be used to verify constraint surface condition. */
{
	int i, ix, iy, iz;
	double curl, maxcurl, gmaxcurl;
	int* anchor_neighb;
	
	maxcurl = 0.;
	
	FORALL_INNER_SITES(ix, iy, iz) {
		i = maggs_get_linear_index(ix, iy, iz, lparam.dim); 
		anchor_neighb = neighbor[i];
		curl = Efield[3*i] + Efield[3*anchor_neighb[0]+1] 
		- Efield[3*anchor_neighb[1]] - Efield[3*i+1];
		curl *= maggs.inva;
		if(fabs(curl)>maxcurl) maxcurl = fabs(curl);
		curl = Efield[3*i+2] + Efield[3*anchor_neighb[2]] 
		- Efield[3*anchor_neighb[0]+2] - Efield[3*i];
		curl *= maggs.inva;
		if(fabs(curl)>maxcurl) maxcurl = fabs(curl);
		curl = Efield[3*i+1] + Efield[3*anchor_neighb[1]+2] 
		- Efield[3*anchor_neighb[2]+1] - Efield[3*i+2];
		curl *= maggs.inva;
		if(fabs(curl)>maxcurl) maxcurl = fabs(curl);
	}
	MPI_Allreduce(&maxcurl,&gmaxcurl,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);  
	return gmaxcurl;
}

void maggs_perform_rot_move_inplane(int i, int n)
/** If curl(E) is not zero, update plaquette with corrections. */
/** coordinate "int n" is the normal direction to the plaquette */
{
	int mue, nue;
	int * anchor_neighb;
	double delta;
	double ROUND_ERR = 0.01*ROUND_ERROR_PREC;
	
	mue = 0; nue = 0;
	
	switch(n) {
		case 0 :
			mue = 1;
			nue = 2;
			break;
			
		case 1 :
			mue = 2;
			nue = 0;
			break;
		case 2 :
			mue = 0;
			nue = 1;
			break;
	}
	
	anchor_neighb = &neighbor[i][0];
	
	delta = Efield[3*i+mue] + Efield[3*anchor_neighb[mue]+nue] 
    - Efield[3*anchor_neighb[nue]+mue] - Efield[3*i+nue];
	if(fabs(delta)>=ROUND_ERR) {
		delta = -delta/4.; 
		maggs_update_plaquette(mue, nue, anchor_neighb, i, delta);
	}
}


void maggs_minimize_transverse_field()
/** For energy minimization: */
/** Relax B-fields in all directions. */
{
	int k, l, m;
	int i, d;
	int ind_i, ind_j;
	int size[2]={0,0};
	int index = -1; /* force allocation error */
	
	FOR3D(d) {
		switch(d) {
			case 0 :
				size[0] = lparam.size[2];
				size[1] = lparam.size[1];
				break;
			case 1 :
				size[0] = lparam.size[2];
				size[1] = lparam.size[0];
				break;
			case 2 :
				size[0] = lparam.size[1];
				size[1] = lparam.size[0];
				break;
		}
		for(i=0;i<2;i++) {
			/* at first even sites (i==0) then odd */
			for(k=1;k<=lparam.size[d];k++) {
				/* update every plane in direction d */
				ind_i=0;
				for(l=0; l<=size[1]; l++){
					ind_j=0;
					for(m=0;m<=size[0]; m++) {
						switch(d) {
							case 0 :
								index=maggs_get_linear_index(k,l,m,lparam.dim);
								break;
							case 1 :
								index=maggs_get_linear_index(l,k,m,lparam.dim); 
								break;
							case 2 :
								index=maggs_get_linear_index(l,m,k,lparam.dim); 
								break;
						}
						if((ind_i+ind_j)%2==i)
							maggs_perform_rot_move_inplane(index, d);
						ind_j++;
					}
					ind_i++;
				}   
			}
			/* update boundaries - update halo regions */
			maggs_exchange_surface_patch(Efield, 3, 0);
		}
	}
}


void maggs_calc_init_e_field()
/** calculates initial electric field configuration. */
/** currently uses simple and slow method of plaquettes and links. */
/** energy minimization takes up lots of time. */
{
	int xsizeplus, ysizeplus;
	double localqy, localqz;
	double qplane, qline;
	int    i, k, ix, iy, iz;
	int index = 0;
	double sqrE, invasq, tmp_field=0.0;
	double gsqrE, goldE, gavgEx, gavgEy, gavgEz;
	double qz, qy, qx, avgEx, avgEy, avgEz;
	double Eall[SPACE_DIM], gEall[SPACE_DIM], maxcurl;
	MPI_Status status;
	MPI_Comm zplane, yline;
	int color, rank, dim;
	
	MAGGS_TRACE(fprintf(stderr,"%d:Initialize field\n",this_node));
	
	invasq = maggs.inva*maggs.inva;
	
	xsizeplus = lparam.dim[0];
	ysizeplus = lparam.dim[1];
	
	/* sort particles for the calculation of initial charge distribution */
	MAGGS_TRACE(fprintf(stderr,"%d:Sorting particles...\n",this_node));
	
	cells_resort_particles(CELL_GLOBAL_EXCHANGE);
	/*  fprintf(stderr, "done\n"); */
	maggs_distribute_particle_charges();
	
	dim = node_grid[1]*node_grid[0];
	color = this_node/dim;
	rank  = this_node%dim;
	MPI_Comm_split(MPI_COMM_WORLD, color, rank, &zplane);
	color = rank/node_grid[0];
	rank  = rank%node_grid[0];
	MPI_Comm_split(zplane, color, rank, &yline);
	
	/* calculate initial solution of Poisson equation */
	
	/** CAUTION: the indexing of the neighbor nodes in Espresso
	 *  starts from x left neighbor node
	 */
	
	/* get process coordinates */
	if(node_pos[2]!= 0) {
		MPI_Recv(&tmp_field, 1, MPI_DOUBLE, node_neighbors[4], REQ_MAGGS_EQUIL, MPI_COMM_WORLD, &status);
		for(iy=lparam.in_ld[1];iy<lparam.in_ur[1];iy++) {
			for(ix=lparam.in_ld[0];ix<lparam.in_ur[0];ix++) {  
				index = maggs_get_linear_index(ix, iy, lparam.in_ld[2], lparam.dim);
				Efield[3*neighbor[index][ZMINUS]+ZPLUS] = tmp_field;
			}
		}
	}

	localqz = 0.;
	for(iz=lparam.in_ld[2];iz<lparam.in_ur[2];iz++) { /* loop over z-planes */
		localqz = 0.;
		for(iy=lparam.in_ld[1];iy<lparam.in_ur[1];iy++) {
			for(ix=lparam.in_ld[0];ix<lparam.in_ur[0];ix++) {  
				index = maggs_get_linear_index(ix, iy, iz, lparam.dim);
				/* Sum over the charge of all sides in z-plane */
				localqz += lattice[index].charge / lattice[index].permittivity[2];
			}
		} 
		
		MPI_Allreduce(&localqz, &qz, 1, MPI_DOUBLE, MPI_SUM, zplane);
		qz = qz/(maggs.mesh*maggs.mesh);
		qplane = qz*maggs.prefactor*invasq;
		/*    if(fabs(qplane) >= 0.01*ROUND_ERROR_PREC) { */
		for(iy=lparam.in_ld[1];iy<lparam.in_ur[1];iy++) {
			for(ix=lparam.in_ld[0];ix<lparam.in_ur[0];ix++) {  
				index = maggs_get_linear_index(ix, iy, iz, lparam.dim);
				Efield[3*index+ZPLUS]  = Efield[3*neighbor[index][ZMINUS]+ZPLUS] + qplane;
				/*	    + qz*maggs.prefactor*invasq; */
			}
		}
		/*    } */
		if(iz>=lparam.in_ur[2]-1) {
			if (node_pos[2]<node_grid[2]-1) {
				if(node_grid[2]>1) {
					MPI_Send(&Efield[3*index+ZPLUS], 1, MPI_DOUBLE, node_neighbors[5], REQ_MAGGS_EQUIL, MPI_COMM_WORLD); 
				}
			}
			else 
				if (fabs(Efield[3*index+ZPLUS]) > 10.*ROUND_ERROR_PREC) {
					fprintf(stderr, "%d: Error in the calculation of Ez(%d,%d,%d)=%f!!\n", 
							this_node,lattice[index].r[0], lattice[index].r[1], lattice[index].r[2],
							Efield[3*index+ZPLUS]);
					fflush(stderr);
				}
		}
		
		if(node_pos[1]!= 0) {
			MPI_Recv(&tmp_field, 1, MPI_DOUBLE, node_neighbors[2], REQ_MAGGS_EQUIL, MPI_COMM_WORLD, &status);
			for(ix=lparam.in_ld[0];ix<lparam.in_ur[0];ix++) {  
				index = maggs_get_linear_index(ix, lparam.in_ld[1], iz, lparam.dim);
				Efield[3*neighbor[index][YMINUS]+YPLUS] = tmp_field;
			}
		}
		
		for(iy=lparam.in_ld[1];iy<lparam.in_ur[1];iy++) {
			localqy = 0.;
			for(ix=lparam.in_ld[0];ix<lparam.in_ur[0];ix++) {  
				index = maggs_get_linear_index(ix, iy, iz, lparam.dim);
				localqy += lattice[index].charge / lattice[index].permittivity[1];
			}
			
			MPI_Allreduce(&localqy, &qy, 1, MPI_DOUBLE, MPI_SUM, yline);
			
			qy = qy/maggs.mesh;
			qline = (qy-qz)*maggs.prefactor*invasq;
			/*      if(fabs(qy-qz)>=ROUND_ERROR_PREC) { */
			for(ix=lparam.in_ld[0];ix<lparam.in_ur[0];ix++) {  
				index = maggs_get_linear_index(ix, iy, iz, lparam.dim);
				Efield[3*index+YPLUS]  = Efield[3*neighbor[index][YMINUS]+YPLUS] + qline;
				/*	    (qy-qz)*maggs.prefactor*invasq; */
			}
			/*      } */
			
			if(iy>=lparam.in_ur[1]-1) {
				if(node_pos[1] < node_grid[1]-1) {
					if (node_grid[1]>1)
						MPI_Send(&Efield[3*index+YPLUS], 1, MPI_DOUBLE, node_neighbors[3], REQ_MAGGS_EQUIL, MPI_COMM_WORLD); 
				}
				else
					if (fabs(Efield[3*index+YPLUS]) > 10.*ROUND_ERROR_PREC)
						fprintf(stderr, "%d: Error in the calculation of Ey(%d,%d,%d)=%f!!\n",
								this_node, lattice[index].r[0], lattice[index].r[1], lattice[index].r[2],
								Efield[3*index+YPLUS]);	  
			}
			
			if(node_pos[0]!= 0) {
				MPI_Recv(&tmp_field, 1, MPI_DOUBLE, node_neighbors[0], REQ_MAGGS_EQUIL, MPI_COMM_WORLD, &status);
				index = maggs_get_linear_index(lparam.in_ld[0], iy, iz, lparam.dim);
				Efield[3*neighbor[index][XMINUS]+XPLUS] = tmp_field;
			}
			
			for(ix=lparam.in_ld[0];ix<lparam.in_ur[0];ix++) {  
				index = maggs_get_linear_index(ix, iy, iz, lparam.dim);
				qx = lattice[index].charge / lattice[index].permittivity[0]; 
				Efield[3*index+XPLUS] = Efield[3*neighbor[index][XMINUS]+XPLUS] + 
				(qx-qy)*maggs.prefactor*invasq;
			}
			
			if(ix>=lparam.in_ur[0]-1) {
				if(node_pos[0] < node_grid[0]-1) {
					if(node_grid[0]>1)
						MPI_Send(&Efield[3*index+XPLUS], 1, MPI_DOUBLE, node_neighbors[1], REQ_MAGGS_EQUIL, MPI_COMM_WORLD); 
				}
				else
					if (fabs(Efield[3*index+XPLUS]) > 10.*ROUND_ERROR_PREC)
						fprintf(stderr, "%d: Error in the calculation of Ex(%d,%d,%d)=%f!!\n",
								this_node, lattice[index].r[0], lattice[index].r[1], lattice[index].r[2],
								Efield[3*index+XPLUS]);	  
			}
		}    /*** loop over iy */
	}
	
	/** exchange halo-surfaces */
	maggs_exchange_surface_patch(Efield, 3, 1); 
	
	avgEz = 0.;
	for(iz=lparam.in_ld[2];iz<lparam.in_ur[2];iz++) {
		index = maggs_get_linear_index(lparam.in_ld[0], lparam.in_ld[1], iz, lparam.dim);
		avgEz += Efield[3*index+ZPLUS];
	}
	
	
	/*  MPI_Barrier(MPI_COMM_WORLD); */
	MPI_Allreduce(&avgEz,&gavgEz,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	gavgEz = gavgEz/(maggs.mesh*node_grid[0]*node_grid[1]);
	
	FORALL_INNER_SITES(ix, iy,iz) {
		index = maggs_get_linear_index(ix, iy, iz, lparam.dim);
		Efield[3*index+ZPLUS] -= gavgEz;
	}
	
	for(iz = lparam.in_ld[2];iz<lparam.in_ur[2];iz++) {
		avgEy = 0.;  
		for(iy = lparam.in_ld[1];iy<lparam.in_ur[1];iy++) {
			index = maggs_get_linear_index(lparam.in_ld[0], iy, iz, lparam.dim);
			avgEy += Efield[3*index+YPLUS];
		}    
		
		MPI_Allreduce(&avgEy, &gavgEy, 1, MPI_DOUBLE, MPI_SUM, zplane);
		gavgEy = gavgEy/(maggs.mesh*node_grid[0]);
		
		for(iy=lparam.in_ld[1];iy<lparam.in_ur[1];iy++) {
			for(ix=lparam.in_ld[0];ix<lparam.in_ur[0];ix++)  
				Efield[3*maggs_get_linear_index(ix, iy, iz, lparam.dim)+YPLUS] -= gavgEy;
		}
	}
	
	for(iz=lparam.in_ld[2];iz<lparam.in_ur[2];iz++) {
		for(iy=lparam.in_ld[1];iy<lparam.in_ur[1];iy++) {
			avgEx = 0.;
			for(ix=lparam.in_ld[0];ix<lparam.in_ur[0];ix++) {
				avgEx += Efield[3*maggs_get_linear_index(ix, iy, iz, lparam.dim)+XPLUS];
			}
			
			MPI_Allreduce(&avgEx, &gavgEx, 1, MPI_DOUBLE, MPI_SUM, yline);
			gavgEx = gavgEx/maggs.mesh;
			
			for(ix=lparam.in_ld[0];ix<lparam.in_ur[0];ix++)  
				Efield[3*maggs_get_linear_index(ix, iy, iz, lparam.dim)+XPLUS] -= gavgEx;
		}
	}
	
	
	/** exchange halo-surfaces */
	maggs_exchange_surface_patch(Efield, 3, 1);
	
	
	/* iterative procedure of energy minimization */
/*
	sqrE = 0.;
	FORALL_INNER_SITES(ix, iy, iz) {
		i = maggs_get_linear_index(ix, iy, iz, lparam.dim);
		FOR3D(k) sqrE += SQR(Efield[3*i+k]);
	}

	MPI_Allreduce(&sqrE,&gsqrE,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD); 
	gsqrE = gsqrE/(SPACE_DIM*maggs.mesh*maggs.mesh*maggs.mesh);  

	MAGGS_TRACE( if(!this_node) iteration = 0;);
	do {
		goldE = gsqrE;
		sqrE = 0.;
		maggs_minimize_transverse_field();
		
		FORALL_INNER_SITES(ix, iy, iz) {
			i = maggs_get_linear_index(ix, iy, iz, lparam.dim);
			FOR3D(k) sqrE += SQR(Efield[3*i+k]);
		}
		MPI_Allreduce(&sqrE,&gsqrE,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD); 
		gsqrE = gsqrE/(SPACE_DIM*maggs.mesh*maggs.mesh*maggs.mesh);  
		maxcurl = maggs_check_curl_E();
		
		MAGGS_TRACE(
					if(!this_node) {
						iteration++;
						if(iteration%1==0) {
							fprintf(stderr, "# iteration for field equilibration %d, diff=%9.4e, curlE=%9.4e\n", 
									iteration, fabs(gsqrE-goldE),maxcurl);
							fflush(stderr);    
						}
					}
					);
		
	} while(fabs(maxcurl)>1000000.*ROUND_ERROR_PREC);
	
*/
	
	/** exchange halo-surfaces */
	
	FOR3D(k) Eall[k] = 0.;
	FORALL_INNER_SITES(ix, iy, iz) {
		i = maggs_get_linear_index(ix, iy, iz, lparam.dim);
		FOR3D(k) {
			Eall[k] += Efield[3*i+k];
		}
	}
	
	MPI_Allreduce(Eall,gEall,3,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	
	FOR3D(k) gEall[k] /= (node_grid[0]*node_grid[1]*node_grid[2]*lparam.size[0]*lparam.size[1]*lparam.size[2]);
	
	FORALL_INNER_SITES(ix, iy, iz) {
		i = maggs_get_linear_index(ix, iy, iz, lparam.dim);
		FOR3D(k) Efield[3*i+k] -= gEall[k];
	}
	/** exchange whole glue-patch region */
	maggs_exchange_surface_patch(Efield, 3, 0);
	if(!this_node)
		MAGGS_TRACE(fprintf(stderr, "Ex = %16.12e, Ey = %15.12e, Ez = %15.12e\n", gEall[0], gEall[1], gEall[2]));
}


void maggs_calc_self_influence(Particle* P)
/** For each particle P, calculates self energy influence */
/** with direct Greens function, assuming constant permittivity */
/** on each lattice site. */
{
	int    k, ix, iy, iz;
	double invasq = maggs.inva*maggs.inva;
	double position[SPACE_DIM], left_down_position[SPACE_DIM], relative_position[SPACE_DIM];
	int index = 0;
	int globalindex = 0;
	double local_force[SPACE_DIM];
	double particle_charge = P->p.q;

	/* calculate position in cell, normalized to lattice size: */
	FOR3D(k) {
		position[k]           = (P->r.p[k] - lparam.ld_pos[k]) * maggs.inva;
		left_down_position[k] = floor(position[k]);
		relative_position[k]  = position[k] - left_down_position[k];
		local_force[k] = 0.0;
	}
	
	
	/* Copy permittivity values to the mini-lattice: */
	/*
	for (iz=0;iz<zmax;iz++) {
		for (iy=0;iy<ymax;iy++) {
			for (ix=0;ix<xmax;ix++) {
				index = (iz + zmax*(iy + ymax*ix));
				globalindex = maggs_get_linear_index((left_down_position[0]+ix),
													 (left_down_position[1]+iy),
													 (left_down_position[2]+iz), lparam.dim);
				self_influence_correction[index].permittivity = lattice[globalindex].permittivity;
				self_influence_correction[index].charge = 0.0; //background_charge;
				FOR3D(k) {
					charge_gradient[index][k] = 0.0;
					local_E_field[index][k] = 0.0;
				}
			}
		}
	}
	 */

				
	/* Calculate self-force directly: */	
	FOR3D(k) {
		local_force[k] = pow(((0.5 - relative_position[k])*1.57364595), 3.0) + // 1.57364595
		(0.5 - relative_position[k]) * 1.5078141; // 1.5078141
		local_force[k] *= maggs.bjerrum;
		local_force[k] *= particle_charge * particle_charge;
		local_force[k] *= invasq;
	}
	// maggs.prefactor = 3.544907702 = sqrt(4. * M_PI * maggs.bjerrum * temperature);
	// Fehlender Faktor: 4*pi=12.5663706
	// a^3/c = 2.58448064965803

	/* Correct force: */
	FOR3D(k) {
		P->f.f[k] -= local_force[k];
	}
}


void maggs_interpolate_charge(int *first, double *rel, double q)
/** Does the actual interpolating calculation. */
/** int *first: 3dim-array lattice position, */
/** double *rel: 3dim-array relative position in cube, */
/** double q: charge to interpolate */
{
	int i, k, l, m, index, temp_ind;
	int help_index[3];
	double temp;
	double help[SPACE_DIM];
	
	FOR3D(i) help[i] = 1. - rel[i];     /* relative pos. w.r.t. first */
	
	/* calculate charges at each vertex */
	index = maggs_get_linear_index(first[0],first[1],first[2],lparam.dim);
	
	FOR3D(i) {
		temp_ind = neighbor[index][i];
		if(temp_ind == NOWHERE) help_index[i] = lparam.volume; /* force huge index */
		else { /* incr. for x-neighbor */
			help_index[i] = maggs_get_offset(lattice[neighbor[index][i]].r[i], first[i], i, lparam.dim);
		}
		
	}
	
	for(k=0;k<2;k++){   /* jumps from x- to x+ */
		for(l=0;l<2;l++){  /* jumps from y- to y+ */
			for(m=0;m<2;m++){ /* jumps from z- to z+ */      
				if(index < lparam.volume) {
					temp = q;
					FOR3D(i) temp *= maggs_interpol1D(help[i]);
					lattice[index].charge += temp;
				}
				
				index+=help_index[2];
				help[2]=1.-help[2];
				help_index[2]=-help_index[2];
			}
			index+=help_index[1];
			help[1]=1.-help[1];
			help_index[1]=-help_index[1];
		}
		index+=help_index[0];
		help[0]=1.-help[0];
		help_index[0]=-help_index[0];
	}		
}


void maggs_calc_charge_gradients(double *rel, double q, double *grad)
/** Does the actual calculation of the gradient. Parameters: */
/** double *rel: 3dim-array of relative position in cube, */
/** double q: charge, */
/** double *grad: huge gradient array to write into */
{
	int i,l,m,index, d;
	double help[3];
	int dir1, dir2;
	
	FOR3D(i) help[i] = 1. - rel[i];     /* relative pos. w.r.t. x_int */
	
	index = 0;
	
	FOR3D(d) {
		maggs_calc_directions(d, &dir1, &dir2);
		for(l=0;l<2;l++){  /* jumps from dir2- to dir2+ */
			for(m=0;m<2;m++){ /* jumps from dir1- to dir1+ */          
				
				/* with q!!! */
				grad[index] = - q * help[dir1] * help[dir2];

				index++;
				help[dir1] = 1.-help[dir1];
			}
			help[dir2] = 1.-help[dir2];
		}
	}
}

void maggs_calc_charge_fluxes_1D(double q, double *help, double *flux, int dir)
/** calculate the charge flux in direction "dir". */
/** parameter "help" is the relative position in the cube */
/** to the opposite of left down front lattice site. */
{
	/** at the moment works only for linear interpolation */
	int index, dir1, dir2;
	int l,m; 
	double q_scaled;
	
	q_scaled = q * maggs.prefactor*maggs.inva;
	index = 0;
	
	maggs_calc_directions(dir, &dir1, &dir2);   
	
	for(l=0;l<2;l++){  /* jumps from dir2- to dir2+ */
		for(m=0;m<2;m++){ /* jumps from dir1- to dir1+ */   
			
			flux[index] = q_scaled * help[dir1]*help[dir2];
			index++;
			
			help[dir1] = 1. - help[dir1];
		}
		help[dir2] = 1. - help[dir2]; 
	}
}

short maggs_check_intersect_1D(double delta, double r_new, int dir, int first, double *t_step, int identity)
/** Extend particle trajectories on */
/** the one time step and check if */
/** the trajectory intersects the cell */
/** boundary in direction dir */
{	
	int candidateplane = -1; /* force alloc error */
	short f_crossing; 
	double r_old, temp;
	double ZERO = 0.0;
	double ONE  = 1.0;
	
	f_crossing = 0;
	r_old = r_new - delta;
	f_crossing = f_crossing||(r_old>=ONE||r_old<ZERO);
	
	if(dir==2) temp = 1.;
	else       temp = 0.5;
	
	if(f_crossing) {
		MAGGS_TRACE(
					fprintf(stderr, "Cube crossing in dir %d for particle %d at time = %f:\n", dir, identity, sim_time);
					fprintf(stderr,"  rold[%d]=%f, rnew[%d]=%f\n", dir, (first+r_old-1.)*maggs.a,
							dir, (first+r_new-1.)*maggs.a);
					fflush(stderr);
					);
		
		if(r_old >= ONE) candidateplane = ONE;
		if(r_old < ZERO) candidateplane = ZERO;
		
		/****** Update time step *********************/
		*t_step = temp * fabs((candidateplane-r_new)/delta);
	} /* end if crossing */
	else *t_step = temp;
	return f_crossing;
}

void maggs_add_current_on_segment(Particle *p, int ghost_cell)
{
	int d;
	int icoord, dir;
	int lat_index = -1; /* force alloc error */
	int first[SPACE_DIM];
	int f_crossing, flag_update_flux;
	t_dvector r_temp; 
	double inva_half;
	double delta;
	double flux[4], v;
	double v_inva[SPACE_DIM], v_invasq[SPACE_DIM];
	double pos[SPACE_DIM], help[3];
	double t_step;
	
	inva_half = 0.5*maggs.inva;
	
	FOR3D(d) {
		pos[d]   = (p->r.p[d] - lparam.ld_pos[d])* maggs.inva;
		first[d] = (int) floor(pos[d]);
		r_temp[d]   = pos[d] - first[d]; /* it is the updated coord (we have to go back) */
		help[d]     = 1. - r_temp[d];
		v_inva[d]   = maggs.inva * p->m.v[d];
		v_invasq[d] = maggs.inva * v_inva[d];
	}
	
	flag_update_flux = 1;
	if(ghost_cell) {
		FOR3D(d) if(first[d]<lparam.halo_ld[d] || first[d]>= lparam.halo_ur[d])
		{flag_update_flux = 0;break;}
	}
	
	if(flag_update_flux) {
		lat_index = maggs_get_linear_index(first[0], first[1], first[2], lparam.dim);
	}
	
	/* loop coordinates in order x->y->z->y->x */
	for(dir=0; dir<5; dir++) {
		icoord = dir;
		if(dir>2) icoord = dir%2;
		if(icoord == 2) delta = v_inva[icoord];
		else            delta = 0.5 * v_inva[icoord];
		
		f_crossing = maggs_check_intersect_1D(delta, r_temp[icoord], icoord, first[icoord], &t_step, p->p.identity);
		
		/** calculate flux */
		if(flag_update_flux) {
			maggs_calc_charge_fluxes_1D(p->p.q, help, flux, icoord);
			
			v = t_step * v_invasq[icoord];
			maggs_calc_e_field_on_link_1D(lat_index, flux, v, icoord);
		}
		
		if(f_crossing) {
			if(delta > 0.) {
				first[icoord]--;
				r_temp[icoord] += 1.;
			}
			else {
				first[icoord]++;
				r_temp[icoord] -= 1.;
			}
			if(icoord == 2) t_step = 1.  - t_step;
			else            t_step = 0.5 - t_step;
			
			if(ghost_cell){
				if(flag_update_flux) {
					if(first[icoord]<lparam.halo_ld[icoord] || first[icoord]>= lparam.halo_ur[icoord])
					{flag_update_flux = 0;}
				}
				else {
					flag_update_flux = 1;
					FOR3D(d) if(first[d]<lparam.halo_ld[d] || first[d]>= lparam.halo_ur[d])
					{flag_update_flux = 0;break;}
					if(flag_update_flux) maggs_calc_charge_fluxes_1D(p->p.q, help, flux, icoord); 
				}
			}
			
			if(flag_update_flux) {
				v = t_step * v_invasq[icoord];
				lat_index = maggs_get_linear_index(first[0], first[1], first[2], lparam.dim);
				maggs_calc_e_field_on_link_1D(lat_index, flux, v, icoord);
			}
		}
		r_temp[icoord] -= delta;
		help[icoord]    = 1. - r_temp[icoord];
	}
}


void symplect_couple_current_to_efield()
{
	/**************************************/
	/* Calculate fluxes and couple them   */
	/* with fields symplectically         */
	/* It is assumed that the particle    */
	/* can not cross more than one        */
	/* cell boundary per direction        */
	/**************************************/
	
	Cell *cell;
	Particle* p;
	int i, c, d, np;
	int flag_inner;
	double q;
	double r1, r2;
	
	/* loop over real particles */
	for (c = 0; c < local_cells.n; c++) {
		cell = local_cells.cell[c];
		p  = cell->part;
		np = cell->n;
		for(i = 0; i < np; i++) {
			if((q=p[i].p.q) != 0.) {
				/*	if(sim_time>49.08&&p[i].p.identity==231) */
				/*	  fprintf(stderr,"time=%f, v=(%f,%f,%f)\n",sim_time, p[i].m.v[0], p[i].m.v[1],p[i].m.v[2]); */
				maggs_add_current_on_segment(&p[i], 0);
			}/* if particle.q != ZERO */
		}
	}
	
	/* loop over ghost particles */
	for (c = 0; c < ghost_cells.n; c++) {
		cell = ghost_cells.cell[c];
		p  = cell->part;
		np = cell->n;
		for(i = 0; i < np; i++) {
			if((q=p[i].p.q) != 0.) {
				flag_inner = 1;
				FOR3D(d) {
					r2 = p[i].r.p[d];
					r1 = r2 - p[i].m.v[d];
					if(((r2 < lparam.ld_pos[d])&&(r1 < lparam.ld_pos[d]))
					   ||((r2 >= lparam.ur_pos[d] && r1 >= lparam.ur_pos[d])))
					{flag_inner = 0; break;}
				}
				if(flag_inner) {
					maggs_add_current_on_segment(&p[i], 1);
				}
			}/* if particle.q != ZERO */
		}
	}
}

void  maggs_calc_e_field_on_link_1D(int index, double *flux, double v, int dir)
/* updates field on link coupling it with current.
 * Force is multiplied by the time_step
 */
{  
	int l, m, ind_flux, dir1, dir2;
	int temp_ind;
	int help_index[2];
	int* anchor_neighb;
	t_site* anchor_site;
	
	maggs_calc_directions(dir, &dir1, &dir2);
	
	anchor_neighb = &neighbor[index][0]; 
	anchor_site = &lattice[index];
	
	temp_ind = anchor_neighb[dir1];
	if(temp_ind == NOWHERE) help_index[0] = lparam.volume;
	else
		help_index[0] = maggs_get_offset(lattice[temp_ind].r[dir1], anchor_site->r[dir1], dir1, lparam.dim);
	temp_ind = anchor_neighb[dir2];
	if(temp_ind == NOWHERE) help_index[1] = lparam.volume;
	else
		help_index[1] = maggs_get_offset(lattice[temp_ind].r[dir2], anchor_site->r[dir2], dir2, lparam.dim);
	
	
	ind_flux = 0;
	for(l=0;l<2;l++){  /* jumps from dir2- to dir2+ */
		for(m=0;m<2;m++){ /* jumps from dir1- to dir1+ */  
			
			if(index < lparam.volume){
				Efield[3*index+dir] -= flux[ind_flux] * v;
			}
			
			ind_flux++; 
			
			index+=help_index[0];
			help_index[0]=-help_index[0];	
		}
		index+=help_index[1];
		help_index[1]=-help_index[1];     
	}
}

void calc_part_link_forces(Particle *p, int index, double *grad)
/* Calculate the actual force from the E-Field */
/* This is, where the actual INTERLACING happens! */
{
	static int init = 1;
	static int help_index[SPACE_DIM];
	int ind_grad, j;
	int dir1, dir2;
	/*  int* anchor_neighb; */
	int l,m;
	double local_force[SPACE_DIM];
	
	if(init) {
		t_site* anchor_site;
		anchor_site = &lattice[index];
		FOR3D(j) {
			help_index[j] = maggs_get_offset(lattice[neighbor[index][j]].r[j], anchor_site->r[j], j, lparam.dim);
		}
		init = 0;
	}
	
	FOR3D(j){
		local_force[j] = 0.;
	}
	
	ind_grad = 0; 
	
	FOR3D(j) {
		maggs_calc_directions(j, &dir1, &dir2);
		
		for(l=0;l<2;l++){  /* jumps from dir2- to dir2+ */
			for(m=0;m<2;m++){ /* jumps from dir1- to dir1+ */   
				local_force[j] += -grad[ind_grad]*Efield[3*index+j];
				//fprintf(stderr, "charge_gradient %d: %1.9f\n", ind_grad, grad[ind_grad]);
				ind_grad++;
				index += help_index[dir1];
				help_index[dir1] = -help_index[dir1];
			}
			index += help_index[dir2];
			help_index[dir2] = -help_index[dir2];
		}
	}  
	/* Attention! Here, the INTERLACING is done! */
	FOR3D(j){
		p->f.f[j] += maggs.prefactor * local_force[j];
	}
}

void add_transverse_field_to_e_field(double dt)
{
	int i, index;
	double invasq; 
	int x, y, z;
	int offset, xoffset, yoffset;
	double help;
	
	invasq = SQR(maggs.inva);
	help = dt * invasq * maggs.invsqrt_f_mass;
	
	/***calculate e-field***/ 
	offset = maggs_get_linear_index(1,1,1, lparam.dim);
	yoffset = lparam.dim[2];
	xoffset = 2*lparam.dim[2];
	for(x=0;x<lparam.size[0];x++) {
		for(y=0;y<lparam.size[1];y++) {
			for(z=0;z<lparam.size[2];z++) {
				/*  FORALL_INNER_SITES(x, y, z) { */
				/*    i = maggs_get_linear_index(x, y, z, lparam.dim); */
				i = offset+z;
				index = 3*i;
				Efield[index  ] += help * maggs_calc_curl(2, 1, Bfield, neighbor[i], index);
				Efield[index+1] += help * maggs_calc_curl(0, 2, Bfield, neighbor[i], index);
				Efield[index+2] += help * maggs_calc_curl(1, 0, Bfield, neighbor[i], index);
			}
			offset += yoffset;
		}
		offset += xoffset;
	} 
	
	maggs_exchange_surface_patch(Efield, 3, 0);
}

void maggs_calc_forces()
{ 
	Cell *cell;
	static int init = 1;
	static int Npart_old;
	Particle *p;
	int i, c, np, d, index, Npart, ip; 
	double q;
	/* position of a particle in local lattice units */
	double pos[SPACE_DIM];
	/* index of first assignment lattice point */
	int first[3];
	/* charge gradient (number of neighbor sites X number of dimensions) */
	static double *grad;
	
	if(init) Npart_old = 0;
	
	Npart = cells_get_n_particles();
	if(Npart>Npart_old) {
		grad = (double *) realloc(grad, 12*Npart*sizeof(double));
		Npart_old = Npart;
	}
	
	/* Hopefully only needed for Yukawa: */
	maggs_calc_grad_and_point_forces(grad);
	
	if(!init) {
		MAGGS_TRACE(fprintf(stderr, "running symplectic update\n"));
		symplect_couple_current_to_efield();
		add_transverse_field_to_e_field(time_step);  
	}
	else init = 0;
	
	ip = 0;
	for (c = 0; c < local_cells.n; c++) {
		cell = local_cells.cell[c];
		p  = cell->part;
		np = cell->n;
		for(i=0; i<np; i++) { 
			q = p[i].p.q;
			if( abs(q) > 1.0e-5 ) {
				FOR3D(d) {
					pos[d]   = (p[i].r.p[d] - lparam.ld_pos[d])* maggs.inva;
					first[d] = (int) pos[d];
				}
				
				index = maggs_get_linear_index(first[0],first[1],first[2],lparam.dim);
				calc_part_link_forces(&p[i], index, &grad[ip]);
				maggs_calc_self_influence(&p[i]);

				ip+=12;
			}
		}
	}
}


int maggs_set_parameters(Tcl_Interp *interp, double bjerrum, double f_mass, int mesh)
{
	if (f_mass <=0.) {
		Tcl_AppendResult(interp, "mass of the field is negative", (char *)NULL);
		return TCL_ERROR;
	} 
	if(mesh<0) {
		Tcl_AppendResult(interp, "mesh must be positive", (char *) NULL);
		return TCL_ERROR;
	}
	
	maggs.mesh           = mesh; 
	maggs.bjerrum        = bjerrum;
	maggs.f_mass         = f_mass; 
	maggs.invsqrt_f_mass = 1./sqrt(f_mass); 	
	
	mpi_bcast_coulomb_params();
	
	return TCL_OK;
}


void Maggs_init()
{

	int d, ncharges;
	int max_node_grid;

	if(maggs.bjerrum == 0.) {
		if(this_node==0) 
			MAGGS_TRACE(fprintf(stderr,"0: Maggs_init: Bjerrum length is zero.\n");
						fprintf(stderr,"   Electrostatics switched off!\n"));    
	}
	else {
		MAGGS_TRACE(fprintf(stderr,"\n%d: Maggs_init: \n",this_node));
		if( (box_l[0] != box_l[1]) || (box_l[1] != box_l[2]) ) {
			if(this_node==0) {
				fprintf(stderr,"0: Maggs_init: SERIOUS WARNING:\n"); 
				fprintf(stderr,"   No long range interactions for non cubic box.\n"); 
				fprintf(stderr,"   Switch off long range interactions! \n");
			}
			maggs.bjerrum =  0.0;
			return;
		} 
		
		max_node_grid = 1;
		FOR3D(d) if(node_grid[d] > max_node_grid) 
			max_node_grid = node_grid[d];
		
		if(maggs.mesh%max_node_grid != 0) {
			if(this_node==0) {
				fprintf(stderr,"%d: Maggs_init: SERIOUS WARNING:\n", this_node); 
				fprintf(stderr,"   Number of mesh points is incompatible with number of processes.\n"); 
				errexit();
			}
		}
		
		ncharges = maggs_count_charged_particles();
		if(this_node == 0 && ncharges == 0) {
			fprintf(stderr,"\nFailed to initialize fields - no charged particles in the system!\n");
			errexit();
		}
		
		maggs.inva  = (double) maggs.mesh/box_l[0]; 
		maggs.a     = 1.0/maggs.inva;
		if (temperature>0.0) {
			maggs.prefactor  = sqrt(4. * M_PI * maggs.bjerrum * temperature);
			maggs.pref2      = maggs.bjerrum * temperature;
		} else {
			maggs.prefactor  = sqrt(4. * M_PI * maggs.bjerrum);
			maggs.pref2      = maggs.bjerrum;
		}			
		
		
		maggs_setup_local_lattice();
		
		/* update max_cut */
		integrate_vv_recalc_maxrange();
		on_parameter_change(FIELD_MAXRANGE);
		/* enforce electric field onto the Born-Oppenheimer surface */
		maggs_calc_init_e_field();
		//if(!this_node) fprintf(stderr, "%d: Electric field is initialized\n", this_node);
		MAGGS_TRACE(fprintf(stderr,"%d: maggs initialized\n",this_node));
	}
}

void maggs_propagate_B_field(double dt) {
	int x, y, z, i, offset, index;
	int xoffset, yoffset;
	double help = dt*maggs.invsqrt_f_mass;
	/* B(t+h/2) = B(t-h/2) + h*curlE(t) */ 
	
	offset = maggs_get_linear_index(1,1,1, lparam.dim);
	yoffset = lparam.dim[2];
	xoffset = 2*lparam.dim[2];
	
	for(x=0;x<lparam.size[0];x++) {
		for(y=0;y<lparam.size[1];y++) {
			for(z=0;z<lparam.size[2];z++) {

				i = offset+z;
				index = 3*i;
				Bfield[index+0] += - help*maggs_calc_dual_curl(1,2, Efield, neighbor[i], index); 
				Bfield[index+1] += - help*maggs_calc_dual_curl(2,0, Efield, neighbor[i], index); 
				Bfield[index+2] += - help*maggs_calc_dual_curl(0,1, Efield, neighbor[i], index);  
			}
			offset += yoffset;
		}
		offset += xoffset;
	}
	
	maggs_exchange_surface_patch(Bfield, 3, 0);
}


/************************************************************/
int tclprint_to_result_Maggs(Tcl_Interp *interp)
{
	char buffer[TCL_DOUBLE_SPACE];
	
	Tcl_PrintDouble(interp, maggs.f_mass, buffer);
	Tcl_AppendResult(interp, "maggs ", buffer, " ", (char *) NULL);
	sprintf(buffer,"%d",maggs.mesh);
	Tcl_AppendResult(interp, buffer, " ", (char *) NULL); 
	
	return TCL_OK;
}

int tclcommand_inter_coulomb_parse_maggs(Tcl_Interp * interp, int argc, char ** argv)
{
	int mesh;
	double f_mass;

	
	if(argc < 2) {
		Tcl_AppendResult(interp, "Not enough parameters: inter coulomb maggs <f_mass> <mesh>", (char *) NULL);
		return TCL_ERROR;
	}
	
	if(! ARG_IS_D(0, f_mass))
		return TCL_ERROR;
	
	if(! ARG_IS_I(1, mesh)) {
		Tcl_AppendResult(interp, "integer expected", (char *) NULL);
		return TCL_ERROR;
	}
		
	if(argc > 2) {
		Tcl_AppendResult(interp, "Too many parameters: inter coulomb maggs <f_mass> <mesh>", (char *) NULL);
		return TCL_ERROR;
	}
	
	coulomb.method = COULOMB_MAGGS;
	
	return maggs_set_parameters(interp, coulomb.bjerrum, f_mass, mesh);
}


double maggs_electric_energy()
/** integrates 0.5*D*E over the whole system */
/** public function! */
{
	int x, y, z, i, k;
	double localresult = 0.;
	double globalresult = 0.;
	
	FORALL_INNER_SITES(x, y, z) {
		i = maggs_get_linear_index(x, y, z, lparam.dim);	  
		FOR3D(k){
			localresult += SQR(Efield[i*3+k]) / lattice[i].permittivity[k];
		}
	}
	localresult *= 0.5*maggs.a;
	MPI_Allreduce(&localresult,&globalresult,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);  
	return globalresult;
}


double maggs_magnetic_energy()
{
	int x, y, z, i;
	double result = 0.;
	/*  double invmass = 1./maggs.f_mass; we have B^~=B*c !!!! */
	
	FORALL_INNER_SITES(x, y, z) {
		i = maggs_get_linear_index(x, y, z, lparam.dim);	  
		result += SQR(Bfield[i*3]) + SQR(Bfield[i*3+1]) + SQR(Bfield[i*3+2]);
	}
	/* B is rescaled !!! ATTENTION!!! */
	result *= 0.5*maggs.a;
	return result;
}


void Maggs_exit()
{
	/*  free(send_databuf); */
	/*  free(recv_databuf); */
	free(lattice);
	free(Efield);
	free(Bfield);
}
#endif

