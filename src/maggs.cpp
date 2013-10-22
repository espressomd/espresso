/*
  Copyright (C) 2010,2011,2012,2013 The ESPResSo project
  Copyright (C) 2010,2011 Florian Fahrenberger
  Copyright (C) 2007,2008,2009,2010 
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

/** \file maggs.cpp
 *  Maxwell Equations Molecular Dynamics (MEMD) method for electrostatic
 *  interactions.
 *
 *  We use a local update scheme to propagate artificial B-fields on a
 *  lattice in the system. In principal, the algorithm simulates full
 *  electrodynamics, but with a tunable speed of light.
 *
 *  The method is very usable for large particle numbers or highly
 *  parallel architectures, since it is local and scales linearly.
 *  It is not suited for high-precision calculation of forces, since
 *  the simple interpolation scheme produces errors in the order of
 *  10^-5 in the force.
 *
 *  The chosen mesh should roughly be of the size of the particles.
 *
 *  Further reading on the algorithm:
 *  I. Pasichnyk and B. Dunweg, Coulomb interaction via local
 *  dynamics: a molecular-dynamics algorithm. J. Phys:
 *  Condens. Matter, 16 ,p. 1399-4020, (2004).
 *  
 */


#include <mpi.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "ghosts.hpp"
#include "global.hpp"
#include "grid.hpp"
#include "integrate.hpp"
#include "initialize.hpp"
#include "interaction_data.hpp"
#include "particle_data.hpp"
#include "communication.hpp"
#include "maggs.hpp"
#include "thermostat.hpp"
#include "cells.hpp"
#include "domain_decomposition.hpp"
#include "errorhandling.hpp"

#ifdef ELECTROSTATICS

/* MPI tags for the maggs communications: */
/* Used in maggs_init() -> calc_glue_patch(). */
#define REQ_MAGGS_SPREAD 300
#define REQ_MAGGS_EQUIL  301

/* Factors for self-influence currection */
/* (stemming from epsilon_zero and 4*pi) */
#define SELF_FACTOR_1 1.57364595
#define SELF_FACTOR_2 1.5078141
// from lattice Green's function: 0.5054620197
// and 0.126365505
// e_0 = 8.85418781762 * 10^-12

/* Define numbers for directions and dimensions: */
#define SPACE_DIM 3                 /* number of dimensions */
#define NOWHERE -1                  /* not a direction */
#define NDIRS 6				        /* number of directions */
#define XPLUS 0                     /* add for specific direction */
#define YPLUS 1
#define ZPLUS 2
#define ZMINUS 3
#define YMINUS 4
#define XMINUS 5
#define OPP_DIR(dir)	(5-(dir))	/* Opposite direction */


/* Three often used macros for looping over 3D */
#define FOR3D(dir) for(dir=0; dir<SPACE_DIM; dir++)

#define FORALL_INNER_SITES(i,j,k)					\
  for(i=lparams.inner_left_down[0];i<lparams.inner_up_right[0];i++)	\
    for(j=lparams.inner_left_down[1];j<lparams.inner_up_right[1];j++)	\
      for(k=lparams.inner_left_down[2];k<lparams.inner_up_right[2];k++) 

#define FORALL_SITES(i,j,k)			\
  for(i=0;i<lparams.dim[0];i++)			\
    for(j=0;j<lparams.dim[1];j++)		\
      for(k=0;k<lparams.dim[2];k++) 
/* from ifndef MAGGS_H */

/***************************************/
/****** data types and structures ******/
/***************************************/

typedef int      t_ivector [SPACE_DIM]; /* integer vector for position in grid */
typedef double   t_dvector [SPACE_DIM]; /* double vector for fields etc. */
typedef int      t_dirs    [NDIRS];     /* integer vector for directions */

/** Structure of local lattice parameters. */
typedef struct {
  t_dvector left_down_position;        /* spatial position of left down grid point */  
  t_dvector upper_right_position;      /* spatial positon of upper right grid point */ 
  int inner_left_down[3];              /* inner left down grid point    */
  int inner_up_right[3];               /* inner up right grid point + (1,1,1) */
  int halo_left_down[3];               /* halo-region left down grid point  */
  int halo_upper_right[3];             /* halo-region up right global grid point  */
  int margin[SPACE_DIM*2];             /* number of margin mesh points (even index - left, odd - right). */
  t_ivector dim;                       /* grid dimension (size + glue_patch region) of local mesh.  */
  t_ivector size;                      /* dimension of mesh inside node domain. */
  int volume;                          /* number of lattice sites in local domain */
} lattice_parameters;

/** surface_patch structure for communication. */
typedef struct {
  int offset;        /* source offset for the site index */
  int doffset;       /* destination offset for the site index */
  int stride;        /* minimal contiguous block */  
  int skip;          /* gap between two strides (from the first element of one stride to the first elem. of next stride */
  int nblocks;       /* number of strides */
  int coord[2];      /* coordinates of the vector fields which has to be exchanged */
  int volume;        /* number of lattice sites in surface patch */
} t_surf_patch;

/** structure for properties of each lattice site */
typedef struct {
  double    charge;                     /* charge on site */
  double    permittivity[SPACE_DIM];    /* dielectric properties on adjoined bonds. */
  int       r[SPACE_DIM];               /* position of site in space */
} t_site;






/*******************************/
/****** global variables ***** */
/*******************************/

/* create system structure with all zeros. Filled in maggs_set_parameters(); */
MAGGS_struct maggs = { 1, 1.0, 0. , 0. , 0. , 0. , 0. , 0 , 0. , 0. , {{0.},{0.}}};

/* local mesh. */
static lattice_parameters lparams;
/* local lattice */
static t_site* lattice;
/* local D field */
static double* Dfield;
/* local B field */
static double* Bfield;
/* site neighbors */
static t_dirs* neighbor;






/*******************************************************/
/*      Private Functions list with short comment      */
/*******************************************************/
/*
 * This is supposed to give an overview on the contained functions.
 * Declaration is not needed, since the functions are in right order,
 * but it is possible to do so here, if something is added with
 * non-fitting cross-references.
 *
 * Longer descriptions of the functions' tasks can be found in the
 * implementation.
 */

/****** small helper functions: ******/
// int maggs_get_linear_index(int x, int y, int z, int latticedim[SPACE_DIM]); /* linear indexing for speed */
// int maggs_count_charged_particles(); /* global number of charges */
// int maggs_get_offset(int index_shift, int index_base, int axes, int adim[3]); /* offset of currect lattice site */
// void maggs_calc_directions(int j, int* dir1, int*dir2); /* for a given direction, give back the other two */
// double maggs_calc_curl(int mue, int nue, double* field, int* Neighbor, int index); /* calculate curl in real space */
// double maggs_calc_dual_curl(int mue, int nue, double* field, int* Neighbor, int index); /* calculate curl in dual space */
// void maggs_update_plaquette(int mue, int nue, int* Neighbor, int index, double delta); /* update fields on plaquette */

/****** setup everything: ******/
// int maggs_set_parameters(double bjerrum, double f_mass, int mesh); /* set the system parameters */
// void maggs_setup_neighbors(); /* setup the site neighbors */
// void maggs_setup_local_lattice(); /* set lattice parameters */
// void maggs_calc_surface_patches(t_surf_patch* surface_patch); /* prepare for communication */
// void maggs_prepare_surface_planes(int dim, MPI_Datatype *xy, MPI_Datatype *xz, MPI_Datatype *yz, t_surf_patch *surface_patch); /* prepare for communication */

/****** communication function: ******/
// void maggs_exchange_surface_patch(double *field, int dim, int e_equil); /* communicate */

/****** interpolate charges on lattice: ******/
// double maggs_interpol1D(double x); /* interpolate charge linearly in one direction */
// void maggs_interpolate_charge(int *first, double *rel, double q); /* interpolate all charges */
// void maggs_accumulate_charge_from_ghosts(); /* interplate ghost charges */
// void maggs_distribute_particle_charges(); /* find cube and call interpolation */
// void maggs_calc_charge_gradients(double *rel, double q, double *grad); /* calculate gradients */
// void maggs_update_charge_gradients(double *grad); /* find cube and call calculate gradients */

/****** initialization procedure: ******/
// double maggs_check_curl_E(); /* check Gauss law (maximum deviation) */
// void maggs_perform_rot_move_inplane(int i, int n); /* calculate correction for plaquette */
// void maggs_minimize_transverse_field(); /* B field minimization routine */
// void maggs_calc_init_e_field(); /* initial solution */

/****** calculate currents and E-fields ******/
// void maggs_calc_charge_fluxes_1D(double q, double *help, double *flux, int dir); /* charge flux for current */
// int maggs_check_intersect_1D(double delta, double r_new, int dir, int first, double *t_step, int identity); /* check if particles change cubes */
// void maggs_calc_e_field_on_link_1D(int index, double *flux, double v, int dir); /* update field from current */
// void maggs_add_current_on_segment(Particle *p, int ghost_cell); /* loop, find cube, add current */
// void maggs_couple_current_to_Dfield() /* update field from current */

/****** calculate B-fields and forces ******/
// void maggs_propagate_B_field(double dt); /* propagate the B-field */
// void maggs_add_transverse_field(double dt) /* calculate E-field from B-field */
// void maggs_calc_self_influence(Particle* P); /* correct self influence */
// void maggs_calc_part_link_forces(Particle *p, int index, double *grad) /* calculate force from D-field */
// void maggs_calc_forces(); /* calculate all forces correctly */

/****** get energy and print out stuff ******/
// double maggs_electric_energy(); /* calculate electric (E-field) energy */
// double maggs_magnetic_energy(); /* calculate magnetic (B-field) energy */

/****** init and exit ******/
// maggs_init(); /* initialize: check parameters, assign memory, calculate initial field */
// maggs_exit(); /* free memory */







/*************************************/
/****** small helper functions: ******/
/*************************************/

/** Turns the 3D index into a linear index.
    adim contains the dimensions of the lattice in
    the 3 directions: adim[0] = x_max, ..., adim[2] = z_max
    z is the first index!!!
    @return linear index
    @param x           index in x direction
    @param y           index in y direction
    @param z           index in z direction
    @param latticedim  dimensions of the lattice
*/
int maggs_get_linear_index(int x, int y, int z, int latticedim[SPACE_DIM])
{
  return (z + latticedim[ZPLUS]*(y + latticedim[YPLUS]*x));
}

/** Counts the total number of charged particles on 
    all processors.
    @return total number of charged particles in the system
*/
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
	
  MPI_Reduce(&node_sum, &tot_sum, 1, MPI_DOUBLE, MPI_SUM, 0, comm_cart);
	
  return tot_sum;
}


/** Index shift is calculated for moving in various
    directions with the linear index.
    @return Offset of current node from base
    @param index_shift     amount to move in direction
    @param index_base      index of base on current processor
    @param axes            in which direction to move
    @param adim            dimensions of the local lattice
*/
int maggs_get_offset(int index_shift, int index_base, int axes, int adim[3])
{
  int dif;
  dif = index_shift - index_base;
  if(axes <= 1) dif *= adim[2];
  if(axes == 0) dif *= adim[1]; 
	
  return (dif);
}


/** For any direction j, write the two other directions
    into the pointers in circular permutation.
    @param j    given first direction
    @param dir1 write second direction into
    @param dir2 write third direction into
*/
void maggs_calc_directions(int j, int* dir1, int*dir2)
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

/** Calculates the finite differences rotation in real space in mue-nue plane:
    \f$\frac{\partial}{\partial t}{D} = \nabla \times B\f$ (and prefactors plus current)
    The given "double* field" should be a B-Field!!
    @return rotation result
    @param mue       direction 1 of plane to rotate in
    @param nue       direction 2 of plane to rotate in
    @param field     input B-field
    @param Neighbor  neighbor lattice site
    @param index     index of current lattice site
*/
double maggs_calc_curl(int mue, int nue, double* field, int* Neighbor, int index)
{
  double result;
	
  result = field[index+mue] + field[3*Neighbor[OPP_DIR(mue)]+nue] -
    field[3*Neighbor[OPP_DIR(nue)]+mue] - field[index+nue];

  return result;
}

/** Calculates the finite differences rotation in dual space in mue-nue plane:
    \f$\frac{\partial}{\partial t}{B} = - \nabla \times D / (\epsilon)\f$
    The given "double* field" should be a D-Field!!
    @return rotation result
    @param mue       direction 1 of plane to rotate in
    @param nue       direction 2 of plane to rotate in
    @param field     input D-field
    @param Neighbor  neighbor lattice site
    @param index     index of current lattice site
*/
double maggs_calc_dual_curl(int mue, int nue, double* field, int* Neighbor, int index)
{
  double res;
	
  res = field[index+mue] + field[3*Neighbor[mue]+nue] -
    field[3*Neighbor[nue]+mue] - field[index+nue];
	
  return res;
}


/** updates all D-fields on links of the plaquette
    and the surroundings. delta was calculated before
    in function "maggs_perform_rot_move_inplane".
    @param mue        direction 1 of update
    @param nue        direction 2 of update
    @param Neighbor   neighbor lattice site
    @param index      index of current lattice site
    @param delta      by which amount to update field
*/
void maggs_update_plaquette(int mue, int nue, int* Neighbor, int index, double delta)
{
  int i = 3*index;
  Dfield[i+mue]             += delta;
  Dfield[3*Neighbor[mue]+nue] += delta;
  Dfield[3*Neighbor[nue]+mue] -= delta;
  Dfield[i+nue]             -= delta;  
}


/** Basic sanity checks to see if the code will run.
    @return zero if everything is fine. -1 otherwise.
*/
int maggs_sanity_checks()
{
  char *errtxt;
  int ret = 0;
  int d;
  int max_node_grid = 1;

  FOR3D(d) if(node_grid[d] > max_node_grid) max_node_grid = node_grid[d];
	
  if (maggs.bjerrum == 0.) {
    errtxt = runtime_error(128);
    ERROR_SPRINTF(errtxt, "{301 MEMD: bjerrum length is zero.} ");
    ret = -1;
  }
  else if ( (box_l[0] != box_l[1]) || (box_l[1] != box_l[2]) ) {
    errtxt = runtime_error(128);
    ERROR_SPRINTF(errtxt, "{302 MEMD needs cubic box.} ");
    ret = -1;
  }
  if (!PERIODIC(0) || !PERIODIC(1) || !PERIODIC(2)) {
    errtxt = runtime_error(128);
    ERROR_SPRINTF(errtxt, "{303 MEMD requires periodicity 1 1 1} ");
    ret = 1;
  }
  else if ( maggs.mesh%max_node_grid != 0 ) {
    errtxt = runtime_error(128);
    ERROR_SPRINTF(errtxt, "{304 MEMD: meshsize is incompatible with number of processes.} ");
    ret = -1;
  }
  /*
  else if ( maggs_count_charged_particles() == 0 ) {
      errtxt = runtime_error(128);
      ERROR_SPRINTF(errtxt, "{30? MEMD: No charges in the system.} ");
      ret = -1;
  }
  */
  else if (cell_structure.type != CELL_STRUCTURE_DOMDEC) {
    errtxt = runtime_error(128);
    ERROR_SPRINTF(errtxt, "{305 MEMD requires domain-decomposition cellsystem.} ");
    ret = -1;
  }
  else if (dd.use_vList) {
    errtxt = runtime_error(128);
    ERROR_SPRINTF(errtxt, "{306 MEMD requires no Verlet Lists.} ");
    ret = -1;
  }
  /** check if speed of light parameter makes sense */
  else if (maggs.f_mass < ( 2. * time_step * time_step / maggs.a / maggs.a ) ) {
    errtxt = runtime_error(128);
    ERROR_SPRINTF(errtxt, "{307 MEMD: Speed of light is set too high. Increase f_mass.} ");
    ret = -1;      
  }
  else if (maggs.a < skin) {
    errtxt = runtime_error(128);
    ERROR_SPRINTF(errtxt, "{308 MEMD: Skin should be smaller than MEMD mesh size.} ");
    ret = -1;
  }
#ifdef EXTERNAL_FORCES
  /** check for fixed particles */
  for (int cellnumber = 0; cellnumber < local_cells.n; cellnumber++) {
    Cell *cell = local_cells.cell[cellnumber];
    Particle *p  = cell->part;
    int np = cell->n;
    for(int i = 0; i < np; i++) {
      if ( (p[i].p.q != 0.0) & p[i].l.ext_flag & COORDS_FIX_MASK) {
	errtxt = runtime_error(128);
	ERROR_SPRINTF(errtxt, "{309 MEMD does not work with fixed particles.} ");
	ret = -1;
      }
    }
  }  
#endif
  
  return ret;
}

void maggs_compute_dipole_correction()
{
/*
maggs.prefactor  = sqrt(4. * M_PI * maggs.bjerrum * temperature);
                 = sqrt(1/epsilon)
maggs.pref2      = maggs.bjerrum * temperature;
                 = 1/4*pi*epsilon
*/
 
    /* Local dipole moment */
    double local_dipole_moment[3] = {0.0, 0.0, 0.0};
    /* Global dipole moment */
    double dipole_moment[3];
    
    int dim,c,np,i;
    Particle* p;
    Cell* cell;

    /* Compute the global dipole moment */
    for (c = 0; c < local_cells.n; c++) {
        cell = local_cells.cell[c];
        p  = cell->part;
        np = cell->n;
        for(i=0; i<np; i++)
            for (dim=0;dim<3;dim++)
                 local_dipole_moment[dim] += p[i].r.p[dim] * p[i].p.q;
    }

    MPI_Allreduce(local_dipole_moment, dipole_moment, 3, MPI_DOUBLE, MPI_SUM, comm_cart);
    
    double volume = box_l[0] * box_l[1] * box_l[2];
    double dipole_prefactor = 4.0*M_PI / (2.0*volume*(maggs.epsilon_infty + 1.0));

    /* apply correction to all particles: */
    for (c = 0; c < local_cells.n; c++) {
        cell = local_cells.cell[c];
        p  = cell->part;
        np = cell->n;
        for(i=0; i<np; i++)
            for (dim=0;dim<3;dim++)
                p[i].f.f[dim] += p[i].p.q * dipole_prefactor * dipole_moment[dim];
    }
}





/*******************************/
/****** setup everything: ******/
/*******************************/

int maggs_set_parameters(double bjerrum, double f_mass, int mesh, int finite_epsilon_flag, double epsilon_infty)
{
  if (f_mass <=0.) {
    return -1;
  } 
  if(mesh<0) {
    return -2;
  }
	
  maggs.mesh           = mesh;
  maggs.finite_epsilon_flag  = finite_epsilon_flag;
  maggs.epsilon_infty  = epsilon_infty;
  maggs.bjerrum        = bjerrum;
  maggs.f_mass         = f_mass; 
  maggs.invsqrt_f_mass = 1./sqrt(f_mass); 	
	
  mpi_bcast_coulomb_params();
	
  return ES_OK;
}


/** sets up nearest neighbors for each site in linear index */
void maggs_setup_neighbors()
{
  int ix = 0;
  int iy = 0;
  int iz = 0;
	
  int xsize = lparams.dim[0];
  int ysize = lparams.dim[1];
  int zsize = lparams.dim[2];
	
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
				
	      kount         = maggs_get_linear_index(ix,      iy,      iz,      lparams.dim);
	      kountzplus    = maggs_get_linear_index(ix,      iy,      izplus,  lparams.dim);
	      kountzminus   = maggs_get_linear_index(ix,      iy,      izminus, lparams.dim);
	      kountyplus    = maggs_get_linear_index(ix,      iyplus,  iz,      lparams.dim);
	      kountyminus   = maggs_get_linear_index(ix,      iyminus, iz,      lparams.dim);
	      kountxplus    = maggs_get_linear_index(ixplus,  iy,      iz,      lparams.dim);
	      kountxminus   = maggs_get_linear_index(ixminus, iy,      iz,      lparams.dim);
				
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


/** Set up lattice, calculate dimensions and lattice parameters
    Allocate memory for lattice sites and fields */
void maggs_setup_local_lattice()
{
  int i;
  int ix = 0;
  int iy = 0;
  int iz = 0;
  int linearindex = 0;
  int xyzcube;	
	
  xyzcube = 1;
  FOR3D(i) {
    /** inner left down grid point (global index) */
    lparams.inner_left_down[i] = (int)ceil(my_left[i]*maggs.inva); 
    /** inner up right grid point (global index) */
    lparams.inner_up_right[i] = (int)floor(my_right[i]*maggs.inva); 
    /** correct roundof errors at boundary */
    if(my_right[i]*maggs.inva-lparams.inner_up_right[i]<ROUND_ERROR_PREC) lparams.inner_up_right[i]--;
    if(1.0+my_left[i]*maggs.inva-lparams.inner_left_down[i]<ROUND_ERROR_PREC) lparams.inner_left_down[i]--;
    /** inner grid dimensions */
    lparams.size[i] = lparams.inner_up_right[i] - lparams.inner_left_down[i] + 1;
    /** spacial position of left down grid point */
    lparams.left_down_position[i] = my_left[i] - maggs.a;  
    /** spacial position of upper right grid point */
    lparams.upper_right_position[i] = my_right[i] + maggs.a;  
    /** left down margin */
    lparams.margin[i*2] = 1;
    /** up right margin */
    lparams.margin[(i*2)+1] = 1;
		
    lparams.dim[i] = lparams.size[i] + lparams.margin[i*2] + lparams.margin[i*2+1];
    xyzcube *= lparams.dim[i];
    /** reduce inner grid indices from global to local */
    lparams.inner_left_down[i] = lparams.margin[i*2];
    lparams.inner_up_right[i] = lparams.margin[i*2]+lparams.size[i];
    lparams.halo_left_down[i] = 0;
    lparams.halo_upper_right[i] = lparams.inner_up_right[i];      
  }
	
	
  lparams.volume    = xyzcube;
  /** allocate memory for sites and neighbors */
  lattice  = (t_site*) malloc(xyzcube*sizeof(t_site));
  neighbor = (t_dirs*) malloc(xyzcube*sizeof(t_dirs));
	
  Bfield   = (double*) malloc(3*xyzcube*sizeof(double));
  Dfield   = (double*) malloc(3*xyzcube*sizeof(double));
	
  /** set up lattice sites */
  FORALL_SITES(ix, iy, iz) {
    linearindex = maggs_get_linear_index(ix, iy, iz, lparams.dim);
		
    lattice[linearindex].r[0] = ix;
    lattice[linearindex].r[1] = iy;
    lattice[linearindex].r[2] = iz;

    FOR3D(i) {
      Bfield[3*linearindex+i]  = 0.;
      Dfield[3*linearindex+i]  = 0.;
    }

    lattice[linearindex].charge = 0.;
		
    /* set relative permittivity to 1 for all sites (default) */
    FOR3D(i) lattice[linearindex].permittivity[i] = 1.;
  }
	
  maggs_setup_neighbors();
}


/** sets up surface patches for all domains.
    @param surface_patch the local surface patch
*/
void maggs_calc_surface_patches(t_surf_patch* surface_patch)
{
  /* x=lparams.size[0] plane */
  surface_patch[0].offset   = lparams.dim[2]*lparams.dim[1]*lparams.size[0];    /*(size[0],0,0) point */
  surface_patch[0].doffset  = 0;                                             /*(0,0,0) point */
  surface_patch[0].stride   = lparams.dim[2]*lparams.dim[1];
  surface_patch[0].skip     = 0;
  surface_patch[0].nblocks  = 1;
  surface_patch[0].coord[0] = 2;
  surface_patch[0].coord[1] = 1;
  surface_patch[0].volume   = lparams.dim[2]*lparams.dim[1];
	
  /* x=1 plane */
  surface_patch[1].offset   = lparams.dim[2]*lparams.dim[1];                    /*(1,0,0) point */
  surface_patch[1].doffset  = lparams.dim[2]*lparams.dim[1]*lparams.inner_up_right[0];    /*(halo[0],0,0) point */
  surface_patch[1].stride   = lparams.dim[2]*lparams.dim[1];
  surface_patch[1].skip     = 0;
  surface_patch[1].nblocks  = 1;
  surface_patch[1].coord[0] = 2;
  surface_patch[1].coord[1] = 1;
  surface_patch[1].volume   = lparams.dim[2]*lparams.dim[1]; 
	
  /* y=lparams.size[1] plane */
  surface_patch[2].offset   = lparams.dim[2]*lparams.size[1];               /*(0,size[1],0) point */
  surface_patch[2].doffset  = 0;                                          /*(0,0,0) point */
  surface_patch[2].stride   = lparams.dim[2];
  surface_patch[2].skip     = lparams.dim[2]*lparams.dim[1];
  surface_patch[2].nblocks  = lparams.dim[0];  
  surface_patch[2].coord[0] = 2;
  surface_patch[2].coord[1] = 0;
  surface_patch[2].volume   = lparams.dim[2]*lparams.dim[0];
	
  /* y=1 plane */
  surface_patch[3].offset   = lparams.dim[2];                             /*(0,1,0) point */
  surface_patch[3].doffset  = lparams.dim[2]*lparams.inner_up_right[1];             /*(0,inner_up_right[1],0) point */
  surface_patch[3].stride   = lparams.dim[2];
  surface_patch[3].skip     = lparams.dim[2]*lparams.dim[1];
  surface_patch[3].nblocks  = lparams.dim[0];
  surface_patch[3].coord[0] = 2;
  surface_patch[3].coord[1] = 0;
  surface_patch[3].volume   = lparams.dim[2]*lparams.dim[0];
	
  /* z=lparams.size[2] plane */
  surface_patch[4].offset   = lparams.size[2];    /*(0,0,size[2]) point */
  surface_patch[4].doffset  = 0;                 /*(0,0,0) point */
  surface_patch[4].stride   = 1;
  surface_patch[4].skip     = lparams.dim[2];
  surface_patch[4].nblocks  = lparams.dim[0]*lparams.dim[1];
  surface_patch[4].coord[0] = 1;
  surface_patch[4].coord[1] = 0;
  surface_patch[4].volume   = lparams.dim[0]*lparams.dim[1];
	
  /* z=1 plane for z it must be higher*/
  surface_patch[5].offset   = 1;                   /*(0,0,1) point */
  surface_patch[5].doffset  = lparams.inner_up_right[2];     /*(0,0,inner_up_right[2]) point */
  surface_patch[5].stride   = 1;
  surface_patch[5].skip     = lparams.dim[2];
  surface_patch[5].nblocks  = lparams.dim[0]*lparams.dim[1];
  surface_patch[5].coord[0] = 1;
  surface_patch[5].coord[1] = 0;
  surface_patch[5].volume   = lparams.dim[0]*lparams.dim[1];
}


/** sets up MPI communications for domain surfaces */
void maggs_prepare_surface_planes(int dim, MPI_Datatype *xy, MPI_Datatype *xz, MPI_Datatype *yz, 
				  t_surf_patch *surface_patch)
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

/** get lattice size in one dimension
 @return mesh in 1D
 */
int maggs_get_mesh_1D()
{
    return maggs.mesh;
}

/** set permittivity for single lattice links
 @param node_x              index of the node in x direction
 @param node_y              index of the node in y direction
 @param node_z              index of the node in z direction
 @param direction           direction in which the link points from the node. 0 is for x, 1 is for y, 2 is for z
 @param relative_epsilon    permittivity to set, relative to the background permittivity set by the bjerrum length
 */
double maggs_set_permittivity(int node_x, int node_y, int node_z, int direction, double relative_epsilon)
{
    int node_index = maggs_get_linear_index(node_x, node_y, node_z, lparams.dim);
    /* save and return old epsilon for information purposes */
    double epsilon_before = lattice[node_index].permittivity[direction];
    /* set relative epsilon */
    lattice[node_index].permittivity[direction] = relative_epsilon;

    return epsilon_before;
}






/*****************************************/
/****** Surface patch communication ******/
/*****************************************/

/** MPI communication of surface region.
    works for D- and B-fields.
    @param field   Field to communicate. Can be B- or D-field.
    @param dim     Dimension in which to communicate
    @param e_equil Flag if field is already equilibated
*/
void maggs_exchange_surface_patch(double *field, int dim, int e_equil)
{
  static int init = 1;
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
	
		
  /** direction loop */
  for(s_dir=0; s_dir < 6; s_dir++) { 
    offset = dim * surface_patch[s_dir].offset;
    doffset= dim * surface_patch[s_dir].doffset;
		
    if(s_dir%2==0) r_dir = s_dir+1;
    else           r_dir = s_dir-1;
    /** pack send halo-plane data */
    if(node_neighbors[s_dir] != this_node) {
      /** communication */
      switch(s_dir) {
      case 0 :
      case 1 :
	if(e_equil || dim == 1) {
/* MPI_Sendrecv( sendbuf, sendcount, sendtype, dest, sendtag, recvbuf, recvcount, recvtype, source, recvtag, comm, status )
 
 if (rank%2) {
 MPI_Send(&data, 1, MPI_INT, next, 42, comm);
 MPI_Recv(&data, 1, MPI_INT, prev, 42, comm, 0);
 } else {
 MPI_Recv(&data, 1, MPI_INT, prev, 42, comm, 0);
 MPI_Send(&data, 1, MPI_INT, next, 42, comm);
 }
 
 MPI_Sendrecv(&send_data, 1, MPI_FLOAT, next, 42, &recv_data, 1, MPI_FLOAT, prev, 42, comm, 0);
 MPI_Sendrecv_replace(&data, 1, MPI_FLOAT, next, 42, prev, 42, comm, 0);
 */
	  MPI_Irecv (&field[doffset],1,yzPlane,node_neighbors[s_dir],REQ_MAGGS_SPREAD,comm_cart,&request[0]);
	  MPI_Isend(&field[offset],1,yzPlane,node_neighbors[r_dir],REQ_MAGGS_SPREAD,comm_cart,&request[1]);
	}
	else {
	  MPI_Irecv (&field[doffset+1],1,yzPlane2D,node_neighbors[s_dir],REQ_MAGGS_SPREAD,comm_cart,&request[0]);
	  MPI_Isend(&field[offset+1],1,yzPlane2D,node_neighbors[r_dir],REQ_MAGGS_SPREAD,comm_cart,&request[1]);
	}	  
					
	MPI_Waitall(2,request,status);
	break;
      case 2 :
      case 3 :
	if(e_equil || dim == 1) {
	  MPI_Irecv (&field[doffset],1,xzPlane,node_neighbors[s_dir],REQ_MAGGS_SPREAD,comm_cart,&request[0]);
	  MPI_Isend(&field[offset],1,xzPlane,node_neighbors[r_dir],REQ_MAGGS_SPREAD,comm_cart,&request[1]);
	}
	else {
	  MPI_Irecv (&field[doffset],1,xzPlane2D,node_neighbors[s_dir],REQ_MAGGS_SPREAD,comm_cart,&request[0]);
	  MPI_Isend(&field[offset],1,xzPlane2D,node_neighbors[r_dir],REQ_MAGGS_SPREAD,comm_cart,&request[1]);
	}	  
	MPI_Waitall(2,request,status);
	break;
      case 4 :
      case 5 : 
	if(e_equil || dim == 1) {
	  MPI_Irecv (&field[doffset],1,xyPlane,node_neighbors[s_dir],REQ_MAGGS_SPREAD,comm_cart,&request[0]);
	  MPI_Isend(&field[offset],1,xyPlane,node_neighbors[r_dir],REQ_MAGGS_SPREAD,comm_cart,&request[1]);
	}
	else {
	  MPI_Irecv (&field[doffset],1,xyPlane2D,node_neighbors[s_dir],REQ_MAGGS_SPREAD,comm_cart,&request[0]);
	  MPI_Isend(&field[offset],1,xyPlane2D,node_neighbors[r_dir],REQ_MAGGS_SPREAD,comm_cart,&request[1]);
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







/*********************************************/
/****** interpolate charges on lattice: ******/
/*********************************************/

/** Interpolation function in one dimension.
    Currently only linear interpolation conserves charge.
    @return interpolation value
    @param x relative position of particle
*/
double maggs_interpol1D(double x)
{	
  return x;
  /* Cosine interpolation would be: */
  /* return sqr(sin(M_PI_2*x)); */
}


/** Does the actual interpolating calculation.
    @param first  3dim-array lattice position,
    @param rel    3dim-array relative position in cube,
    @param q      charge to interpolate
*/
void maggs_interpolate_charge(int *first, double *rel, double q)
{
  int i, k, l, m, index, temp_ind;
  int help_index[3];
  double temp;
  double help[SPACE_DIM];
	
  FOR3D(i) help[i] = 1. - rel[i];     /** relative pos. w.r.t. first */
  //  printf("first: %d %d %d\n", first[0], first[1], first[2]);
  /** calculate charges at each vertex */
  index = maggs_get_linear_index(first[0],first[1],first[2],lparams.dim);
	
  FOR3D(i) {
    temp_ind = neighbor[index][i];
    if(temp_ind == NOWHERE) help_index[i] = lparams.volume; /* force huge index */
    else { /* incr. for x-neighbor */
      help_index[i] = maggs_get_offset(lattice[neighbor[index][i]].r[i], first[i], i, lparams.dim);
    }
		
  }
	
  for(k=0;k<2;k++){   /* jumps from x- to x+ */
    for(l=0;l<2;l++){  /* jumps from y- to y+ */
      for(m=0;m<2;m++){ /* jumps from z- to z+ */      
	if(index < lparams.volume) {
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

/** add charges from ghost cells to lattice sites. */
void maggs_accumulate_charge_from_ghosts()
{
  Cell *cell;
  Particle* p;
  int i, c, d;
  int np;
  int flag_inner=0;
  int first[SPACE_DIM];
  double q;
  double pos[SPACE_DIM], rel[SPACE_DIM];
	
  /** loop over ghost cells */
  for (c = 0; c < ghost_cells.n; c++) {
    cell = ghost_cells.cell[c];
    p  = cell->part;
    np = cell->n;
    for(i = 0; i < np; i++) {
      if( (q=p[i].p.q) != 0.0 ) {
	flag_inner=1;
	FOR3D(d) {
	  if(p[i].r.p[d]<lparams.left_down_position[d]||p[i].r.p[d]>=lparams.upper_right_position[d])
	    {flag_inner=0; break;}
	}
      }
      if(flag_inner) {
	FOR3D(d) {
	  pos[d]        = (p[i].r.p[d] - lparams.left_down_position[d])* maggs.inva;
	  first[d]      = (int) pos[d];
	  rel[d]        = pos[d] - first[d];
	}
	//      	fprintf(stderr,"pos: %f %f %f\n", p[i].r.p[0], p[i].r.p[1], p[i].r.p[2]);
	maggs_interpolate_charge(first, rel, q);
      }
    }      
  }  
	
}


/** finds current lattice site of each particle.
    calculates charge interpolation on cube. */
void maggs_distribute_particle_charges()
{	
  Cell *cell;
  Particle* p;
  int i, c, d;
  int np;
  int first[SPACE_DIM];
  double q;
  double pos[SPACE_DIM], rel[SPACE_DIM];
		
  for(i=0;i<lparams.volume;i++) lattice[i].charge = 0.;
  /** === charge assignment === */ 
  /** loop over inner cells */
  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    p  = cell->part;
    np = cell->n;
    for(i = 0; i < np; i++) {
      if( (q=p[i].p.q) != 0.0 ) {
	FOR3D(d) {
	  pos[d]        = (p[i].r.p[d] - lparams.left_down_position[d])* maggs.inva;
	  first[d]      = (int) pos[d];
	  rel[d]        = pos[d] - first[d];
	}
	maggs_interpolate_charge(first, rel, q);
      }
    }      
  }
  maggs_accumulate_charge_from_ghosts();	
}
	
/** Does the actual calculation of the gradient. Parameters:
    @param rel  3dim-array of relative position in cube,
    @param q    charge,
    @param grad huge gradient array to write into
*/
void maggs_calc_charge_gradients(double *rel, double q, double *grad)
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



/** finds correct lattice cube for particles.
    calculates charge gradients on this cube.
    @param grad gradient array to write into
*/
void maggs_update_charge_gradients(double *grad)
{
  Cell *cell;
  Particle* p;
  int i, c, d, ip;
  int np;
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
	  pos[d]        = (p[i].r.p[d] - lparams.left_down_position[d])* maggs.inva;
	  first[d]      = (int) pos[d];
	  rel[d]        = pos[d] - first[d];
	}
	maggs_calc_charge_gradients(rel, q, &grad[ip]);
	ip += 12;
      }
    }      
  }  
  
} 







/***************************************/
/****** initialization procedure: ******/
/***************************************/

/** Calculate the self energy coefficients for the system, if
    corrected with Lattice Green's function. */
void maggs_calc_self_energy_coeffs()
{
  double factor, prefac;
  int px = 0;
  int py = 0;
  int pz = 0;
  int i, j, k, l, m, index;
  double sx = 0.;
  double sy = 0.;
  double sz = 0.;
  double sxy = 0.;
  double sxyz = 0.;
  double nomx, nomy, nomz;
  double inva = maggs.inva;
  double invasq = inva*inva;
  double n[8][SPACE_DIM];// = {{0.,0.,0.}, {1.,0.,0.}, {0.,1.,0.}, {1.,1.,0.}, 
  // {0.,0.,1.}, {1.,0.,1.}, {0.,1.,1.}, {1.,1.,1.}};

  m=0;
  l=0;
  index = 0;
  for(i=0;i<2;i++)
    for(j=0;j<2;j++)
      for(k=0;k<2;k++) {
	n[index][2] = m; 
	n[index][1] = l; 
	n[index][0] = k; 
	index++;
  }

  factor = M_PI / maggs.mesh;
  prefac = 1. / maggs.mesh;
  prefac = 0.5 * prefac * prefac * prefac * SQR(maggs.prefactor);
  
  for(i=0;i<8;i++) 
    {
      for(j=0;j<8;j++) 
	{
	  maggs.alpha[i][j] = 0.;
	  for(px = 0; px < maggs.mesh; ++px)
	    {
	      sx = sin( factor * px );
	      sx = sx * sx;
	      nomx = 2.*factor*px*(n[i][0] - n[j][0]);
	      for(py = 0; py < maggs.mesh; ++py)
		{
		  sy = sin( factor * py );
		  sy = sy * sy; 
		  nomy = 2.*factor*py*(n[i][1] - n[j][1]);
		  sxy = sx + sy;
		  for(pz = 0; pz < maggs.mesh; ++pz)
		    {
		      sz = sin( factor * pz );
		      sz = sz * sz;
		      nomz = 2.*factor*pz*(n[i][2] - n[j][2]);
		      sxyz = sxy + sz;
		      sxyz *= 4.;
		      if(sxyz > 0)
			{
			  maggs.alpha[i][j] += cos(nomx + nomy + nomz) / sxyz;
 			}
		    }
		}
	    }
	  /* invasq is needed for the calculation of forces */
	  maggs.alpha[i][j] = invasq * prefac * maggs.alpha[i][j];
	}
    }
		      
}

/** For energy minimization:
    @return maximum of curl(E) from all sites.
    May also be used to verify constraint surface condition. */
double maggs_check_curl_E()
{
  int i, ix, iy, iz;
  double curl, maxcurl, gmaxcurl;
  int* anchor_neighb;
	
  maxcurl = 0.;
	
  FORALL_INNER_SITES(ix, iy, iz) {
    i = maggs_get_linear_index(ix, iy, iz, lparams.dim); 
    anchor_neighb = neighbor[i];
    curl = Dfield[3*i] + Dfield[3*anchor_neighb[0]+1] 
      - Dfield[3*anchor_neighb[1]] - Dfield[3*i+1];
    curl *= maggs.inva;
    if(fabs(curl)>maxcurl) maxcurl = fabs(curl);
    curl = Dfield[3*i+2] + Dfield[3*anchor_neighb[2]] 
      - Dfield[3*anchor_neighb[0]+2] - Dfield[3*i];
    curl *= maggs.inva;
    if(fabs(curl)>maxcurl) maxcurl = fabs(curl);
    curl = Dfield[3*i+1] + Dfield[3*anchor_neighb[1]+2] 
      - Dfield[3*anchor_neighb[2]+1] - Dfield[3*i+2];
    curl *= maggs.inva;
    if(fabs(curl)>maxcurl) maxcurl = fabs(curl);
  }
  MPI_Allreduce(&maxcurl,&gmaxcurl,1,MPI_DOUBLE,MPI_MAX,comm_cart);  
  return gmaxcurl;
}

/** If curl(E) is not zero, update plaquette with corrections.
    @param i index of the current lattice site
    @param n coordinate, is the normal direction to the plaquette
*/
void maggs_perform_rot_move_inplane(int i, int n)
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
	
  delta = Dfield[3*i+mue] + Dfield[3*anchor_neighb[mue]+nue] 
    - Dfield[3*anchor_neighb[nue]+mue] - Dfield[3*i+nue];
  if(fabs(delta)>=ROUND_ERR) {
    delta = -delta/4.; 
    maggs_update_plaquette(mue, nue, anchor_neighb, i, delta);
  }
}


/** For energy minimization:
    Relax B-fields in all directions.*/
void maggs_minimize_transverse_field()
{
  int k, l, m;
  int i, d;
  int ind_i, ind_j;
  int size[2]={0,0};
  int index = -1; /* force allocation error */
	
  FOR3D(d) {
    switch(d) {
    case 0 :
      size[0] = lparams.size[2];
      size[1] = lparams.size[1];
      break;
    case 1 :
      size[0] = lparams.size[2];
      size[1] = lparams.size[0];
      break;
    case 2 :
      size[0] = lparams.size[1];
      size[1] = lparams.size[0];
      break;
    }
    for(i=0;i<2;i++) {
      /* at first even sites (i==0) then odd */
      for(k=1;k<=lparams.size[d];k++) {
	/* update every plane in direction d */
	ind_i=0;
	for(l=0; l<=size[1]; l++){
	  ind_j=0;
	  for(m=0;m<=size[0]; m++) {
	    switch(d) {
	    case 0 :
	      index=maggs_get_linear_index(k,l,m,lparams.dim);
	      break;
	    case 1 :
	      index=maggs_get_linear_index(l,k,m,lparams.dim); 
	      break;
	    case 2 :
	      index=maggs_get_linear_index(l,m,k,lparams.dim); 
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
      maggs_exchange_surface_patch(Dfield, 3, 0);
    }
  }
}


/** calculates initial electric field configuration.
    currently uses simple and slow method of plaquettes and links.
    energy minimization takes up lots of time. */
void maggs_calc_init_e_field()
{
  double localqy, localqz;
  double qplane, qline;
  int    i, k, ix, iy, iz;
  int index = 0;
  double invasq, tmp_field=0.0;
  double sqrE, gsqrE, gavgEx, gavgEy, gavgEz;
  double qz, qy, qx, avgEx, avgEy, avgEz;
  double Eall[SPACE_DIM], gEall[SPACE_DIM];
  double maxcurl;
  MPI_Status status;
  MPI_Comm zplane, yline;
  int color, rank, dim;
	
  MAGGS_TRACE(fprintf(stderr,"%d:Initialize field\n",this_node));
	
  invasq = maggs.inva*maggs.inva;
	
  /* sort particles for the calculation of initial charge distribution */
  MAGGS_TRACE(fprintf(stderr,"%d:Sorting particles...\n",this_node));
	
  cells_resort_particles(CELL_GLOBAL_EXCHANGE);
  /*  fprintf(stderr, "done\n"); */
  maggs_distribute_particle_charges();
	
  dim = node_grid[1]*node_grid[0];
  color = node_pos[2];
  rank  = this_node%dim;
  MPI_Comm_split(comm_cart, color, rank, &zplane);
  color = node_pos[1];
  rank  = rank%node_grid[0];
  MPI_Comm_split(zplane, color, rank, &yline);
	
  /* calculate initial solution of Poisson equation */
	
  /* CAUTION: the indexing of the neighbor nodes in Espresso
   *  starts from x left neighbor node
   */
	
  /* get process coordinates */
  if(node_pos[2]!= 0) {
    MPI_Recv(&tmp_field, 1, MPI_DOUBLE, node_neighbors[4], REQ_MAGGS_EQUIL, comm_cart, &status);
    for(iy=lparams.inner_left_down[1];iy<lparams.inner_up_right[1];iy++) {
      for(ix=lparams.inner_left_down[0];ix<lparams.inner_up_right[0];ix++) {  
	index = maggs_get_linear_index(ix, iy, lparams.inner_left_down[2], lparams.dim);
	Dfield[3*neighbor[index][ZMINUS]+ZPLUS] = tmp_field;
      }
    }
  }

  localqz = 0.;
  for(iz=lparams.inner_left_down[2];iz<lparams.inner_up_right[2];iz++) { /* loop over z-planes */
    localqz = 0.;
    for(iy=lparams.inner_left_down[1];iy<lparams.inner_up_right[1];iy++) {
      for(ix=lparams.inner_left_down[0];ix<lparams.inner_up_right[0];ix++) {  
	index = maggs_get_linear_index(ix, iy, iz, lparams.dim);
	/* Sum over the charge of all sides in z-plane */
	localqz += lattice[index].charge / lattice[index].permittivity[2];
      }
    } 

    MPI_Allreduce(&localqz, &qz, 1, MPI_DOUBLE, MPI_SUM, zplane);
    qz = qz/(maggs.mesh*maggs.mesh);
    qplane = qz*maggs.prefactor*invasq;
    /*    if(fabs(qplane) >= 0.01*ROUND_ERROR_PREC) { */
    for(iy=lparams.inner_left_down[1];iy<lparams.inner_up_right[1];iy++) {
      for(ix=lparams.inner_left_down[0];ix<lparams.inner_up_right[0];ix++) {  
	index = maggs_get_linear_index(ix, iy, iz, lparams.dim);
	Dfield[3*index+ZPLUS]  = Dfield[3*neighbor[index][ZMINUS]+ZPLUS] + qplane;
	/*	    + qz*maggs.prefactor*invasq; */
      }
    }
    /*    } */
    if(iz>=lparams.inner_up_right[2]-1) {
      if (node_pos[2]<node_grid[2]-1) {
	if(node_grid[2]>1) {
	  MPI_Send(&Dfield[3*index+ZPLUS], 1, MPI_DOUBLE, node_neighbors[5], REQ_MAGGS_EQUIL, comm_cart); 
	}
      }
      else 
	if (fabs(Dfield[3*index+ZPLUS]) > 100.*ROUND_ERROR_PREC) {
	  fprintf(stderr, "%d: Error in the calculation of Ez(%d,%d,%d)=%f!!\n", 
		  this_node,lattice[index].r[0], lattice[index].r[1], lattice[index].r[2],
		  Dfield[3*index+ZPLUS]);
	  fflush(stderr);
	}
    }

    if(node_pos[1]!= 0) {
      MPI_Recv(&tmp_field, 1, MPI_DOUBLE, node_neighbors[2], REQ_MAGGS_EQUIL, comm_cart, &status);
      for(ix=lparams.inner_left_down[0];ix<lparams.inner_up_right[0];ix++) {  
	index = maggs_get_linear_index(ix, lparams.inner_left_down[1], iz, lparams.dim);
	Dfield[3*neighbor[index][YMINUS]+YPLUS] = tmp_field;
      }
    }
		
    for(iy=lparams.inner_left_down[1];iy<lparams.inner_up_right[1];iy++) {
      localqy = 0.;
      for(ix=lparams.inner_left_down[0];ix<lparams.inner_up_right[0];ix++) {  
	index = maggs_get_linear_index(ix, iy, iz, lparams.dim);
	localqy += lattice[index].charge / lattice[index].permittivity[1];
      }

      MPI_Allreduce(&localqy, &qy, 1, MPI_DOUBLE, MPI_SUM, yline);
			
      qy = qy/maggs.mesh;
      qline = (qy-qz)*maggs.prefactor*invasq;
      /*      if(fabs(qy-qz)>=ROUND_ERROR_PREC) { */
      for(ix=lparams.inner_left_down[0];ix<lparams.inner_up_right[0];ix++) {  
	index = maggs_get_linear_index(ix, iy, iz, lparams.dim);
	Dfield[3*index+YPLUS]  = Dfield[3*neighbor[index][YMINUS]+YPLUS] + qline;
	/*	    (qy-qz)*maggs.prefactor*invasq; */
      }
      /*      } */
			
      if(iy>=lparams.inner_up_right[1]-1) {
	if(node_pos[1] < node_grid[1]-1) {
	  if (node_grid[1]>1)
	    MPI_Send(&Dfield[3*index+YPLUS], 1, MPI_DOUBLE, node_neighbors[3], REQ_MAGGS_EQUIL, comm_cart); 
	}
	else
	  if (fabs(Dfield[3*index+YPLUS]) > 100.*ROUND_ERROR_PREC)
	    fprintf(stderr, "%d: Error in the calculation of Ey(%d,%d,%d)=%f!!\n",
		    this_node, lattice[index].r[0], lattice[index].r[1], lattice[index].r[2],
		    Dfield[3*index+YPLUS]);	  
      }
			
      if(node_pos[0]!= 0) {
	MPI_Recv(&tmp_field, 1, MPI_DOUBLE, node_neighbors[0], REQ_MAGGS_EQUIL, comm_cart, &status);
	index = maggs_get_linear_index(lparams.inner_left_down[0], iy, iz, lparams.dim);
	Dfield[3*neighbor[index][XMINUS]+XPLUS] = tmp_field;
      }
			
      for(ix=lparams.inner_left_down[0];ix<lparams.inner_up_right[0];ix++) {  
	index = maggs_get_linear_index(ix, iy, iz, lparams.dim);
	qx = lattice[index].charge / lattice[index].permittivity[0]; 
	Dfield[3*index+XPLUS] = Dfield[3*neighbor[index][XMINUS]+XPLUS] + 
	  (qx-qy)*maggs.prefactor*invasq;
      }
			
      if(ix>=lparams.inner_up_right[0]-1) {
	if(node_pos[0] < node_grid[0]-1) {
	  if(node_grid[0]>1)
	    MPI_Send(&Dfield[3*index+XPLUS], 1, MPI_DOUBLE, node_neighbors[1], REQ_MAGGS_EQUIL, comm_cart); 
	}
	else
	  if (fabs(Dfield[3*index+XPLUS]) > 100.*ROUND_ERROR_PREC)
	    fprintf(stderr, "%d: Error in the calculation of Ex(%d,%d,%d)=%f!!\n",
		    this_node, lattice[index].r[0], lattice[index].r[1], lattice[index].r[2],
		    Dfield[3*index+XPLUS]);	  
      }
    }    /*** loop over iy */
  }
	
  /* exchange halo-surfaces */
  maggs_exchange_surface_patch(Dfield, 3, 1); 
	
  avgEz = 0.;
  for(iz=lparams.inner_left_down[2];iz<lparams.inner_up_right[2];iz++) {
    index = maggs_get_linear_index(lparams.inner_left_down[0], lparams.inner_left_down[1], iz, lparams.dim);
    avgEz += Dfield[3*index+ZPLUS];
  }
	
	
  /*  MPI_Barrier(comm_cart); */
  MPI_Allreduce(&avgEz,&gavgEz,1,MPI_DOUBLE,MPI_SUM,comm_cart);
  gavgEz = gavgEz/(maggs.mesh*node_grid[0]*node_grid[1]);
	
  FORALL_INNER_SITES(ix, iy,iz) {
    index = maggs_get_linear_index(ix, iy, iz, lparams.dim);
    Dfield[3*index+ZPLUS] -= gavgEz;
  }
	
  for(iz = lparams.inner_left_down[2];iz<lparams.inner_up_right[2];iz++) {
    avgEy = 0.;  
    for(iy = lparams.inner_left_down[1];iy<lparams.inner_up_right[1];iy++) {
      index = maggs_get_linear_index(lparams.inner_left_down[0], iy, iz, lparams.dim);
      avgEy += Dfield[3*index+YPLUS];
    }    
		
    MPI_Allreduce(&avgEy, &gavgEy, 1, MPI_DOUBLE, MPI_SUM, zplane);
    gavgEy = gavgEy/(maggs.mesh*node_grid[0]);
		
    for(iy=lparams.inner_left_down[1];iy<lparams.inner_up_right[1];iy++) {
      for(ix=lparams.inner_left_down[0];ix<lparams.inner_up_right[0];ix++)  
	Dfield[3*maggs_get_linear_index(ix, iy, iz, lparams.dim)+YPLUS] -= gavgEy;
    }
  }
	
  for(iz=lparams.inner_left_down[2];iz<lparams.inner_up_right[2];iz++) {
    for(iy=lparams.inner_left_down[1];iy<lparams.inner_up_right[1];iy++) {
      avgEx = 0.;
      for(ix=lparams.inner_left_down[0];ix<lparams.inner_up_right[0];ix++) {
	avgEx += Dfield[3*maggs_get_linear_index(ix, iy, iz, lparams.dim)+XPLUS];
      }
			
      MPI_Allreduce(&avgEx, &gavgEx, 1, MPI_DOUBLE, MPI_SUM, yline);
      gavgEx = gavgEx/maggs.mesh;
			
      for(ix=lparams.inner_left_down[0];ix<lparams.inner_up_right[0];ix++)  
	Dfield[3*maggs_get_linear_index(ix, iy, iz, lparams.dim)+XPLUS] -= gavgEx;
    }
  }
	
	
  /* exchange halo-surfaces */
  maggs_exchange_surface_patch(Dfield, 3, 1);
	

	
  MPI_Comm_free(&zplane);
  MPI_Comm_free(&yline);
	
	
	
	
  /* iterative procedure of energy minimization */

  sqrE = 0.;
  FORALL_INNER_SITES(ix, iy, iz) {
    i = maggs_get_linear_index(ix, iy, iz, lparams.dim);
    FOR3D(k) sqrE += SQR(Dfield[3*i+k]);
  }

  MPI_Allreduce(&sqrE,&gsqrE,1,MPI_DOUBLE,MPI_SUM,comm_cart); 
  gsqrE = gsqrE/(SPACE_DIM*maggs.mesh*maggs.mesh*maggs.mesh);  

#ifdef MAGGS_DEBUG
  int iteration; 
  if(!this_node) iteration = 0;
#endif
  do {
#ifdef MAGGS_DEBUG
    double goldE = gsqrE;
#endif
    sqrE = 0.;
    maggs_minimize_transverse_field();
		
    FORALL_INNER_SITES(ix, iy, iz) {
      i = maggs_get_linear_index(ix, iy, iz, lparams.dim);
      FOR3D(k) sqrE += SQR(Dfield[3*i+k]);
    }
    MPI_Allreduce(&sqrE,&gsqrE,1,MPI_DOUBLE,MPI_SUM,comm_cart); 
    gsqrE = gsqrE/(SPACE_DIM*maggs.mesh*maggs.mesh*maggs.mesh);  
    maxcurl = maggs_check_curl_E();
		
#ifdef MAGGS_DEBUG
    if(!this_node) {
      iteration++;
      if(iteration%1==0) {
	fprintf(stderr, "# iteration for field equilibration %d, diff=%9.4e, curlE=%9.4e\n", 
		iteration, fabs(gsqrE-goldE),maxcurl);
	fflush(stderr);    
      }
    }
#endif
		
  } while(fabs(maxcurl)>1000000.*ROUND_ERROR_PREC);
	

	
	
	
	
	
	
	
  /* exchange halo-surfaces */
	
  FOR3D(k) Eall[k] = 0.;
  FORALL_INNER_SITES(ix, iy, iz) {
    i = maggs_get_linear_index(ix, iy, iz, lparams.dim);
    FOR3D(k) {
      Eall[k] += Dfield[3*i+k];
    }
  }
	
  MPI_Allreduce(Eall,gEall,3,MPI_DOUBLE,MPI_SUM,comm_cart);
	
  FOR3D(k) gEall[k] /= (n_nodes*lparams.size[0]*lparams.size[1]*lparams.size[2]);
	
  FORALL_INNER_SITES(ix, iy, iz) {
    i = maggs_get_linear_index(ix, iy, iz, lparams.dim);
    FOR3D(k) Dfield[3*i+k] -= gEall[k];
  }
  /* exchange whole glue-patch region */
  maggs_exchange_surface_patch(Dfield, 3, 0);
  if(!this_node)
    MAGGS_TRACE(fprintf(stderr, "Ex = %16.12e, Ey = %15.12e, Ez = %15.12e\n", gEall[0], gEall[1], gEall[2]));
}








/*********************************************/
/****** calculate currents and E-fields ******/
/*********************************************/

/** calculate the charge flux in direction "dir".
    @param q     charge of the moving particle
    @param help  the relative position in the cube to the opposite of left down front lattice site.
    @param flux  flux variable to write into
    @param dir   direction in which to calculate the flux
*/
void maggs_calc_charge_fluxes_1D(double q, double *help, double *flux, int dir)
{
  /* at the moment works only for linear interpolation */
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

/** Extend particle trajectories on the one time step and check if the
    trajectory intersects the cell boundary in direction dirx#
    @return zero if no particle crosses intersection
    @param delta    amount by which the particle moves
    @param r_new    current position of particle
    @param dir      direction in which to check for intersection
    @param first    position of particle within lattice cube
    @param t_step   time step
    @param identity particle ID
*/
int maggs_check_intersect_1D(double delta, double r_new, int dir, int first, double *t_step, int identity)
{	
  int candidateplane = -1; /* force alloc error */
  int f_crossing; 
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

/** updates field on link coupling it with current.
    Force is multiplied by the time_step
    @param index   index of lattice site
    @param flux    charge flux in lattice site
    @param v       speed of particle
    @param dir     first direction of flux calculation
*/
void maggs_calc_e_field_on_link_1D(int index, double *flux, double v, int dir)
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
  if(temp_ind == NOWHERE) help_index[0] = lparams.volume;
  else
    help_index[0] = maggs_get_offset(lattice[temp_ind].r[dir1], anchor_site->r[dir1], dir1, lparams.dim);
  temp_ind = anchor_neighb[dir2];
  if(temp_ind == NOWHERE) help_index[1] = lparams.volume;
  else
    help_index[1] = maggs_get_offset(lattice[temp_ind].r[dir2], anchor_site->r[dir2], dir2, lparams.dim);
	
	
  ind_flux = 0;
  for(l=0;l<2;l++){  /* jumps from dir2- to dir2+ */
    for(m=0;m<2;m++){ /* jumps from dir1- to dir1+ */  
			
      if(index < lparams.volume){
	Dfield[3*index+dir] -= flux[ind_flux] * v;
      }
			
      ind_flux++; 
			
      index+=help_index[0];
      help_index[0]=-help_index[0];	
    }
    index+=help_index[1];
    help_index[1]=-help_index[1];     
  }
}


/** loop over all cells and call charge current functions
    @param p           Particle pointer
    @param ghost_cell  flag if cell is ghost cell
*/
void maggs_add_current_on_segment(Particle *p, int ghost_cell)
{
  int d;
  int icoord, dir;
  int lat_index = -1; /* force alloc error */
  int first[SPACE_DIM];
  int f_crossing, flag_update_flux;
  t_dvector r_temp; 
  double delta;
  double flux[4], v;
  double v_inva[SPACE_DIM], v_invasq[SPACE_DIM];
  double pos[SPACE_DIM], help[3];
  double t_step;
	
  FOR3D(d) {
    pos[d]   = (p->r.p[d] - lparams.left_down_position[d])* maggs.inva;
    first[d] = (int) floor(pos[d]);
    r_temp[d]   = pos[d] - first[d]; /* it is the updated coord (we have to go back) */
    help[d]     = 1. - r_temp[d];
    v_inva[d]   = maggs.inva * p->m.v[d];
    v_invasq[d] = maggs.inva * v_inva[d];
  }
	
  flag_update_flux = 1;
  if(ghost_cell) {
    FOR3D(d) if(first[d]<lparams.halo_left_down[d] || first[d]>= lparams.halo_upper_right[d])
      {flag_update_flux = 0;break;}
  }
	
  if(flag_update_flux) {
    lat_index = maggs_get_linear_index(first[0], first[1], first[2], lparams.dim);
  }
	
  /* loop coordinates in order x->y->z->y->x */
  for(dir=0; dir<5; dir++) {
    icoord = dir;
    if(dir>2) icoord = dir%2;
    if(icoord == 2) delta = v_inva[icoord];
    else            delta = 0.5 * v_inva[icoord];
		
    f_crossing = maggs_check_intersect_1D(delta, r_temp[icoord], icoord, first[icoord], &t_step, p->p.identity);
		
    /* calculate flux */
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
	  if(first[icoord]<lparams.halo_left_down[icoord] || first[icoord]>= lparams.halo_upper_right[icoord])
	    {flag_update_flux = 0;}
	}
	else {
	  flag_update_flux = 1;
	  FOR3D(d) if(first[d]<lparams.halo_left_down[d] || first[d]>= lparams.halo_upper_right[d])
	    {flag_update_flux = 0;break;}
	  if(flag_update_flux) maggs_calc_charge_fluxes_1D(p->p.q, help, flux, icoord); 
	}
      }
			
      if(flag_update_flux) {
	v = t_step * v_invasq[icoord];
	lat_index = maggs_get_linear_index(first[0], first[1], first[2], lparams.dim);
	maggs_calc_e_field_on_link_1D(lat_index, flux, v, icoord);
      }
    }
    r_temp[icoord] -= delta;
    help[icoord]    = 1. - r_temp[icoord];
  }
}


/** Calculate fluxes and couple them with fields symplectically.  It
    is assumed that the particle can not cross more than one cell
    boundary per direction */
void maggs_couple_current_to_Dfield()
{
  Cell *cell;
  Particle* p;
  int i, c, d, np;
  int flag_inner;
  double q;
  double r1, r2;
	
  /** loop over real particles */
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
	
  /** loop over ghost particles */
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
	  if(((r2 < lparams.left_down_position[d])&&(r1 < lparams.left_down_position[d]))
	     ||((r2 >= lparams.upper_right_position[d] && r1 >= lparams.upper_right_position[d])))
	    {flag_inner = 0; break;}
	}
	if(flag_inner) {
	  maggs_add_current_on_segment(&p[i], 1);
	}
      }/* if particle.q != ZERO */
    }
  }
}









/*******************************************/
/****** calculate B-fields and forces ******/
/*******************************************/

/** propagate the B-field via \f$\frac{\partial}{\partial t}{B} = \nabla\times D\f$ (and prefactor)
    CAREFUL: Usually this function is called twice, with dt/2 each time
    to ensure a time reversible integration scheme!
    @param dt time step for update. Should be half the MD time step
*/
void maggs_propagate_B_field(double dt)
{
  int x, y, z, i, offset, index;
  int xoffset, yoffset;
  double help = dt*maggs.invsqrt_f_mass;
  /* B(t+h/2) = B(t-h/2) + h*curlE(t) */ 
	
  offset = maggs_get_linear_index(1,1,1, lparams.dim);
  yoffset = lparams.dim[2];
  xoffset = 2*lparams.dim[2];
	
  for(x=0;x<lparams.size[0];x++) {
    for(y=0;y<lparams.size[1];y++) {
      for(z=0;z<lparams.size[2];z++) {
				
	i = offset+z;
	index = 3*i;
	Bfield[index+0] += - help*maggs_calc_dual_curl(1,2, Dfield, neighbor[i], index); 
	Bfield[index+1] += - help*maggs_calc_dual_curl(2,0, Dfield, neighbor[i], index); 
	Bfield[index+2] += - help*maggs_calc_dual_curl(0,1, Dfield, neighbor[i], index);  
      }
      offset += yoffset;
    }
    offset += xoffset;
  }
	
  maggs_exchange_surface_patch(Bfield, 3, 0);
}

/** calculate D-field from B-field according to
    \f$\frac{\partial}{\partial t}{D} = \nabla\times B\f$ (and prefactors)
    @param dt MD time step
*/
void maggs_add_transverse_field(double dt)
{
  int i, index;
  double invasq; 
  int x, y, z;
  int offset, xoffset, yoffset;
  double help;
	
  invasq = SQR(maggs.inva);
  help = dt * invasq * maggs.invsqrt_f_mass;
	
  /***calculate e-field***/ 
  offset = maggs_get_linear_index(1,1,1, lparams.dim);
  yoffset = lparams.dim[2];
  xoffset = 2*lparams.dim[2];
  for(x=0;x<lparams.size[0];x++) {
    for(y=0;y<lparams.size[1];y++) {
      for(z=0;z<lparams.size[2];z++) {
	/*  FORALL_INNER_SITES(x, y, z) { */
	/*    i = maggs_get_linear_index(x, y, z, lparams.dim); */
	i = offset+z;
	index = 3*i;
	Dfield[index  ] += help * maggs_calc_curl(2, 1, Bfield, neighbor[i], index);
	Dfield[index+1] += help * maggs_calc_curl(0, 2, Bfield, neighbor[i], index);
	Dfield[index+2] += help * maggs_calc_curl(1, 0, Bfield, neighbor[i], index);
      }
      offset += yoffset;
    }
    offset += xoffset;
  } 
	
  maggs_exchange_surface_patch(Dfield, 3, 0);
}


/** interpolation function, solely for new self force correction. */
void maggs_interpolate_part_charge_from_grad(double *rel, double *grad, double *rho)
{
  int i, k, l, m, index;
  int grad_ind;

  int help_index[3];
  double help[3];

  FOR3D(i) help[i] = 1. - rel[i];     /* relative pos. w.r.t. first */  

  help_index[0] = 4;
  help_index[1] = 2; 
  help_index[2] = 1;

  grad_ind = 0;
  index = 0;
  for(i=0;i<8;i++) rho[i] = 0.;

  for(k=0;k<2;k++){   /* jumps from x- to x+ */
    for(l=0;l<2;l++){  /* jumps from y- to y+ */
      for(m=0;m<2;m++){ /* jumps from z- to z+ */
	// without q!!!
	if(k==0) rho[index] += - help[0] * grad[grad_ind];
	else {
	  rho[index] += - rel[0] * grad[grad_ind%4];
	}

	grad_ind ++;
	index+=help_index[2];
	help_index[2]=-help_index[2];
      }
      index+=help_index[1];
      help_index[1]=-help_index[1];
    }
    index+=help_index[0];
    help_index[0]=-help_index[0];
  }
  //  for(i=0;i<8;i++) printf("rho: %f\n", rho[i]);
}

/** new self force correction with lattice Green's function.
    @param p Particle pointer
*/
void maggs_calc_interpolated_self_force(Particle *p)
{
// TODO : FIX WARNING ABOUT UNUSED VARIABLES
// FUCTION CURRENTLY NOT USED IN CODE.
/*

  int ix,iy,iz,k,index,globalindex;
  int xmax=2, ymax=2, zmax=2;
  double self_force[SPACE_DIM];
  double position[SPACE_DIM];
  double relative_position[SPACE_DIM];
  int left_down_position[SPACE_DIM];

  double local_rho[8];
  double local_permittivity[12];
  double local_D_field[12];
  double *force = p->f.f;

  // calculate position in cell, normalized to lattice size:
  FOR3D(k) {
    position[k]           = (p->r.p[k] - lparams.left_down_position[k]) * maggs.inva;
    left_down_position[k] = floor(position[k]);
    relative_position[k]  = position[k] - left_down_position[k];
    self_force[k] = 0.0;
    local_D_field[k] = 0.0;
  }

  // Copy permittivity values to the mini-lattice:

  for (iz=0;iz<zmax;iz++) {
    for (iy=0;iy<ymax;iy++) {
      for (ix=0;ix<xmax;ix++) {
	index = (iz + zmax*(iy + ymax*ix));
	globalindex = maggs_get_linear_index((left_down_position[0]+ix),
					     (left_down_position[1]+iy),
					     (left_down_position[2]+iz), lparams.dim);
	local_permittivity[index] = lattice[globalindex].permittivity[0];
	local_rho[index] = 0.0;
      }
    }
  }


  FOR3D(k) {
    self_force[k] = 0.0;
  }

  FOR3D(k) {
    force[k] += self_force[k];
  }
*/
}


/** For each particle P, calculates self energy influence
    with direct Greens function, assuming constant permittivity
    on each lattice site.
    @param P Particle pointer
*/
void maggs_calc_self_influence(Particle* P)
{
  int k;
  //	int ix, iy, iz;
  double invasq = maggs.inva*maggs.inva;
  double position[SPACE_DIM], left_down_position[SPACE_DIM], relative_position[SPACE_DIM];
  //	int index = 0;
  //	int globalindex = 0;
  double local_force[SPACE_DIM];
  double particle_charge = P->p.q;
	
  /* calculate position in cell, normalized to lattice size: */
  FOR3D(k) {
    position[k]           = (P->r.p[k] - lparams.left_down_position[k]) * maggs.inva;
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
    (left_down_position[2]+iz), lparams.dim);
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
    local_force[k] = pow(((0.5 - relative_position[k])*SELF_FACTOR_1), 3.0) +
      (0.5 - relative_position[k]) * SELF_FACTOR_2;
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

/** Calculate the actual force from the E-Field
    @param p      Particle pointer
    @param index  index of the lattice site
    @param grad   charge gradient
 */
void maggs_calc_part_link_forces(Particle *p, int index, double *grad)
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
      help_index[j] = maggs_get_offset(lattice[neighbor[index][j]].r[j], anchor_site->r[j], j, lparams.dim);
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
	local_force[j] += -grad[ind_grad]*Dfield[3*index+j];
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


/** Public function.
    Calculates the actual force on each particle
    by calling all other needed functions (except
    for maggs_propagate_B_field) */
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
  maggs_update_charge_gradients(grad);
	
  if(!init) {
    MAGGS_TRACE(fprintf(stderr, "running symplectic update\n"));
    maggs_couple_current_to_Dfield();
    maggs_add_transverse_field(time_step);  
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
	  pos[d]   = (p[i].r.p[d] - lparams.left_down_position[d])* maggs.inva;
	  first[d] = (int) pos[d];
	}
				
	index = maggs_get_linear_index(first[0],first[1],first[2],lparams.dim);
	maggs_calc_part_link_forces(&p[i], index, &grad[ip]);
	maggs_calc_self_influence(&p[i]);
	//      	printf("before: %f\n", p->f.f[0]);
	//	maggs_calc_interpolated_self_force(&p[i]);
	//	printf("after: %f\n", p->f.f[0]);
	ip+=12;
      }
    }
  }

  if (maggs.finite_epsilon_flag) maggs_compute_dipole_correction();

}








/********************************************/
/****** get energy and print out stuff ******/
/********************************************/

/** integrates 0.5*D*E over the whole system
    public function!
    @return returns electric energy
*/
double maggs_electric_energy()
{
  int x, y, z, i, k;
  double localresult = 0.;
  double globalresult = 0.;
	
  FORALL_INNER_SITES(x, y, z) {
    i = maggs_get_linear_index(x, y, z, lparams.dim);	  
    FOR3D(k){
      localresult += SQR(Dfield[i*3+k]) / lattice[i].permittivity[k];
    }
  }
  localresult *= 0.5*maggs.a;
  MPI_Allreduce(&localresult,&globalresult,1,MPI_DOUBLE,MPI_SUM,comm_cart);  
  return globalresult;
}


/** Public funxtion.
    Integrates the B-field over the whole system to get the
    energy of the magnetic field.
    @return returns magnetic energy
*/
double maggs_magnetic_energy()
{
  int x, y, z, i;
  double result = 0.;
  /*  double invmass = 1./maggs.f_mass; we have B^~=B*c !!!! */
	
  FORALL_INNER_SITES(x, y, z) {
    i = maggs_get_linear_index(x, y, z, lparams.dim);	  
    result += SQR(Bfield[i*3]) + SQR(Bfield[i*3+1]) + SQR(Bfield[i*3+2]);
  }
  /* B is rescaled !!! ATTENTION!!! */
  result *= 0.5*maggs.a;
  return result;
}

/***************************/
/****** init and exit ******/
/***************************/

/** Initialization function.
    Sets maggs structure variables.
    Calls calculation of initial D-field. */
void maggs_init()
{
    maggs.inva  = (double) maggs.mesh/box_l[0]; 
    maggs.a     = 1.0/maggs.inva;
    if (temperature>0.0) {
      maggs.prefactor  = sqrt(4. * M_PI * maggs.bjerrum * temperature);
      maggs.pref2      = maggs.bjerrum * temperature;
    } else {
      maggs.prefactor  = sqrt(4. * M_PI * maggs.bjerrum);
      maggs.pref2      = maggs.bjerrum;
    }			
    
    if (maggs_sanity_checks()) {
      maggs.bjerrum = 0.0;
      coulomb.bjerrum = 0.0;
      return;
    }
    		
    maggs_setup_local_lattice();
		
    /* enforce electric field onto the Born-Oppenheimer surface */
    maggs_calc_init_e_field();
    //    if(!this_node) fprintf(stderr, "%d: Electric field is initialized\n", this_node);
    maggs_calc_self_energy_coeffs();
}

/** Frees the dynamically allocated memory
    Currently not called from anywhere. */
void maggs_exit()
{
  free(lattice);
  free(neighbor);
  free(Dfield);
  free(Bfield);
}



#endif // ELECTROSTATICS
