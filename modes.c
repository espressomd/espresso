// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2004; all rights reserved unless otherwise stated.

/** \file modes.c
    Implementation of \ref modes.h "modes.h"
*/

#include "modes.h"
#include "communication.h"
/** fftw plan for calculating the 2d mode analysis */
rfftwnd_plan mode_analysis_plan;
/** Flag to indicate when the grid size has changed*/
int mode_grid_changed = 1;

/** An array containing the size of the mode grid along each
    coordinate axis.  Note that at present this should be a symmetric
    square */
int mode_grid_3d[3] = {0,0,0};
/** Integer labels for grid axes compared to real axes*/
int xdir = -1;
int ydir = -1;
int zdir = -1;




/** Numerical tolerance to be used only in modes2d*/
#define MODES2D_NUM_TOL 0.00001

/** Default value for the distance beyond which a particle is said to
    have escaped the bilayer.  The default is set very large so that
    all lipids are counted*/
double stray_cut_off = 10000000.0;

/** This function takes a given grid a supplied by the user and
    determines the correct orientation in which to do the fourier
    transform.  In this regard one of the grid dimensions must be 0
    and the other two must be integer multiples of two and equal to
    each other.  The dimension that is 0 will be assigned an internal
    reference called zdir and will be used to calculate the height
    function used in the fft
*/
void map_to_2dgrid() {
  int i;
  STAT_TRACE(fprintf(stderr,"%d,executing map_to_2dgrid \n",this_node));
  /* Reset values of mapping */
  xdir = -1;
  ydir = -1;
  zdir = -1;

  /* Find the grid normal */
  for ( i = 0 ; i < 3 ; i++) {
    if ( mode_grid_3d[i] == 0 ) {
      if (zdir != -1 ) { /* grid normal must be unique */ 
	fprintf(stderr,
		"%d, fft_modes_init: grid dimensions are <%d,%d,%d> one and only one must be = 0 \n",
		this_node, mode_grid_3d[0],mode_grid_3d[1],mode_grid_3d[2]);
	errexit();
      } else {
	zdir = i;
      }
    } 
    else if ( mode_grid_3d[i] < 0 ) {
      fprintf(stderr,
	      "%d, fft_modes_init: grid dimensions are <%d,%d,%d> and all must be >= 0 \n",
	      this_node, mode_grid_3d[0],mode_grid_3d[1],mode_grid_3d[2]);
	errexit();
    }
    else {
      if (  xdir == -1 ) {xdir = i;}
      else {ydir = i;}      
    }    
  }
  /* Check that grid normal was found */
  if ( zdir == -1 ) {
    fprintf(stderr,
	    "%d, fft_modes_init: grid dimensions are <%d,%d,%d>. One and only one must be 0 \n",
	    this_node, mode_grid_3d[0],mode_grid_3d[1],mode_grid_3d[2]);
    errexit();
  }
  STAT_TRACE(fprintf(stderr,
		     "%d,map_to_2dgrid found the following mapping: xdir = %d, ydir = %d, zdir = %d \n",
		     this_node, xdir, ydir, zdir));


  /* Now that we know the grid normal check that the other two dimensions are equal and multiples of 2 */
  if ( mode_grid_3d[xdir] != mode_grid_3d[ydir] ) {
        fprintf(stderr,
		"%d, fft_modes_init: grid dimensions are <%d,%d,%d>. Two must be equal and the other 0 \n",
		this_node, mode_grid_3d[xdir],mode_grid_3d[ydir],mode_grid_3d[zdir]);
	errexit();
  }

  if ( (mode_grid_3d[xdir]/2.0 - floor(mode_grid_3d[xdir]/2.0) > MODES2D_NUM_TOL) 
       || (mode_grid_3d[ydir]/2.0 - floor(mode_grid_3d[ydir]/2.0) > MODES2D_NUM_TOL) ) {
    fprintf(stderr,
	    "%d, fft_modes_init: grid dimensions are <%d,%d,%d>. All non zero values must be integer multiples of 2 \n",
	    this_node, mode_grid_3d[xdir],mode_grid_3d[ydir],mode_grid_3d[zdir]);
    errexit();
  }


}

/** 
    Initialising routine for modes2d.  Called every time a change is
    made to the grid. This routine calculates the fftw plan for the
    subsequent fft and destroys any existing plan that might exist
*/
void fft_modes_init() {
  STAT_TRACE(fprintf(stderr,"%d,initializing fftw for mode analysis \n",this_node);)
    if ( xdir + ydir + zdir == -3 ) {
      fprintf(stderr,"%d,attempt to perform mode analysis with uninitialized grid \n",this_node);
      errexit();
    }

  STAT_TRACE(fprintf(stderr,"%d,destroying old fftw plan \n",this_node);)
  rfftwnd_destroy_plan(mode_analysis_plan);


  mode_analysis_plan = rfftw2d_create_plan(mode_grid_3d[xdir], mode_grid_3d[ydir], FFTW_REAL_TO_COMPLEX,FFTW_MEASURE);
  STAT_TRACE(fprintf(stderr,"%d,created new fftw plan \n",this_node);)
  mode_grid_changed = 0;  
}

/**
   This routine calculates the orientational order parameter for a
   lipid bilayer as defined in Brannigan and Brown 2004
*/
int orient_order(double* result)
{
  double dir[3];
  double sumdir[3] = {0,0,0};
  double* stored_dirs;
  double zref;
  int bilayer_cnt;
  int i,atom,tmpzdir;
  double dp;


  bilayer_cnt = 0;
  *result = 0;

  IntList l_orient;
  init_intlist(&l_orient);
  realloc_intlist(&l_orient, n_molecules);

  stored_dirs = malloc(sizeof(double)*n_molecules*3);

  if ( xdir + ydir + zdir == -3 ) {
    tmpzdir = 2;
  } else { 
    tmpzdir = zdir;
  };

  /* Update particles */
  updatePartCfg(WITHOUT_BONDS);
  /* Make sure particles are sorted */
  if (!sortPartCfg()) {
    fprintf(stderr,"%d,could not sort partCfg \n",this_node);
    errexit();
  }


  /* Find the mean z position and fold x y coordinates but not z*/
  zref = 0;
  for (i = 0 ; i < n_total_particles ; i++) {
    zref += partCfg[i].r.p[tmpzdir];
  }
  zref = zref/(double)(n_total_particles);
  

  for ( i = 0 ; i < n_molecules ; i++) {
    atom = topology[i].part.e[0];
    l_orient.e[i] = lipid_orientation(atom,partCfg,zref,dir);
    stored_dirs[i*3] = dir[0];
    stored_dirs[i*3+1] = dir[1];
    stored_dirs[i*3+2] = dir[2];

    if ( l_orient.e[i] == LIPID_UP ) {
      sumdir[0] += dir[0];
      sumdir[1] += dir[1];
      sumdir[2] += dir[2];
      bilayer_cnt++;
    }
    if ( l_orient.e[i] == LIPID_DOWN ) {
      sumdir[0] -= dir[0];
      sumdir[1] -= dir[1];
      sumdir[2] -= dir[2];
      bilayer_cnt++;
    }
  }



  for ( i = 0 ; i < 3 ; i++) {
    sumdir[i] = sumdir[i]/(double)(bilayer_cnt);
  }

  for ( i = 0 ; i < n_molecules ; i++ ) {
    dir[0] = stored_dirs[i*3];
    dir[1] = stored_dirs[i*3+1];
    dir[2] = stored_dirs[i*3+2];

    if ( l_orient.e[i] != LIPID_STRAY ) {
      dp = scalar(dir,sumdir);
      *result += dp*dp*1.5-0.5;      
    }

  }

  free(stored_dirs);
  realloc_intlist(&l_orient, 0);

  *result = *result/(double)(bilayer_cnt);
  return TCL_OK;

}





/** 
    This routine performs a simple check to see whether a lipid is
    oriented up or down or if it has escaped the bilayer. 

    \param id The particle identifier
    \param partCfg An array of sorted particles
    \param zref The average z position of all particles
 */
int lipid_orientation( int id, Particle* partCfg , double zref, double director[3]) {
  int mol_size, head_id, tail_id, mol_id, mol_type;
  int i, tmpzdir;
  double distance;


  if ( xdir + ydir + zdir == -3 ) {
    tmpzdir = 2;
  } else { 
    tmpzdir = zdir;
  };


  /* check molecule information */
  if ( n_molecules < 0 ) return (TCL_ERROR);
  mol_id = partCfg[id].p.mol_id ;
  mol_size = topology[mol_id].part.n;
  mol_type = topology[mol_id].type;

  /*
    printf("mol_type %d \n", mol_type);  
    printf("ids: %d %d %d N: %d parts: %d %d %d \n",id,partCfg[id].p.identity,mol_id,topology[mol_id].part.n,topology[mol_id].part.e[0],topology[mol_id].part.e[1],topology[mol_id].part.e[2]);
  fflush(stdout);
  */

  // Search the molecule for head and tail beads
  head_id = 0 ;
  tail_id = 0;
  for ( i = 0 ; i < mol_size ; i++ ) {
    /*    
	  printf("type: %d %d \n", partCfg[topology[mol_id].part.e[i]].p.type, LIPID_HEAD_TYPE);
	  fflush(stdout);
    */
    if ( partCfg[topology[mol_id].part.e[i]].p.type == LIPID_HEAD_TYPE ) {
      head_id = topology[mol_id].part.e[i];
      if ( i == 0 ) { 
	tail_id = topology[mol_id].part.e[mol_size-1]; 
      } else {
	tail_id = topology[mol_id].part.e[0];
      }
      break;
    }
  }

  /*  printf("head: %d, tail %d", head_id, tail_id); */

  /*
  if ( mol_type == 0 ) { 
    head_id = topology[mol_id].part.e[mol_size-1];
    tail_id = topology[mol_id].part.e[0];
    } else {
    tail_id = topology[mol_id].part.e[mol_size-1];
    head_id = topology[mol_id].part.e[0];
    }
  */

  distance = sqrt(pow((partCfg[tail_id].r.p[tmpzdir] - zref),2));
  /*  printf("zdir: %d \n", zdir);
      printf("head pos %f : tail pos %f \n", partCfg[head_id].r.p[tmpzdir], partCfg[tail_id].r.p[tmpzdir] );
      printf("dist %f : vect %f \n", distance,(partCfg[head_id].r.p[tmpzdir] - partCfg[tail_id].r.p[tmpzdir]) );
      fflush(stdout);
  */

  double len = 0;
  for ( i = 0 ; i < 3 ; i++ ) {
    director[i] = (partCfg[head_id].r.p[i] - 
		   partCfg[tail_id].r.p[i]);
    len += director[i]*director[i];
  }
  for ( i = 0 ; i < 3 ; i++ ) {
    director[i] = director[i]/sqrt(len);
  }
  


  if ( (partCfg[head_id].r.p[tmpzdir] - partCfg[tail_id].r.p[tmpzdir]) > 0.0 ) {
    /* Lipid is oriented up */
    if (  distance  > stray_cut_off ) {
      return LIPID_STRAY;
    } else {
      return LIPID_UP;
    }
  } else {
    if (  distance  > stray_cut_off ) {
      return LIPID_STRAY;
    } else {
      return LIPID_DOWN;
    }
  }

  return -1;
}

/** This routine performs must of the work involved in the analyze
    modes2d command.  A breakdown of what the routine does is as
    follows

    \li Calculates the average bead position in the height dimension \ref zdir

    \li Calculates the average height of each bilayer leaflet above or
    below \ref zdir separately.  These averages are then averaged
    together to create a height function over the 2d grid. In
    calculating this height function, checks are made for grid cells
    that contain no lipids. In such cases a value equal to the mean of
    surrounding cells is substituted.

    \li The height function is fourier transformed using the fftw library.

*/
int modes2d(fftw_complex* modes) {
  STAT_TRACE(fprintf(stderr,"%d,executing modes2d \n",this_node);)
  int i,j, gi, gj;
  double grid_size[2];
  double direction[3];
  fftw_real* height_grid;
  double* height_grid_up;
  double* height_grid_down;
  int* grid_parts;
  int* grid_parts_up;
  int* grid_parts_down;
  double zref;
  int nup;
  int ndown;
  int nstray;
  int l_orient;
  double norm;
  int xi, yi;
  double meanval ;
  int nonzerocnt, gapcnt;

  if ( mode_grid_changed ) {    
    fft_modes_init();
  }

  /* Allocate memory for height grid arrays and initialize these arrays */
  height_grid = malloc((mode_grid_3d[xdir])*sizeof(fftw_real)*mode_grid_3d[ydir]);
  height_grid_up = malloc((mode_grid_3d[xdir])*sizeof(double)*mode_grid_3d[ydir]);
  height_grid_down = malloc((mode_grid_3d[xdir])*sizeof(double)*mode_grid_3d[ydir]);
  grid_parts_up = malloc((mode_grid_3d[xdir])*sizeof(int)*mode_grid_3d[ydir]);
  grid_parts_down = malloc((mode_grid_3d[xdir])*sizeof(int)*mode_grid_3d[ydir]);
  grid_parts = malloc((mode_grid_3d[xdir])*sizeof(int)*mode_grid_3d[ydir]);
  for ( i = 0 ; i < mode_grid_3d[xdir] ; i++) {
    for ( j = 0 ; j < mode_grid_3d[ydir] ; j++) {
      height_grid[j+i*mode_grid_3d[xdir]] = 0;
      grid_parts[j+i*mode_grid_3d[xdir]] = 0;
      height_grid_up[j+i*mode_grid_3d[xdir]] = 0;
      grid_parts_up[j+i*mode_grid_3d[xdir]] = 0;
      height_grid_down[j+i*mode_grid_3d[xdir]] = 0;
      grid_parts_down[j+i*mode_grid_3d[xdir]] = 0;
    }
  }

  /* Calculate physical size of grid mesh */
  grid_size[xdir] = box_l[xdir]/(double)mode_grid_3d[xdir];
  grid_size[ydir] = box_l[ydir]/(double)mode_grid_3d[ydir];

  /* Update particles */
  updatePartCfg(WITHOUT_BONDS);
  //Make sure particles are sorted
    if (!sortPartCfg()) {
      fprintf(stderr,"%d,could not sort partCfg \n",this_node);
      return -1;
    }
  

  /* Find the mean z position and fold x y coordinates but not z*/
  zref = 0;
  for (i = 0 ; i < n_total_particles ; i++) {
    fold_coordinate(partCfg[i].r.p,partCfg[i].l.i,xdir);
    fold_coordinate(partCfg[i].r.p,partCfg[i].l.i,ydir);
    zref += partCfg[i].r.p[zdir];
  }
  zref = zref/(double)(n_total_particles);

  /* Calculate the non normalized height function of head lipids */
  nup = ndown = nstray = 0;
  for (i = 0 ; i < n_total_particles ; i++) {
    if ( (partCfg[i].p.type == 1)) {
      gi = floor( partCfg[i].r.p[xdir]/grid_size[xdir] );
      gj = floor( partCfg[i].r.p[ydir]/grid_size[ydir] );

      l_orient = lipid_orientation(i,partCfg,zref,direction);

      if ( l_orient != LIPID_STRAY ) {
	if ( l_orient == LIPID_UP ) {
	  nup++;
	  height_grid_up[gj + gi*mode_grid_3d[xdir]] += partCfg[i].r.p[zdir] - zref;
	  grid_parts_up[gj + gi*mode_grid_3d[xdir]] += 1;
	  //	  printf("up lipid \n");
	} else if ( l_orient == LIPID_DOWN ) {
	  ndown++;
	  height_grid_down[gj + gi*mode_grid_3d[xdir]] += partCfg[i].r.p[zdir] - zref;
	  grid_parts_down[gj + gi*mode_grid_3d[xdir]] += 1;
	  //  printf("down lipid \n");
	}
      } else {
	nstray++;
	//	printf("stray lipid \n");
      }
    }
    fflush(stdout);
    unfold_position(partCfg[i].r.p,partCfg[i].l.i);    
  }
  if ( nstray > 0 ) {
    printf("Warning: there were %d stray lipids in height calculation \n",nstray);
  }
  //  printf(" Lipids up = %d , Lipids down = %d \n",nup, ndown);

  STAT_TRACE(fprintf(stderr,"%d, Lipids up = %d , Lipids down = %d \n",this_node, nup, ndown));



  /*
  // Now for debugging purposes we impose a height function with known
  // power spectrum
  double qx = 0;
  double qy = 0;
  int qxi = 0;
  int qyi = 0;
  int NGrid = 8;
  double Norm = 1/8.0;
  //  int NGrid2 = 4;
  double BoxLength = 28.75;
  for ( i = 0 ; i < 8 ; i++) {
    for ( j = 0 ; j < 8 ; j++) {
      height_grid[j+i*8] = 0.0;
      for ( qxi = -4 ; qxi < 4 ; qxi++) {
	for ( qyi = -4 ; qyi < 4 ; qyi++) {
	  double dqx = qxi*2.0*3.141592/BoxLength;
	  double dqy = qyi*2.0*3.141592/BoxLength;
	  double dx = qxi*2.0*3.141592/(double)(NGrid);
	  double dy = qyi*2.0*3.141592/(double)(NGrid);
	  //	  int index = (qx+NGrid2)*NGrid + (qy+NGrid2);
	  height_grid[j+i*NGrid] += (1/(1 + dqx*dqx + dqy*dqy ))*cos((double)(dx*i + dy*j));
	}
      }
      height_grid[j+i*NGrid] *= Norm;
    }
  }
  
  // End debugging code
  */

  /* Norm we normalize the height function according the number of
     points in each grid cell */
  //  norm = 1.0/(double)(mode_grid_3d[xdir]);

  
  norm = 1.0;
  for ( i = 0 ; i < mode_grid_3d[xdir] ; i++) {
    for ( j = 0 ; j < mode_grid_3d[ydir] ; j++) {
      
      if ( ( grid_parts_up[j + i*mode_grid_3d[xdir]] > 0 ) && ( grid_parts_down[j + i*mode_grid_3d[xdir]] > 0 ) ) {
	height_grid[j+i*mode_grid_3d[xdir]] = 
	  0.5*norm*((height_grid_up[j+i*mode_grid_3d[xdir]])/(double)(grid_parts_up[j + i*mode_grid_3d[xdir]]) + 
		    (height_grid_down[j+i*mode_grid_3d[xdir]])/(double)(grid_parts_down[j + i*mode_grid_3d[xdir]]));
	grid_parts[j+i*mode_grid_3d[xdir]] = grid_parts_up[j + i*mode_grid_3d[xdir]] + grid_parts_down[j + i*mode_grid_3d[xdir]];
      } else {
	// Either upper or lower layer has no lipids
	height_grid[j+i*mode_grid_3d[xdir]] = 0.0;
	grid_parts[j+i*mode_grid_3d[xdir]] = 0;
      }
    }
  }
  

  /* Check height grid for zero values and substitute mean of surrounding cells */
  gapcnt = 0;
  for ( i = 0 ; i < mode_grid_3d[xdir] ; i++) {
    for ( j = 0 ; j < mode_grid_3d[ydir] ; j++) {
      if ( grid_parts[j + i*mode_grid_3d[xdir]] ==  0) {
	meanval = 0.0;
	nonzerocnt = 0;
	for ( xi = (i-1) ; xi <= (i+1) ; xi++) {
	  for ( yi = (j-1) ; yi <= (j+1) ; yi++) {
	    if ( height_grid[yi+xi*mode_grid_3d[xdir]] != 0 ) {
	      meanval += height_grid[yi+xi*mode_grid_3d[xdir]];
	      nonzerocnt++;
	    }
	  }
	}
	if ( nonzerocnt == 0 ) { 
	  fprintf(stderr,"Error: hole in membrane \n ");
	  fflush(stdout);
	  return -1;
	}
	gapcnt++;
      }      
    }
  }
  if ( gapcnt != 0 ) { 
    fprintf(stderr,"Warning: %d, gridpoints with no particles \n",gapcnt);
    fflush(stdout);
  }



  STAT_TRACE(fprintf(stderr,"%d,calling fftw \n",this_node));
  rfftwnd_one_real_to_complex(mode_analysis_plan, height_grid, modes);
  STAT_TRACE(fprintf(stderr,"%d,called fftw \n",this_node));


  free(height_grid);
  free(grid_parts);
  free(height_grid_up);
  free(height_grid_down);
  free(grid_parts_up);
  free(grid_parts_down);

  return 1;

}

#undef LIPID_UP 
#undef LIPID_DOWN 
#undef LIPID_STRAY 
#undef MODES2D_NUM_TOL

