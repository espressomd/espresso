// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2003; all rights reserved unless otherwise stated.

/** \file modes.c
    Implementation of \ref modes.h "modes.h"
*/

#include "modes.h"
#include "communication.h"
/** fftw plan for calculating the 2d mode analysis */
rfftwnd_plan mode_analysis_plan;
/** Flag to indicate when the grid size has changed*/
int mode_grid_changed = 1;

/** The full 3d grid for mode analysis */
int mode_grid_3d[3] = {0,0,0};
/** Integers labels for grid axes compared to real axes*/
int xdir = -1;
int ydir = -1;
int zdir = -1;




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
	fprintf(stderr,"%d, fft_modes_init: grid dimensions are <%d,%d,%d> one and only one must be = 0 \n",this_node, mode_grid_3d[0],mode_grid_3d[1],mode_grid_3d[2]);
	errexit();
      } else {
	zdir = i;
      }
    } 
    else if ( mode_grid_3d[i] < 0 ) {
      fprintf(stderr,"%d, fft_modes_init: grid dimensions are <%d,%d,%d> and all must be >= 0 \n",this_node, mode_grid_3d[0],mode_grid_3d[1],mode_grid_3d[2]);
	errexit();
    }
    else {
      if (  xdir == -1 ) {xdir = i;}
      else {ydir = i;}      
    }    
  }
  /* Check that grid normal was found */
  if ( zdir == -1 ) {
    fprintf(stderr,"%d, fft_modes_init: grid dimensions are <%d,%d,%d>. One and only one must be 0 \n",this_node, mode_grid_3d[0],mode_grid_3d[1],mode_grid_3d[2]);
    errexit();
  }
  STAT_TRACE(fprintf(stderr,"%d,map_to_2dgrid found the following mapping: xdir = %d, ydir = %d, zdir = %d \n",this_node, xdir, ydir, zdir);)
}

void fft_modes_init() {
  STAT_TRACE(fprintf(stderr,"%d,initializing fft for mode analysis \n",this_node);)
    if ( xdir + ydir + zdir == -3 ) {
      fprintf(stderr,"%d,attempt to perform mode analysis with uninitialized grid \n",this_node);
      errexit();
    }

  fprintf(stderr,"%d,destroying plan \n",this_node);
  rfftwnd_destroy_plan(mode_analysis_plan);


  mode_analysis_plan = rfftw2d_create_plan(mode_grid_3d[xdir], mode_grid_3d[ydir], FFTW_REAL_TO_COMPLEX,FFTW_MEASURE);
  STAT_TRACE(fprintf(stderr,"%d,created plan \n",this_node);)
  mode_grid_changed = 0;  
}

int modes2d(IntList *ptype, fftw_complex* modes) {
  int i,j, gi, gj;
  double grid_size[2];
  fftw_real* height_grid;

  STAT_TRACE(fprintf(stderr,"%d,executing modes2d \n",this_node);)
  if ( mode_grid_changed ) {    
    fft_modes_init();
  }

  grid_size[xdir] = box_l[xdir]/(double)mode_grid_3d[xdir];
  grid_size[ydir] = box_l[ydir]/(double)mode_grid_3d[ydir];
  STAT_TRACE(fprintf(stderr,"%d,grid mesh dimensions are: %f times %f \n",this_node,grid_size[xdir], grid_size[ydir]));    
  
  STAT_TRACE(fprintf(stderr,"%d, initializing height_grid \n",this_node));    

  height_grid = malloc((mode_grid_3d[xdir])*sizeof(fftw_complex)*mode_grid_3d[ydir]);
  for ( i = 0 ; i < mode_grid_3d[xdir] ; i++) {
    for ( j = 0 ; j < mode_grid_3d[ydir] ; j++) {
      //     printf("ihg %d, %d \n", i, j);
      height_grid[i+j*mode_grid_3d[ydir]] = 0;
    }
  }

  STAT_TRACE(fprintf(stderr,"%d, updating partconfig \n",this_node);)    
  updatePartCfg(WITHOUT_BONDS);


  STAT_TRACE(fprintf(stderr,"%d, creating height grid \n",this_node);)    
  for (i = 0 ; i < n_total_particles ; i++) {
    if ( !ptype || intlist_contains(ptype,partCfg[i].p.type)) {      
      fold_position(partCfg[i].r.p,partCfg[i].l.i);

      //      printf("i %d, %f, %f \n", i,partCfg[i].r.p[xdir], partCfg[i].r.p[ydir]);
      gi = floor( partCfg[i].r.p[xdir]/grid_size[xdir] );
      gj = floor( partCfg[i].r.p[ydir]/grid_size[ydir] ); 
      //      printf("hg %d, %d, %d \n", i, gi, gj);
      //      printf("hg %d, %d, %d, %f \n", i, gi, gj,height_grid[gi][gj]);
      // Use unfolded position for this
      unfold_position(partCfg[i].r.p,partCfg[i].l.i);
      height_grid[gi + gj*mode_grid_3d[ydir]] += partCfg[i].r.p[zdir];

    }
  }

  // Check array for zero values
  STAT_TRACE(fprintf(stderr,"%d,calling fftw \n",this_node));
  rfftwnd_one_real_to_complex(mode_analysis_plan, height_grid, modes);
  STAT_TRACE(fprintf(stderr,"%d,called fftw \n",this_node));

  free(height_grid);
  return 1;

}

