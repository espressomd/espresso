/*
  Copyright (C) 2010,2011,2012,2013 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
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

/** \file modes.cpp
    Implementation of \ref modes.hpp "modes.h"
*/

#include "modes.hpp"
#include "communication.hpp"
#include "errorhandling.hpp"
#include "grid.hpp"

#ifdef MODES

/** fftw plan for calculating the 2d mode analysis */

#ifdef FFTW_ENABLE_FLOAT
typedef float fftw_real;
#else
typedef double fftw_real;
#endif

/** Flag to indicate when the grid size has changed*/
int mode_grid_changed = 1;

/** An array containing the size of the mode grid along each
    coordinate axis.  Note that at present this should be a symmetric
    square */
int mode_grid_3d[3] = {0,0,0};
/** Integer label x  for grid axes compared to real axes */
int xdir = -1;
/** Integer label y  for grid axes compared to real axes */
int ydir = -1;
/** Integer label z  for grid axes compared to real axes */
int zdir = -1;

/** Numerical tolerance to be used only in modes2d*/
#define MODES2D_NUM_TOL 0.00001

/** Default value for the distance beyond which a particle is said to
    have escaped the bilayer.  The default is set very large so that
    all lipids are counted*/
double stray_cut_off = 10000000.0;


void fold_all ( void ) {
  int i ;
  for (i = 0 ; i < n_part ; i++) {
    fold_coordinate(partCfg[i].r.p,partCfg[i].l.i,xdir);
    fold_coordinate(partCfg[i].r.p,partCfg[i].l.i,ydir);
    fold_coordinate(partCfg[i].r.p,partCfg[i].l.i,zdir);
  }
}

/* Simple helper function for calculating the reference
   zposition. Usually the bilayer midplane. Also fold the
   particles. Not exported */
double calc_zref ( int tmpzdir ) {
  float zref;
  int i;

  /* Fold all the particles */
  fold_all ();

 /* Find the mean z position of folded particles*/
  zref = 0;
  for (i = 0 ; i < n_part ; i++) {
    zref += partCfg[i].r.p[zdir];
  }
  zref = zref/(double)(n_part);
  return zref;
}


/** This function takes a given grid supplied by the user and
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
	char *errtxt = runtime_error(128 + 3*ES_INTEGER_SPACE);
	ERROR_SPRINTF(errtxt, "{029 fft_modes_init: grid dimensions are <%d,%d,%d>, but one and only one must be = 0} ",
		mode_grid_3d[0],mode_grid_3d[1],mode_grid_3d[2]);
	return;
      } else {
	zdir = i;
      }
    } 
    else if ( mode_grid_3d[i] < 0 ) {
      char *errtxt = runtime_error(128 + 3*ES_INTEGER_SPACE);
      ERROR_SPRINTF(errtxt, "{030 fft_modes_init: grid dimensions are <%d,%d,%d>, but all must be >= 0} ",
	      mode_grid_3d[0],mode_grid_3d[1],mode_grid_3d[2]);
      return;
    }
    else {
      if (  xdir == -1 ) {xdir = i;}
      else {ydir = i;}      
    }    
  }
  /* Check that grid normal was found */
  if ( zdir == -1 ) {
    char *errtxt = runtime_error(128 + 3*ES_INTEGER_SPACE);
    ERROR_SPRINTF(errtxt, "{031 fft_modes_init: grid dimensions are <%d,%d,%d>, but one and only one must be = 0} ",
	    mode_grid_3d[0],mode_grid_3d[1],mode_grid_3d[2]);
    return;
  }
  STAT_TRACE(fprintf(stderr,
		     "%d,map_to_2dgrid found the following mapping: xdir = %d, ydir = %d, zdir = %d \n",
		     this_node, xdir, ydir, zdir));


  /* Now that we know the grid normal check that the other two dimensions are equal and multiples of 2 */
  if ( mode_grid_3d[xdir] != mode_grid_3d[ydir] ) {
    char *errtxt = runtime_error(128 + 3*ES_INTEGER_SPACE);
    ERROR_SPRINTF(errtxt, "{032 fft_modes_init: grid dimensions are <%d,%d,%d>, but two must be equal and the other 0} ",
	    mode_grid_3d[xdir],mode_grid_3d[ydir],mode_grid_3d[zdir]);
    return;
  }

  if ( (mode_grid_3d[xdir]/2.0 - floor(mode_grid_3d[xdir]/2.0) > MODES2D_NUM_TOL) 
       || (mode_grid_3d[ydir]/2.0 - floor(mode_grid_3d[ydir]/2.0) > MODES2D_NUM_TOL) ) {
    char *errtxt = runtime_error(128 + 3*ES_INTEGER_SPACE);
    ERROR_SPRINTF(errtxt, "{033 fft_modes_init: grid dimensions are <%d,%d,%d>. All non zero values must be integer multiples of 2} ",
	    mode_grid_3d[xdir],mode_grid_3d[ydir],mode_grid_3d[zdir]);
    return;
  }
}



/**
   This routine calculates the orientational order parameter for a
   lipid bilayer.  
*/
int orient_order(double* result, double* stored_dirs)
{
  double dir[3];
  double refdir[3] = {0,0,1};
  double sumdir[3] = {0,0,0};
  double zref;
  int bilayer_cnt;
  int i,atom,tmpzdir;
  double dp;
  double len;

  IntList l_orient;
  init_intlist(&l_orient);
  realloc_intlist(&l_orient, n_molecules);

  bilayer_cnt = 0;
  *result = 0;


  // Check to see if the grid has been set and if not then interpret
  if ( xdir + ydir + zdir == -3 ) {
    tmpzdir = 2;
  } else { 
    tmpzdir = zdir;
  };

  /* Update particles */
  updatePartCfg(WITHOUT_BONDS);
  /* Make sure particles are sorted */
  if (!sortPartCfg()) {
    char *errtxt = runtime_error(128);
    ERROR_SPRINTF(errtxt, "{035 could not sort partCfg, particles have to start at 0 and have consecutive identities} ");
    return ES_ERROR;
  }

  /* Calculate the reference z position as its mean.*/
  zref = calc_zref( tmpzdir );
 

  /* Calculate the orientation of all the lipids in turn and include
   only UP or DOWN lipids in a calculation of the overall orientation
   direction of the bilayer .. ie the reference vector from which we
   can calculate the orientational order. */

  for ( i = 0 ; i < n_molecules ; i++) {
    atom = topology[i].part.e[0];
    l_orient.e[i] = lipid_orientation(atom,partCfg,zref,dir,refdir);
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

  /* Normalise the bilayer normal vector */
  len = 0.0;
  for ( i = 0 ; i < 3 ; i++) {
    sumdir[i] = sumdir[i]/(double)(bilayer_cnt);
    len += sumdir[i]*sumdir[i];
  }
  for ( i = 0 ; i < 3 ; i++) {
    sumdir[i] = sumdir[i]/sqrt(len);
  }

  /* Calculate the orientational order */
  for ( i = 0 ; i < n_molecules ; i++ ) {
    dir[0] = stored_dirs[i*3];
    dir[1] = stored_dirs[i*3+1];
    dir[2] = stored_dirs[i*3+2];

    if ( l_orient.e[i] != LIPID_STRAY && l_orient.e[i] != REAL_LIPID_STRAY ) {
      dp = scalar(dir,sumdir);
      *result += dp*dp*1.5-0.5;      
    }

  }

  
  realloc_intlist(&l_orient, 0);

  *result = *result/(double)(bilayer_cnt);
  return ES_OK;
}



/** 
    This routine performs a simple check to see whether a lipid is
    oriented up or down or if it has escaped the bilayer. 

    \param id The particle identifier
    \param partCfg An array of sorted particles

    \param zref The average z position of all particles. This is used
    to check for stray lipids

    \param director director
    \param refdir is a vector indicating the direction indicating
    up. This is usually the zaxis. If it is not the z axis then lipids
    will not be returned as stray.
 */
int lipid_orientation( int id, Particle* partCfg , double zref, double director[3], double refdir[3]) {
  int mol_size, head_id, tail_id, mol_id;
  int i;
  int tmpxdir,tmpydir,tmpzdir;
  double distance;
  double tailz;
  double len;
  double proj;

  if ( zdir + ydir + zdir != -3 ) {
    tmpxdir = xdir;
    tmpydir = ydir;
    tmpzdir = zdir;
  } else {
    tmpxdir = 0;
    tmpydir = 1;
    tmpzdir = 2;
  }


  /* check molecule information exists */
  if ( n_molecules < 0 ) return ES_ERROR;

  /* Get basic molecule parameters */
  mol_id = partCfg[id].p.mol_id ;
  mol_size = topology[mol_id].part.n;
 
  /* If the head and tail id's were not found above then assume the
     head atom is the first and tail is the last in the molecule */
  head_id = topology[mol_id].part.e[0];
  tail_id = topology[mol_id].part.e[mol_size -1];

  /* Check for stray lipids only if the system is a flat bilayer and the director is in the zdirection */
  if ( ( refdir[tmpxdir] == 0.0 ) && ( refdir[tmpydir] == 0 ) ) {
    tailz = partCfg[tail_id].r.p[tmpzdir];
    distance = sqrt(pow((tailz - zref),2));
    /* Check if the lipid is a stray based on this distance */
    if ( distance > stray_cut_off ) {
      return REAL_LIPID_STRAY;
    }
  }

  /* Calculate a normalised vector called director that is in the
     direction of the lipid from tail to head */
  len = 0;
  for ( i = 0 ; i < 3 ; i++ ) {
    director[i] = (partCfg[head_id].r.p[i] - 
		   partCfg[tail_id].r.p[i]);
    len += director[i]*director[i];
  }
  for ( i = 0 ; i < 3 ; i++ ) {
    director[i] = director[i]/sqrt(len);
  }
  

  proj = director[0]*refdir[0] + director[1]*refdir[1] + director[2]*refdir[2];
  if ( proj > 0.0 ) {
    /* Lipid is oriented up */
    return LIPID_UP;
  } else {
    return LIPID_DOWN;
  }
}

/* Get a complete list of the orientations of every lipid assuming a
   bilayer structure.  Requires grid*/
int get_lipid_orients(IntList* l_orient) {
  int i,gi,gj, atom;
  double zreflocal,zref;  
  double dir[3];
  double refdir[3] = {0,0,1};
  double grid_size[2];

  double* height_grid;

  if ( xdir + ydir + zdir == -3 || mode_grid_3d[xdir] <= 0 || mode_grid_3d[ydir] <= 0 ) {
    char *errtxt = runtime_error(128);
    ERROR_SPRINTF(errtxt,"{036 cannot calculate lipid orientations with uninitialized grid} ");
    return ES_ERROR;
  }

  /* Allocate memory for height grid arrays and initialize these arrays */
  height_grid = (double*) malloc((mode_grid_3d[xdir])*sizeof(double)*mode_grid_3d[ydir]);


  /* Calculate physical size of grid mesh */
  grid_size[xdir] = box_l[xdir]/(double)mode_grid_3d[xdir];
  grid_size[ydir] = box_l[ydir]/(double)mode_grid_3d[ydir];


  /* Update particles */
  updatePartCfg(WITHOUT_BONDS);
  //Make sure particles are sorted
  
  if (!sortPartCfg()) {
    char *errtxt = runtime_error(128);
    ERROR_SPRINTF(errtxt,"{094 could not sort partCfg} ");
    return -1;
  }
  if ( !calc_fluctuations(height_grid, 1) ) {
    char *errtxt = runtime_error(128);
    ERROR_SPRINTF(errtxt,"{034 calculation of height grid failed } ");
    return -1;
  }

  zref = calc_zref( zdir );

  for ( i = 0 ; i < n_molecules ; i++) {
    atom = topology[i].part.e[0];
    gi = floor( partCfg[atom].r.p[xdir]/grid_size[xdir] );
    gj = floor( partCfg[atom].r.p[ydir]/grid_size[ydir] );
    zreflocal = height_grid[gj+gi*mode_grid_3d[xdir]] + zref;
    l_orient->e[i] = lipid_orientation(atom,partCfg,zreflocal,dir,refdir);
  }

  free(height_grid);

  return 1;
}


/** This routine performs must of the work involved in the analyze
    modes2d command.  A breakdown of what the routine does is as
    follows

    \li fftw plans and in / out arrays are initialized as required

    \li calculate height function is called

    \li The height function is fourier transformed using the fftw library.

    Note: argument switch_fluc
    switch_fluc == 1 for height grid
    switch_fluc == 0 for thickness
*/
int modes2d(fftw_complex* modes, int switch_fluc) {
  /* All these variables need to be static so that the fftw3 plan can
     be initialised and reused */
  static  fftw_plan mode_analysis_plan; // height grid
  /** Input values for the fft */
  static  double* height_grid;
  /** Output values for the fft */
  static  fftw_complex* result;

/** 
    Every time a change is made to the grid calculate the fftw plan
    for the subsequent fft and destroy any existing plans
*/
  if ( mode_grid_changed ) {
    STAT_TRACE(fprintf(stderr,"%d,initializing fftw for mode analysis \n",this_node));
    if ( xdir + ydir + zdir == -3 ) {
      char *errtxt = runtime_error(128);
      ERROR_SPRINTF(errtxt,"{092 attempt to perform mode analysis with uninitialized grid} ");
      return -1;
    }

    STAT_TRACE(fprintf(stderr,"%d,destroying old fftw plan \n",this_node));

    /* Make sure all memory is free and old plan is destroyed. It's ok
       to call these functions on uninitialised pointers I think */
    fftw_free(result);
    fftw_free(height_grid);
    fftw_destroy_plan(mode_analysis_plan);
    fftw_cleanup(); 
    /* Allocate memory for input and output arrays */
    height_grid = (double*) malloc((mode_grid_3d[xdir])*sizeof(double)*mode_grid_3d[ydir]);
    result      = (fftw_complex*) malloc((mode_grid_3d[ydir]/2+1)*(mode_grid_3d[xdir])*sizeof(fftw_complex)); 
    mode_analysis_plan = fftw_plan_dft_r2c_2d(mode_grid_3d[xdir],mode_grid_3d[ydir],height_grid, result,FFTW_ESTIMATE);

    STAT_TRACE(fprintf(stderr,"%d,created new fftw plan \n",this_node));
    mode_grid_changed = 0;  
    
  }

  /* Update particles */
  updatePartCfg(WITHOUT_BONDS);
  //Make sure particles are sorted
  
  if (!sortPartCfg()) {
    char *errtxt = runtime_error(128);
    ERROR_SPRINTF(errtxt,"{094 could not sort partCfg} ");
    return -1;
  }
  if ( !calc_fluctuations(height_grid, switch_fluc)) {
    char *errtxt = runtime_error(128);
    ERROR_SPRINTF(errtxt,"{034 calculation of height grid failed } ");
    return -1;
  }

  STAT_TRACE(fprintf(stderr,"%d,calling fftw \n",this_node));

  fftw_execute(mode_analysis_plan);
  /* Copy result to modes */
  memcpy(modes, result, mode_grid_3d[xdir]*(mode_grid_3d[ydir]/2 + 1)*sizeof(fftw_complex));
  
  
  STAT_TRACE(fprintf(stderr,"%d,called fftw \n",this_node));    
    
  return 1;
    
}


/** This routine calculates density profiles for given bead types as a
    function of height relative to the bilayer midplane where the
    bilayer is assumed to wrap around a sphere of radius r and center
    cx,cy,cz.

    \li The radius value of each bead is then calculated relative to
    the colloid radius and the population of the relevant bin
    is increased.
**/
int bilayer_density_profile_sphere (IntList *beadids, double rrange , DoubleList *density_profile, double radius, double center[3]) {
  int i,j;
  int thisbin,nbins;
  double binwidth;
  int nbeadtypes,l_orient;
  double relativeradius;
  double rvec[3];
  double direction[3]; 
  double innerr;
  double outerr;
  double binvolume;
  double piconst;
  double rs;


  nbins = density_profile[0].max;
  nbeadtypes=beadids->max;
  binwidth = 2*rrange/(double)nbins;


  /* Update particles */
  updatePartCfg(WITHOUT_BONDS);
  //Make sure particles are sorted
  if (!sortPartCfg()) {
    char *errtxt = runtime_error(128);
    ERROR_SPRINTF(errtxt,"{094 could not sort partCfg} ");
    return -1;
  }

  if ( density_profile == NULL ) {
    char *errtxt = runtime_error(128);
    ERROR_SPRINTF(errtxt,"{095 density_profile not initialized in calc_bilayer_density_profile } ");
    return -1;
  }

  // Do this to fold the particles
  fold_all( );


   for (i = 0 ; i < n_part ; i++) {
    for ( j = 0 ; j < nbeadtypes ; j++ ) {
      if ( beadids->e[j] == partCfg[i].p.type ) {
	/* What is the relative height compared to the grid */
	get_mi_vector(rvec,partCfg[i].r.p,center);
	relativeradius = sqrt(sqrlen(rvec)) - radius;

	/* If the particle is within our zrange then add it to the profile */
	if ( ( -rrange < relativeradius) && ( relativeradius < rrange) ) {
	  thisbin = (int)(floor((relativeradius+rrange)/binwidth));
	  if ( thisbin < 0 || thisbin >= density_profile[j].max ) {
	        char *errtxt = runtime_error(128);
		ERROR_SPRINTF(errtxt,"{095 bin is outside range } ");
		return -1;
	  }

	  l_orient = lipid_orientation(i,partCfg,0.0,direction,rvec);
	  /* Distinguish between lipids that are in the top and bottom layers */
	  if ( l_orient == LIPID_UP ) {
	    density_profile[j].e[thisbin] += 1.0;
	  }
	  if ( l_orient == LIPID_DOWN ) {
	    density_profile[2*nbeadtypes-j-1].e[thisbin] += 1.0;	    
	  }
	}
	    
      }
    }
  }
  /* Normalize the density profile */
  piconst = (4.0*3.141592)/(3.0);
  rs = radius - rrange;
  for ( i = 0 ; i < 2*nbeadtypes ; i++ ) {
    for ( j = 0 ; j < nbins ; j++ ) {

      innerr = j*binwidth+rs;
      outerr = (j+1)*binwidth + rs;
      binvolume = piconst*(outerr*outerr*outerr - innerr*innerr*innerr);
      density_profile[i].e[j] = density_profile[i].e[j]/binvolume;
    }
  }

   return 1;

}


/** This routine calculates density profiles for given bead types as a
    function of height relative to the bilayer midplane.

    \li First the height function is calculated

    \li The height value of each bead is then calculated relative to
    the overall height function and the population of the relevant bin
    is increased.

**/
int bilayer_density_profile ( IntList *beadids, double hrange , DoubleList *density_profile, int usegrid) {
  int i,j, gi,gj;
  double* tmp_height_grid, zref,zreflocal;
  int thisbin,nbins;
  double grid_size[2];
  double binwidth;
  int nbeadtypes,l_orient;
  double relativeheight;
  double direction[3];  
  double refdir[3] = {0,0,1};
  int tmpzdir;
  double binvolume;
  nbins = density_profile[0].max;
  nbeadtypes=beadids->max;

  /*        Check to see that there is a mode grid to work with    */
  if ( xdir + ydir + zdir == -3 ) {
    char *errtxt = runtime_error(128);
    ERROR_SPRINTF(errtxt,"{092 attempt to perform mode analysis with uninitialized grid} ");
    return -1;
  }
  
  /* Allocate memory for the grid if we are going to need it */
  tmp_height_grid = (double*) malloc((mode_grid_3d[xdir])*sizeof(double)*mode_grid_3d[ydir]);
  /* Calculate physical size of grid mesh */
  grid_size[xdir] = box_l[xdir]/(double)mode_grid_3d[xdir];
  grid_size[ydir] = box_l[ydir]/(double)mode_grid_3d[ydir];

  
  /* Calculate the height grid which also ensures that particle config is updated */
  if ( !calc_fluctuations(tmp_height_grid, 1) ) {
    char *errtxt = runtime_error(128);
    ERROR_SPRINTF(errtxt,"{034 calculation of height grid failed } ");
    return -1;
  } 
  
  tmpzdir = zdir;


  binwidth = hrange*2.0/(double)(nbins);

  if ( density_profile == NULL ) {
    char *errtxt = runtime_error(128);
    ERROR_SPRINTF(errtxt,"{095 density_profile not initialized in calc_bilayer_density_profile } ");
    return -1;
  }

  zref = calc_zref( tmpzdir );

  for (i = 0 ; i < n_part ; i++) {
    for ( j = 0 ; j < nbeadtypes ; j++ ) {
      if ( beadids->e[j] == partCfg[i].p.type ) {

	if ( usegrid ) {
	  /* Where are we on the height grid */
	  gi = (int)(floor( partCfg[i].r.p[xdir]/grid_size[xdir] ));
	  gj = (int)(floor( partCfg[i].r.p[ydir]/grid_size[ydir] ));
	  zreflocal = tmp_height_grid[gj + gi*mode_grid_3d[xdir]] + zref;
	} else {
	  zreflocal = zref;
	}

	/* What is the relative height compared to the grid */
	relativeheight = partCfg[i].r.p[tmpzdir] - zreflocal;							       
	/* If the particle is within our zrange then add it to the profile */
	if ( (relativeheight*relativeheight - hrange*hrange) <= 0 ) {
	  thisbin = (int)(floor((relativeheight + hrange)/binwidth));
	  l_orient = lipid_orientation(i,partCfg,zreflocal,direction,refdir);
	  /* Distinguish between lipids that are in the top and bottom layers */
	  if ( l_orient == LIPID_UP ) {
	    density_profile[j].e[thisbin] += 1.0;
	  }
	  if ( l_orient == LIPID_DOWN ) {
	    density_profile[2*nbeadtypes-j-1].e[thisbin] += 1.0;	    
	  }
	}	    
      }
    }
  }

  /* Normalize the density profile */
  binvolume = binwidth*box_l[xdir]*box_l[ydir];
  for ( i = 0 ; i < 2*nbeadtypes ; i++ ) {
    for ( j = 0 ; j < nbins ; j++ ) {
      density_profile[i].e[j] = density_profile[i].e[j]/binvolume;
    }
  }


  free(tmp_height_grid);

  return 1;
}

/**
   This routine calculates an average height for the bilayer in each
   grid square as follows;

    \li Calculates the average bead position in the height dimension
    \ref modes::zdir

    \li Calculates the average height of each bilayer leaflet above or
    below \ref modes::zdir separately.  These averages are then averaged
    together to create a height function over the 2d grid. In
    calculating this height function, checks are made for grid cells
    that contain no lipids. In such cases a value equal to the mean of
    surrounding cells is substituted.

    \li Also can calculate the thickness of a flat bilayer.
*/



int calc_fluctuations ( double* height_grid, int switch_fluc ) {
  if (switch_fluc == 1){
    STAT_TRACE(fprintf(stderr,"%d,calculating height grid \n",this_node));
  } else { if (switch_fluc == 0) {
    STAT_TRACE(fprintf(stderr,"%d,calculating thickness \n",this_node));
    } else {
       char *errtxt = runtime_error(128);
       ERROR_SPRINTF(errtxt,"{097 Wrong argument in calc_fluctuations function} ");
       return -1;
    }
  }
    
  int i, j, gi, gj;
  double grid_size[2];
  double direction[3];  
  double refdir[3] = {0,0,1};
  double* height_grid_up;
  double* height_grid_down;
  int* grid_parts;
  int* grid_parts_up;
  int* grid_parts_down;
  double zreflocal, zref;
  int nup;
  int ndown;
  int nstray;
  int nrealstray;
  int l_orient;
  double norm;
  int xi, yi;
  double meanval;
  int nonzerocnt, gapcnt;






  if ( xdir + ydir + zdir == -3 ) {
      char *errtxt = runtime_error(128);
      ERROR_SPRINTF(errtxt,"{092 attempt to calculate height grid / thickness with uninitialized grid} ");
      return -1;
  }
  

  if ( height_grid == NULL ) {
    char *errtxt = runtime_error(128);
    ERROR_SPRINTF(errtxt,"{093 you must allocate memory for the height grid / thickness first} ");
    return -1;
  }



  /* Allocate memory for height grid / thickness arrays and initialize these arrays */

  height_grid_up = (double*) malloc((mode_grid_3d[xdir])*sizeof(double)*mode_grid_3d[ydir]);
  height_grid_down = (double*) malloc((mode_grid_3d[xdir])*sizeof(double)*mode_grid_3d[ydir]);
  grid_parts_up = (int*) malloc((mode_grid_3d[xdir])*sizeof(int)*mode_grid_3d[ydir]);
  grid_parts_down = (int*) malloc((mode_grid_3d[xdir])*sizeof(int)*mode_grid_3d[ydir]);
  grid_parts = (int*) malloc((mode_grid_3d[xdir])*sizeof(int)*mode_grid_3d[ydir]);
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
    char *errtxt = runtime_error(128);
    ERROR_SPRINTF(errtxt,"{094 could not sort partCfg} ");
    return -1;
  }
  
  
  
  /* Find the mean z position of folded coordinates*/ 
  zref = calc_zref( zdir );
  
  /* Calculate an initial height function of all particles */
  for (i = 0 ; i < n_part ; i++) {
    gi = floor( partCfg[i].r.p[xdir]/grid_size[xdir] );
    gj = floor( partCfg[i].r.p[ydir]/grid_size[ydir] );
    height_grid[gj + gi*mode_grid_3d[xdir]] += partCfg[i].r.p[zdir];
    grid_parts[gj + gi*mode_grid_3d[xdir]] += 1;
  }
  
  /* Normalise the initial height function */
  for ( i = 0 ; i < mode_grid_3d[xdir] ; i++) {
    for ( j = 0 ; j < mode_grid_3d[ydir] ; j++) {
      if ( grid_parts[j+i*mode_grid_3d[xdir]] > 0 ) {
	height_grid[j+i*mode_grid_3d[xdir]] = height_grid[j+i*mode_grid_3d[xdir]]/(double)(grid_parts[j+i*mode_grid_3d[xdir]]);
      } else {
	height_grid[j+i*mode_grid_3d[xdir]] = zref;
      }
    }
  }
  
  /* We now use this initial height function to calculate the
     lipid_orientation and thereby calculate populations in upper and
     lower leaflets */
  
  
  /* Calculate the non normalized height function based on all lipids */
  nup = ndown = nstray = nrealstray = 0;
  for (i = 0 ; i < n_part ; i++) {
    gi = floor( partCfg[i].r.p[xdir]/grid_size[xdir] );
    gj = floor( partCfg[i].r.p[ydir]/grid_size[ydir] );
    
    zreflocal = height_grid[gj+gi*mode_grid_3d[xdir]];
    
    l_orient = lipid_orientation(i,partCfg,zreflocal,direction,refdir);
    
    if ( l_orient != LIPID_STRAY && l_orient != REAL_LIPID_STRAY) {
      if ( l_orient == LIPID_UP ) {
	nup++;
	height_grid_up[gj + gi*mode_grid_3d[xdir]] += partCfg[i].r.p[zdir] - zref;
	grid_parts_up[gj + gi*mode_grid_3d[xdir]] += 1;
      } else if ( l_orient == LIPID_DOWN ) {
	ndown++;
	height_grid_down[gj + gi*mode_grid_3d[xdir]] += partCfg[i].r.p[zdir] - zref;
	grid_parts_down[gj + gi*mode_grid_3d[xdir]] += 1;
      }
    } else {
      if ( l_orient == REAL_LIPID_STRAY ) {
	nrealstray++;
      }
    }
  }



  /*
    if ( nrealstray > 0 || nstray > 0) {
    printf("Warning: there were %d stray lipid particles in height calculation \n", nrealstray);
    }
    printf(" Lipid particles up = %d , Lipid particles down = %d \n",nup, ndown); */
  
  STAT_TRACE(fprintf(stderr,"%d, Lipids up = %d , Lipids down = %d \n",this_node, nup, ndown));
  
  /* Reinitialise the height grid */
  for ( i = 0 ; i < mode_grid_3d[xdir] ; i++) {
    for ( j = 0 ; j < mode_grid_3d[ydir] ; j++) {
      height_grid[j+i*mode_grid_3d[xdir]] = 0.0;
      grid_parts[j+i*mode_grid_3d[xdir]] = 0;
    }
  }
  
  /* Norm we normalize the height function according the number of
     points in each grid cell */
  norm = 1.0;
  for ( i = 0 ; i < mode_grid_3d[xdir] ; i++) {
    for ( j = 0 ; j < mode_grid_3d[ydir] ; j++) {
      
      if ( ( grid_parts_up[j + i*mode_grid_3d[xdir]] > 0 ) && ( grid_parts_down[j + i*mode_grid_3d[xdir]] > 0 ) ) {
	/* 
	   This is where we distinguish between height_grid and thickness:
	   h = .5*(h_up + h_down)
	   t = h_up - h_down
	*/
	if (switch_fluc == 1)
	  height_grid[j+i*mode_grid_3d[xdir]] = 
	    0.5*norm*((height_grid_up[j+i*mode_grid_3d[xdir]])/(double)(grid_parts_up[j + i*mode_grid_3d[xdir]]) + 
		      (height_grid_down[j+i*mode_grid_3d[xdir]])/(double)(grid_parts_down[j + i*mode_grid_3d[xdir]]));
	else
	  height_grid[j+i*mode_grid_3d[xdir]] = 
	    norm*((height_grid_up[j+i*mode_grid_3d[xdir]])/(double)(grid_parts_up[j + i*mode_grid_3d[xdir]]) - 
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
	  char *errtxt = runtime_error(128);
	  ERROR_SPRINTF(errtxt,"{095 hole in membrane} ");
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
  
  free(grid_parts);
  free(height_grid_up);
  free(height_grid_down);
  free(grid_parts_up);
  free(grid_parts_down);
  
  return 1;
  
}


#undef MODES2D_NUM_TOL

#endif
