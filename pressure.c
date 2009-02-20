// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2006; all rights reserved unless otherwise stated.
/** \file pressure.c
    Implementation of \ref pressure.h "pressure.h".
*/
#include "pressure.h"
#include "parser.h"
#include "cells.h"
#include "integrate.h"
#include "initialize.h"
#include "domain_decomposition.h"
#include "nsquare.h"
#include "layered.h"

Observable_stat virials  = {0, {NULL,0,0}, 0,0,0,0};
Observable_stat total_pressure = {0, {NULL,0,0}, 0,0,0,0};
Observable_stat p_tensor = {0, {NULL,0,0},0,0,0,0};
Observable_stat total_p_tensor = {0, {NULL,0,0},0,0,0,0};

/* Observables used in the calculation of intra- and inter- molecular non-bonded contributions
   to pressure and to stress tensor */
Observable_stat_non_bonded virials_non_bonded  = {0, {NULL,0,0}, 0,0,0};
Observable_stat_non_bonded total_pressure_non_bonded = {0, {NULL,0,0}, 0,0,0};
Observable_stat_non_bonded p_tensor_non_bonded = {0, {NULL,0,0},0,0,0};
Observable_stat_non_bonded total_p_tensor_non_bonded = {0, {NULL,0,0},0,0,0};

nptiso_struct   nptiso   = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,{0.0,0.0,0.0},{0.0,0.0,0.0},1, 0 ,{NPTGEOM_XDIR, NPTGEOM_YDIR, NPTGEOM_ZDIR},0,0,0};

/************************************************************/
/* callbacks for setmd                                      */
/************************************************************/

int piston_callback(Tcl_Interp *interp, void *_data) {
  double data = *(double *)_data;
  if (data < 0.0) { Tcl_AppendResult(interp, "the piston's mass must be positive.", (char *) NULL); return (TCL_ERROR); }
  nptiso.piston = data;
  mpi_bcast_parameter(FIELD_NPTISO_PISTON);
  return (TCL_OK);
}

int p_ext_callback(Tcl_Interp *interp, void *_data) {
  double data = *(double *)_data;
  nptiso.p_ext = data;
  mpi_bcast_parameter(FIELD_NPTISO_PEXT);
  return (TCL_OK);
}

int p_diff_callback(Tcl_Interp *interp, void *_data) {
  double data = *(double *)_data;
  nptiso.p_diff = data;
  mpi_bcast_parameter(FIELD_NPTISO_PDIFF);
  return (TCL_OK);
}


/************************************************************/
/* local prototypes                                         */
/************************************************************/

/** Calculate long range virials (P3M, MMM2d...). */
void calc_long_range_virials();

/** Initializes a virials Observable stat. */
void init_virials(Observable_stat *stat);

/** Initializes a virials Observable stat. */
void init_virials_non_bonded(Observable_stat_non_bonded *stat_nb);

/** on the master node: calc energies only if necessary
    @param v_comp flag which enables (1) compensation of the velocities required 
		  for deriving a pressure reflecting \ref nptiso_struct::p_inst
		  (hence it only works with domain decomposition); naturally it
		  therefore doesn't make sense to use it without NpT. */
void master_pressure_calc(int v_comp);

/** Initializes stat to be used by \ref calc_p_tensor. */
void init_p_tensor(Observable_stat *stat);

/** Initializes stat_nb to be used by \ref calc_p_tensor. */
void init_p_tensor_non_bonded(Observable_stat_non_bonded *stat_nb);

/*********************************/
/* Scalar and Tensorial Pressure */
/*********************************/

void pressure_calc(double *result, double *result_t, double *result_nb, double *result_t_nb, int v_comp)
{
  int n, i;
  double volume = box_l[0]*box_l[1]*box_l[2];

  if (!check_obs_calc_initialized())
    return;

  init_virials(&virials);

  init_p_tensor(&p_tensor);
  
  init_virials_non_bonded(&virials_non_bonded);

  init_p_tensor_non_bonded(&p_tensor_non_bonded);

  on_observable_calc();

  switch (cell_structure.type) {
  case CELL_STRUCTURE_LAYERED:
    layered_calculate_virials();
    break;
  case CELL_STRUCTURE_DOMDEC:
    if (rebuild_verletlist) build_verlet_lists();
    calculate_verlet_virials(v_comp);
    break;
  case CELL_STRUCTURE_NSQUARE:
    nsq_calculate_virials();
  }
  /* rescale kinetic energy (=ideal contribution) */
#ifdef ROTATION
  virials.data.e[0] /= (6.0*volume*time_step*time_step);
#else
  virials.data.e[0] /= (3.0*volume*time_step*time_step);
#endif

  calc_long_range_virials();

  for (n = 1; n < virials.data.n; n++)
    virials.data.e[n] /= 3.0*volume;

    
 /* stress tensor part */
 /* ROTATION option does not effect stress tensor calculations since rotational
    energy is not included in the ideal term (unlike for the pressure) */
  for(i=0; i<9; i++)
    p_tensor.data.e[i] /= (volume*time_step*time_step);
  
  for(i=9; i<p_tensor.data.n; i++)
    p_tensor.data.e[i]  /= volume;
  
  /* Intra- and Inter- part of nonbonded interaction */
  for (n = 0; n < virials_non_bonded.data_nb.n; n++)
    virials_non_bonded.data_nb.e[n] /= 3.0*volume;
  
  for(i=0; i<p_tensor_non_bonded.data_nb.n; i++)
    p_tensor_non_bonded.data_nb.e[i]  /= volume;
  
  /* gather data */
  MPI_Reduce(virials.data.e, result, virials.data.n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(p_tensor.data.e, result_t, p_tensor.data.n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  
  MPI_Reduce(virials_non_bonded.data_nb.e, result_nb, virials_non_bonded.data_nb.n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(p_tensor_non_bonded.data_nb.e, result_t_nb, p_tensor_non_bonded.data_nb.n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); 
}

/************************************************************/

void calc_long_range_virials()
{
#ifdef ELECTROSTATICS
  /* calculate k-space part of electrostatic interaction. */
  switch (coulomb.method) {
#ifdef ELP3M
  case COULOMB_ELC_P3M:
    fprintf(stderr, "WARNING: pressure calculated, but ELC pressure not implemented\n");
    break;
  case COULOMB_P3M: {
    int k;
    P3M_charge_assign();
    virials.coulomb[1] = P3M_calc_kspace_forces_for_charges(0,1);
    
    for(k=0;k<3;k++)
      p_tensor.coulomb[9+ k*3 + k] = virials.coulomb[1]/3.;
    fprintf(stderr, "WARNING: stress tensor calculated, but P3M stress tensor not implemented\n");
    break;
  }
#endif
  case COULOMB_MMM2D:
    fprintf(stderr, "WARNING: pressure calculated, but MMM2D pressure not implemented\n");
    break;
  case COULOMB_MMM1D:
    fprintf(stderr, "WARNING: pressure calculated, but MMM1D pressure not implemented\n");
    break;
  }
#endif /*ifdef ELECTROSTATICS */  
  
#ifdef MAGNETOSTATICS
  /* calculate k-space part of magnetostatic interaction. */
  switch (coulomb.Dmethod) {
#ifdef DAWAANR 
     case DIPOLAR_ALL_WITH_ALL_AND_NO_REPLICA:
    fprintf(stderr, "WARNING: pressure calculated, but DAWAANR pressure not implemented\n");
    break;
#endif 
#ifdef MAGNETIC_DIPOLAR_DIRECT_SUM
#ifdef MDLC
 case DIPOLAR_MDLC_DS:
    fprintf(stderr, "WARNING: pressure calculated, but DLC pressure not implemented\n");
    break;
#endif    
    case DIPOLAR_DS:
    fprintf(stderr, "WARNING: pressure calculated, but  MAGNETIC DIRECT SUM pressure not implemented\n");
    break;
#endif 
   
#ifdef ELP3M
#ifdef MDLC
  case DIPOLAR_MDLC_P3M:
    fprintf(stderr, "WARNING: pressure calculated, but DLC pressure not implemented\n");
    break;
#endif    
  case DIPOLAR_P3M: {
    int k;
    P3M_dipole_assign();
    virials.dipolar[1] = P3M_calc_kspace_forces_for_dipoles(0,1);
     
    for(k=0;k<3;k++)
      p_tensor.coulomb[9+ k*3 + k] = virials.dipolar[1]/3.;
    fprintf(stderr, "WARNING: stress tensor calculated, but dipolar P3M stress tensor not implemented\n");
    fprintf(stderr, "WARNING: things have been added to the coulomb virial and p_tensor arrays !!!!!!!\n");
    
    break;
  }
#endif
 } 
#endif /*ifdef MAGNETOSTATICS */
}

/* Initialize the virials used in the calculation of the scalar pressure */
/************************************************************/
void init_virials(Observable_stat *stat)
{
    int n_pre, n_non_bonded, n_coulomb, n_dipolar;

  n_pre        = 1;
  n_non_bonded = (n_particle_types*(n_particle_types+1))/2;

  n_coulomb    = 0;
  n_dipolar    = 0;
#ifdef ELECTROSTATICS
  switch (coulomb.method) {
  case COULOMB_NONE: n_coulomb = 0; break;
  case COULOMB_P3M:  n_coulomb = 2; break;
  default: n_coulomb  = 1;
  }
#endif
#ifdef MAGNETOSTATICS
  switch (coulomb.Dmethod) {
  case DIPOLAR_NONE:  n_dipolar = 0; break;
  case DIPOLAR_ALL_WITH_ALL_AND_NO_REPLICA:  n_dipolar = 0; break;
  case DIPOLAR_DS:  n_dipolar = 0; break;
  case DIPOLAR_P3M:   n_dipolar = 2; break;
  }
#endif


  
  obsstat_realloc_and_clear(stat, n_pre, n_bonded_ia, n_non_bonded, n_coulomb, n_dipolar, 1);
  stat->init_status = 0;
}

/************************************************************/
void init_virials_non_bonded(Observable_stat_non_bonded *stat_nb)
{
  int n_non_bonded;

  n_non_bonded = (n_particle_types*(n_particle_types+1))/2;

  obsstat_realloc_and_clear_non_bonded(stat_nb, n_non_bonded, 1);
  stat_nb->init_status_nb = 0;
}

/* Initialize the p_tensor */
/***************************/
void init_p_tensor(Observable_stat *stat)
{
    int n_pre, n_non_bonded, n_coulomb, n_dipolar;

  n_pre        = 1;
  n_non_bonded = (n_particle_types*(n_particle_types+1))/2;

  n_coulomb = 0;
  n_dipolar = 0;

#ifdef ELECTROSTATICS
  switch (coulomb.method) {
  case COULOMB_NONE: n_coulomb = 0; break;
  case COULOMB_P3M:  n_coulomb = 2; break;
  default: n_coulomb  = 1;
  }
#endif
  
#ifdef MAGNETOSTATICS
  switch (coulomb.Dmethod) {
  case DIPOLAR_NONE: n_dipolar = 0; break;
  case DIPOLAR_ALL_WITH_ALL_AND_NO_REPLICA:  n_dipolar = 0; break;
  case DIPOLAR_DS:  n_dipolar = 0; break;
  case DIPOLAR_P3M:  n_dipolar = 2; break;
  }
#endif

  obsstat_realloc_and_clear(stat, n_pre, n_bonded_ia, n_non_bonded, n_coulomb, n_dipolar, 9);
  stat->init_status = 0;
}

/***************************/
void init_p_tensor_non_bonded(Observable_stat_non_bonded *stat_nb)
{
  int n_nonbonded;
  n_nonbonded = (n_particle_types*(n_particle_types+1))/2;

  obsstat_realloc_and_clear_non_bonded(stat_nb, n_nonbonded, 9);
  stat_nb->init_status_nb = 0;
}

/************************************************************/
void master_pressure_calc(int v_comp) {
  if(v_comp)
    mpi_gather_stats(3, total_pressure.data.e, total_p_tensor.data.e, total_pressure_non_bonded.data_nb.e, total_p_tensor_non_bonded.data_nb.e);
  else
    mpi_gather_stats(2, total_pressure.data.e, total_p_tensor.data.e, total_pressure_non_bonded.data_nb.e, total_p_tensor_non_bonded.data_nb.e);

  total_pressure.init_status = 1+v_comp;
  total_p_tensor.init_status = 1+v_comp;
  total_pressure_non_bonded.init_status_nb = 1+v_comp;
  total_p_tensor_non_bonded.init_status_nb = 1+v_comp;
}


/*****************************************************/
/* Routines for Local Stress Tensor                  */
/*****************************************************/

int getintersection(double pos1[3], double pos2[3],int given, int get, double value, double *answer, double box_size[3])
{
  /*pos1 and pos2 are two particle positions.                                                  */
  /*given and get are integers from 0 to 2. 0 = x direction. 1 = y direction. 2 = z direction  */
  /*there is a point on the line between the two particles p1 and p2 such that r[given]=value  */
  /*this procedure returns the value of r[get] at that point                                   */

  double p2r[3];
  double oldvalue;
  int i;

  for (i=0;i<3;i++) {
    p2r[i] = drem_down((pos2[i]-pos1[i])+box_size[i]/2.0,box_size[i])-box_size[i]/2.0;
  }
  oldvalue =value;
  value = drem_down((value-pos1[given])+box_size[given]/2.0,box_size[given])-box_size[given]/2.0;
  //PTENSOR_TRACE(fprintf(stderr,"%d: getintersection: p1 is %f %f %f p2 is %f %f %f p2r is %f %f %f value is %f newvalue is %f\n",this_node,pos1[0],pos1[1],pos1[2],pos2[0],pos2[1],pos2[2],p2r[0],p2r[1],p2r[2],oldvalue,value););
  
  if ((value)*(p2r[given]) < -0.0001) {
    char *errtxt = runtime_error(128 + 3*TCL_INTEGER_SPACE);
    ERROR_SPRINTF(errtxt, "{analyze stress_profile: getintersection: intersection is not between the two given particles - %e is not between %e and %e and box size is %e, given is %d\n ",value,0.0,p2r[given],box_size[given],given);
    return 0; 
  } else if (given == get) {
    *answer =  drem_down(value + pos1[given],box_size[given]);;
  } else if (0==p2r[given]) {
    char *errtxt = runtime_error(128 + 3*TCL_INTEGER_SPACE);
    ERROR_SPRINTF(errtxt, "{analyze stress_profile: getintersection: intersection is a line, not a point - value is %f same as %f and %f\n",value,0.0,p2r[given]);
    return 0;   
  } else {
    *answer =  drem_down(pos1[get]+p2r[get]/p2r[given]*value,box_size[get]);
  }
  return 1;
}

int getlength(double pos1[3], double pos2[3], int d1, double val1, int d2, double val2, int l, double *answer)
{
  /*p1 and p2 are two particles positions                                                          
    d1 and d2 are integers between 0 and 2 denoting an axis (x, y, or z)                        
    l and k are integers between 0 and 2 denoting an axis (x, y, or z)                          
    two points on the line connecting these particles are defined by r[d1]=val1 and r[d2]=val2  
    call these two points p3 and p4 (not program variables)                                     
    this program returns the distance between p3 and p4 in the l direction
  */
  
  double intersect1, intersect2;

  *answer = 0;
  
  if  (! getintersection(pos1,pos2,d2,l,val2,&intersect1,box_l) || ! getintersection(pos1,pos2,d1,l,val1,&intersect2,box_l)) {
    return 0;
  } else {
    *answer = drem_down(intersect2 - intersect1 + box_l[l]/2.0, box_l[l]) - box_l[l]/2.0;
    return 1;
  }
}

int does_line_go_through_cube(double pos1[3], double pos2[3], double range_start[3], double range[3], int sign[3], double entry[3], double exit[3],int *facein, int *faceout)
     /* p1 and p2 are two particle positions
        there is a cube in the simulation box with one vertex at range_start and the opposite vertex at range_start + range
        this routine calculates where the line connecting p1 and p2 enters and exits this cube
	these points are returned in the variables entry and exit
	the function returns a 0 if the line does not pass through the cube and 1 if it does
     */
{

  double centre[3];             /* centre of the cube */
  double vect2centre1[3];       /* vectors from centre of cube to pos1 and pos2 */
  double vect2centre2[3];
  int i;

  int doesntenter = 0;          /* boolean that indicates whether we have already determined that the line does not enter the cube at all */
  double intersection1 = 0.0, intersection2 = 0.0;
  int found_entry = 0;          /* boolean that indicates whether we have determined where the line enters */
  int found_exit = 0;           /* boolean that indicates whether we have determined where the line exits */
  int i1, i2;
  int inside1, inside2;

   /* find centre of analyzed cube */
  for (i=0;i<3;i++) {
    centre[i] = range_start[i] + range[i]/2.0;
  }

  /* find vectors connecting two particles to centre */
  get_mi_vector(vect2centre1,pos1,centre);
  get_mi_vector(vect2centre2,pos2,centre);

  *facein = -1;
  *faceout = -1;

  /* check if particles are inside cube */
  inside1 = 1;
  inside2 = 1;
  for (i=0;i<3;i++) {
    if (fabs(vect2centre1[i])>range[i]/2.0) inside1 = 0;
    if (fabs(vect2centre2[i])>range[i]/2.0) inside2 = 0;
  }
  PTENSOR_TRACE(fprintf(stderr,"%d: does_line_go_through_cube: Particle1 inside cube = %d Particle2 inside cube = %d \n",this_node,inside1,inside2););

  /* work out which face connecting line enters through */
  if (! inside1) {
    for (i=0;i<3;i++) {
      i1 = (i+1)%3;
      i2 = (i+2)%3;
      /*does the line start outside the cube in direction i?*/
      if ( (! doesntenter ) && (! found_entry) && fabs(vect2centre1[i])>range[i]/2.0 ) { 
	 /*does the bond heads away from the cube or is part2 is before the cube in direction i? */
	if ((vect2centre1[i] * sign[i] > 0) ||  (vect2centre2[i]*sign[i] < -range[i]/2.0)) {
	  doesntenter = 1;
	} else {
	  getintersection(vect2centre1, vect2centre2, i, i1, -range[i]/2.0*sign[i],&intersection1,box_l);
	  if (intersection1 > box_l[i1]/2) intersection1 -= box_l[i1];
	  if (fabs(intersection1) < range[i1]/2.0) {
	    getintersection(vect2centre1, vect2centre2, i, i2, -range[i]/2.0*sign[i],&intersection2,box_l);
	    if (intersection2 > box_l[i2]/2) intersection2 -= box_l[i2];
	    if (fabs(intersection2) < range[i2]/2.0) {
	      found_entry = 1;
	      entry[i] = centre[i] -range[i]/2.0*sign[i];
	      entry[i1] = centre[i1] + intersection1;
	      entry[i2] = centre[i2] + intersection2;
	      *facein = i; 
	    }
	  }
	}
      }
    }
  } else {
    entry[0] = pos1[0];
    entry[1] = pos1[1];
    entry[2] = pos1[2];
  }
  if (! (found_entry) && ! (inside1)) {
    PTENSOR_TRACE(fprintf(stderr,"%d: does_line_go_through_cube: Line does not pass through cube \n",this_node););
    return 0;
  } else {
    PTENSOR_TRACE(fprintf(stderr,"%d: does_line_go_through_cube: Entry is %f %f %f \n",this_node,entry[0], entry[1], entry[2]););
    if (! inside2) {
      /*check which outside faces of box the line exits through */
      for (i=0;i<3;i++) {
	i1 = (i+1)%3;
	i2 = (i+2)%3;
	/* does it enter into the box through the i face */
	if ((found_exit == 0) && (fabs(vect2centre2[i]) > range[i]/2.0)) {  /*starts outside cube in i direction */
	  getintersection(vect2centre1, vect2centre2, i, i1, range[i]/2.0*sign[i],&intersection1,box_l);
	  if (intersection1 > box_l[i1]/2) intersection1 -= box_l[i1];
	  if (fabs(intersection1) < range[i1]/2.0) {
	    getintersection(vect2centre1, vect2centre2, i, i2, range[i]/2.0*sign[i],&intersection2,box_l);
	    if (intersection2 > box_l[i2]/2) intersection2 -= box_l[i2];
	    if (fabs(intersection2) < range[i2]/2.0) {
	      found_exit = 1;
	      exit[i] = centre[i]   + range[i]/2.0*sign[i];
	      exit[i1] = centre[i1] + intersection1;
	      exit[i2] = centre[i2] + intersection2;
	      *faceout = i; 
	    }
	  }
	}
      }
    } else {
      exit[0] = pos2[0];
      exit[1] = pos2[1];
      exit[2] = pos2[2];
    }
    PTENSOR_TRACE(fprintf(stderr,"%d: does_line_go_through_cube: Exit is %f %f %f\n",this_node,exit[0], exit[1], exit[2]););
  }
    PTENSOR_TRACE(fprintf(stderr,"%d: does_line_go_through_cube: facein is %d faceout is %d\n",this_node,*facein,*faceout););
  return 1;
}

int distribute_tensors(DoubleList *TensorInBin, double *force, int bins[3], double range_start[3], double range[3], double pos1[3], double pos2[3])
{
  /*calculates how to distribute a force tensor between two particles between various bins
    the amount distributed to a bin is proportional to the length of the line between the
      two particles p1 and p2 that passes through the bin volume
    we consider a cube of space starting with a corner at {x_range_start, y_range_start, z_range_start} extending to 
      {x_range_start+x_range, y_range_start+y_range, z_range_start+z_range}
    this cube is split into x_bins bins in the x direction - we refer to these as x-bins
    each x-bin into y_bins bins in the y direction - we refer to these as y-bins
    each y_bin into z_bins bins in the z direction - we refer to these as z-bins
  */

  int sign[3],sign10[3];    /* indicates whether the line goes in the positive or negative direction in each dimension */
  double entry[3], exit[3]; /* the positions at which the line enters and exits the cube */
  int startx, endx;         /* x-bins in which the line starts and ends in */
  int occupiedxbins;        /* number of x-bins occuped by the line */
  int *starty;              /* y-bins in which the line starts in for each x-bin.  This array has dimension occupiedxbins+1. */
  int totoccupiedybins;     /* total number of y-bins through which the line passes  */
  int *occupiedybins;       /* number of occupied y-bins for each x-bin */
  int *occupiedzbins;       /* number of occupied z-bins for each y-bin */
  int *startz;              /* z-bins in which the line starts in for each y_bin.  This array has dimension totaloccupiedybins. */
  int xbin, ybin, zbin;     /* counters to keep track of bins x_bin goes from 0 to x_bins-1, y_bins from 0 to y_bins-1, z_bins from 0 to Z-bins-1 */
  int i ,k, l;    
  int counter;              /* keeps track of where we are in the startz array */
  int zi;                           
  double length;            /* length of line between the two points */
  int d1,d2;                /* for each z-bin d1 and d2 are calculated.  they indicate through which faces of the bin the line enters and leaves the bin.  i.e. if the line enters through the face corresponding to y = 0.34 then d1 = 1 (y-direction) and val1 = 0.34 */
  double val1, val2;
  double intersect;
  double segment,segment2;
  double calclength;        
  int xa, ya, za;           /* counters for bins */
  double temp[3];
  double redentry[3], redexit[3];  /* like entry and exit but using a coordinate system where range_start corresponds to (0,0,0) and the length scale in each direction is the bin width */
  double redbox_l[3];       /* box size is reduced units */
  int facein,faceout;

  PTENSOR_TRACE(fprintf(stderr,"%d: distribute_tensors: Distributing tensors for particle p1 %f %f %f and p2 %f %f %f\n",this_node,pos1[0],pos1[1],pos1[2],pos2[0],pos2[1],pos2[2]););
  PTENSOR_TRACE(fprintf(stderr,"%d: distribute_tensors: range_start is %f %f %f range is %f %f %f bins is %d %d %d\n",this_node,range_start[0],range_start[1],range_start[2],range[0],range[1],range[2],bins[0],bins[1],bins[2]););
  /* work out what direction the line joining the particles goes in */
  length = 0;
  get_mi_vector(temp, pos2, pos1);
  for (i=0;i<3;i++) {
    sign[i] = (temp[i] > 0)*2-1;
    sign10[i]  = (sign[i]+1)/2.0;
  }
  PTENSOR_TRACE(fprintf(stderr,"%d: distribute_tensors: sign is %d %d %d\n",this_node,sign[0],sign[1],sign[2]););
  
  calclength = 0;

  /* check if the line passes through the cube and if so where it enters and exits it */
  if (does_line_go_through_cube(pos1, pos2, range_start, range, sign, entry, exit,&facein,&faceout)) {
    
    /* calculate reduced coordinates */
    /* range_start becomes (0,0,0) and range_start + range becomes (bins[0],bins[1],bins[2]) */
    get_mi_vector(redentry, entry, range_start);
    get_mi_vector(redexit, exit, range_start);
    for (i=0;i<3;i++) {
      redentry[i] = drem_down(redentry[i]+box_l[i],box_l[i]);
      redexit[i] = drem_down(redexit[i]+box_l[i],box_l[i]);
    }
    length = 0;
    get_mi_vector(temp, exit, entry);
    for (i=0;i<3;i++) {
      redentry[i] = redentry[i]/range[i]*bins[i];
      redexit[i] = redexit[i]/range[i]*bins[i];
      length += pow(temp[i],2);
    }
    length = pow(length,0.5);
    PTENSOR_TRACE(fprintf(stderr,"%d: distribute_tensors: Reduced Entry is %f %f %f Reduced Exit is %f %f %f \n",this_node,redentry[0],redentry[1],redentry[2],redexit[0], redexit[1], redexit[2]););

    for (i=0;i<3;i++) {
      redbox_l[i] = box_l[i]/range[i]*bins[i];
    }
   
    /* find in which x-bins the line starts and stops */ 
    if (facein == 0) {
      startx = dround(redentry[0]) - 1 + sign10[0];
    } else {
      startx = floor(redentry[0]);
    }
    if ((startx < 0) && (range[0]==box_l[0])) startx += bins[0];
    if (faceout == 0) {
      endx = dround(redexit[0] - sign10[0]);
    } else {
      endx = floor(redexit[0]);
    }
    if ((endx < 0) && (range[0]==box_l[0])) startx += bins[0];
    occupiedxbins = (sign[0]*(endx-startx) + bins[0])%bins[0] +1;
  
    PTENSOR_TRACE(fprintf(stderr,"%d: distribute_tensors: x goes from %d to %d\n",this_node,startx, endx);)
    /* Initialise starty array */
    starty = (int *)malloc(sizeof(int)*(occupiedxbins+1));
    occupiedybins = (int *)malloc(sizeof(int)*occupiedxbins);

    /* find in which y-bins the line starts and stops for each x-bin */
    /* in xbin the line starts in y-bin number starty[xbin-startx] and ends in starty[xbin-startx+1] */
    totoccupiedybins = 0;
    if (facein == 1) {
      starty[0] = dround(redentry[1]) - 1 + sign10[1];
    } else {
      starty[0] = floor(redentry[1]);
    }
    if ((starty[0] < 0) && (range[1]==box_l[1])) starty[0] += bins[1];
    for (xa=0; xa < occupiedxbins; xa++) {
      xbin = (startx + xa*sign[0]+bins[0])%bins[0];
      if (xbin == endx) {
	intersect  = redexit[1];
	if (faceout == 1) intersect -= sign10[1];
	if (( intersect < 0) && (range[1]==box_l[1])) intersect += bins[1];
      } else {
	if (getintersection(redentry,redexit,0,1,xbin+sign10[0],&intersect,redbox_l) != 1) return 0;
      }
      starty[xa+1]=floor(intersect);
      occupiedybins[xa] = ((starty[xa+1]-starty[xa])*sign[1]+bins[1])%bins[1]+1;
      totoccupiedybins += occupiedybins[xa];
      PTENSOR_TRACE(fprintf(stderr,"%d: distribute_tensors: in xbin %d y goes from %d to %d\n",this_node, xbin, starty[xa],starty[xa+1]););
    }

    /* Initialise startz array */
    occupiedzbins = (int *)malloc(sizeof(int)*totoccupiedybins);
    startz = (int *)malloc(sizeof(int)*(totoccupiedybins+1));
    /* find in which z-bins the line starts and stops for each y-bin*/
    counter = 0;
    if (facein == 2) {
      zi = dround(redentry[2]) - 1 + sign10[2];
    } else {
      zi = floor(redentry[2]);
    }
    if (( zi < 0) && (range[2]==box_l[2])) zi += bins[2];
    startz[counter] = zi;

    for (xa=0; xa < occupiedxbins; xa++) {
      xbin = (startx + xa*sign[0]+bins[0])%bins[0];
      for (ya = 0; ya < occupiedybins[xa]-1; ya++) {
	ybin = (starty[xa] + ya*sign[1] + bins[1])%bins[1];
	PTENSOR_TRACE(fprintf(stderr,"%d: distribute_tensors: xbin is %d ya is %d, occupiedybins[xa] is %d, ybin is %d\n",this_node,xbin,ya,occupiedybins[xa],ybin););
	if (getintersection(redentry,redexit,1,2,ybin+sign10[1],&intersect,redbox_l) != 1) return 0;
	zi = floor(intersect);
	startz[counter+1] = zi;
	occupiedzbins[counter]= (sign[2]*(startz[counter+1]-startz[counter]) + bins[2])%bins[2] +1;
	PTENSOR_TRACE(fprintf(stderr,"%d: distribute_tensors: in xbin %d ybin %d zbin goes from %d to %d\n",this_node,xbin,ybin,startz[counter],startz[counter+1]););
	counter ++;
      }
      ybin = starty[xa+1];
      if (xbin == endx) {
	if (faceout == 2) {
	  zi = dround(redexit[2] -sign10[2]);
	} else {
	  zi = floor(redexit[2]);
	}
	if (( zi < 0) && (range[2]==box_l[2])) zi += bins[2];
      
      } else {
	if (getintersection(redentry,redexit,0,2,xbin+sign10[0], &intersect, redbox_l) != 1) return 0;
	zi = floor(intersect);
      }
      startz[counter+1]=zi;
      occupiedzbins[counter]= (sign[2]*(startz[counter+1]-startz[counter]) + bins[2])%bins[2] +1;
      PTENSOR_TRACE(fprintf(stderr,"%d: distribute_tensors: in xbin %d ybin %d zbin goes from %d to %d\n",this_node,xbin,ybin,startz[counter],startz[counter+1]););
      counter ++;
    }
    PTENSOR_TRACE(fprintf(stderr,"%d: distribute_tensors: occupiedxbins is %d,occupiedybins[0] is %d, occupiedzbins[0] is %d\n",this_node,occupiedxbins,occupiedybins[0],occupiedzbins[0]););

    /* find out what length of the line passes through each z-bin */
    counter = 0;
    for (xa = 0; xa < occupiedxbins; xa ++) {
      xbin = (startx + xa*sign[0]+bins[0])%bins[0];
      for (ya = 0; ya < occupiedybins[xa]; ya++) {
	ybin = (starty[xa] + ya*sign[1] + bins[1])%bins[1];
	for (za = 0; za < occupiedzbins[counter]; za++) {
	  zbin = (startz[counter] + za*sign[2] +bins[2])%bins[2];
	  if (zbin == startz[counter]) {
	    if (ybin == starty[xa]) {
	      if (xbin == startx) {
		if (entry[0]-exit[0] != 0) {
		  d1 = 0;
		  val1 = redentry[0];
		} else if (entry[1]-exit[1] != 0){
		  d1 = 1;
		  val1 = redentry[1];
		} else {
		  d1 = 2;
		  val1 = redentry[2];
		}                                                        
		PTENSOR_TRACE(fprintf(stderr,"%d: distribute_tensors: %d %d %d line starts at point p1 inside bin\n",this_node,xbin,ybin,zbin););
	      } else {
		d1 = 0;
		val1  = xbin+1-sign10[0];
		PTENSOR_TRACE(fprintf(stderr,"%d: distribute_tensors: %d %d %d line enters through x = %e face of box\n",this_node,xbin,ybin,zbin,val1););
	      }  
	    } else {
	      d1 = 1;
	      val1 = ybin+1-sign10[1];        
	      PTENSOR_TRACE(fprintf(stderr,"%d: distribute_tensors: %d %d %d line enters through y = %e face of box\n",this_node,xbin,ybin,zbin,val1););
	    }
	  } else {
	    d1 = 2;
	    val1 = zbin+1-sign10[2]; 
	    PTENSOR_TRACE(fprintf(stderr,"%d: distribute_tensors: %d %d %d line enters through z = %e face of box zbin is %d\n",this_node,xbin,ybin,zbin,val1,zbin););
	  }
	  if (zbin == startz[counter+1]) {
	    if (ybin == starty[xa+1]) {
	      if (xbin == endx) {
		if (entry[0]-exit[0] != 0) {
		  d2 = 0;
		  val2 = redexit[0];
		} else if (entry[1]-exit[1]!= 0){
		  d2 = 1;
		  val2 = redexit[1];
		} else {
		  d2 = 2;
		  val2 = redexit[2];
		}                                                      
		PTENSOR_TRACE(fprintf(stderr,"%d: distribute_tensors: %d %d %d line ends at p2 inside bin\n",this_node,xbin,ybin,zbin););
	      }
	      else {
		d2 = 0;
		val2 = xbin+sign10[0];
		PTENSOR_TRACE(fprintf(stderr,"%d: distribute_tensors: %d %d %d line leaves through x = %e face of box\n",this_node,xbin,ybin,zbin,val2););
	      }
	    }
	    else {
	      d2 = 1;
	      val2 = ybin+sign10[1]; 
	      PTENSOR_TRACE(fprintf(stderr,"%d: distribute_tensors: %d %d %d line leaves through y = %e face of box direction %d\n",this_node,xbin,ybin,zbin,val2,sign[1]););
	    }
	  }
        else {
          d2= 2;
          val2 = zbin+sign10[2];
	  PTENSOR_TRACE(fprintf(stderr,"%d: distribute_tensors: %d %d %d line leaves through z = %e face of box\n",this_node,xbin,ybin,zbin,val2););
        }
	  segment2 = 0;
	    PTENSOR_TRACE(fprintf(stderr,"%d: distribute_tensors: entry %f %f %f exit %f %f %f d1 %d val1 %f d2 %d val2 %f\n",this_node,entry[0],entry[1],entry[2],exit[0],exit[1],exit[2],d1,range_start[d1]+val1*range[d1]/(double)bins[d1],d2,range_start[d2]+val2*range[d2]/(double)bins[d2]););
	  for (l=0;l<3;l++) {
	    if (getlength(entry,exit,d1,range_start[d1]+val1*range[d1]/(double)bins[d1],d2,range_start[d2]+val2*range[d2]/(double)bins[d2],l,&segment) != 1) return 0;
	    segment2 += pow(segment,2);
	    for (k=0;k<3;k++) {
	      TensorInBin[xbin*bins[1]*bins[2]+ybin*bins[2]+zbin].e[3*k+l] += force[k] * segment;
	    }
	  }
	  calclength += pow(segment2,0.5);
	}
	counter ++;
      }
    }
    PTENSOR_TRACE(fprintf(stderr,"%d: distribute_tensors: calclength is %e and length is %e\n}",this_node,calclength,length););
    
    if (calclength - length >0.0000000001) {
      char *errtxt = runtime_error(128 + 3*TCL_INTEGER_SPACE);
      ERROR_SPRINTF(errtxt, "{%d: analyze stess_profile: bug in distribute tensor code - calclength is %e and length is %e}",this_node,calclength,length);
      return 0;
    }
    free(occupiedzbins);
    free(occupiedybins);
    free(starty);
    free(startz);
  }
  return 1;
} 

int reducepos(double pos[3], int bins[3], double centre[3], double range[3], int reducedpos[3]) {
  int i;
  double working[3];
  get_mi_vector(working, pos, centre);
  for (i=0;i<3;i++) {
    reducedpos[i] = floor((working[i]+range[i]/2.0)*(double)bins[i]/range[i]);
  }
  return 1; 
}

int incubewithskin(double pos[3], double centre[3], double range[3])
{
  double working[3];
  int i;
  get_mi_vector(working, pos, centre);
  for (i=0; i <3; i++) {
    if (fabs(working[i]) > range[i]/2.0 + skin+max_cut) return 0;
  }
  return 1; 
}

int whichbin(double pos[3], int bins[3], double centre[3], double range[3], int *bin)
/*calculates which bin a particle is in for local_stress_tensor */
{
  int reducedpos[3];
  int i;
  reducepos(pos, bins, centre, range, reducedpos);
  for (i=0;i<3;i++) {
    if ((reducedpos[i] < 0) || (reducedpos[i] >= bins[i])) {
      *bin = -1;
      return 1;
    }
  }
  *bin = reducedpos[0]*bins[1]*bins[2] + reducedpos[1]*bins[2] + reducedpos[2];
  return 1;
}

int get_nonbonded_interaction(Particle *p1, Particle *p2, double *force)
{
  /* returns the non_bonded interaction between two particles */

  double dist2, dist;
  double d[3];
#ifdef ELECTROSTATICS
  int i;
  double eforce[3];
#endif

  force[0]=0; force[1]=0; force[2]=0; 
  

  if ((p1->p.identity != p2->p.identity)&&(checkIfParticlesInteract(p1->p.type, p2->p.type))) {
    /* distance calculation */
    get_mi_vector(d, p1->r.p, p2->r.p);
    dist2 = SQR(d[0]) + SQR(d[1]) + SQR(d[2]);
    dist  = sqrt(dist2);
    calc_non_bonded_pair_force_simple(p1,p2,d,dist,dist2,force);
#ifdef ELECTROSTATICS
    if (coulomb.method != COULOMB_NONE) {
      switch (coulomb.method) {
#ifdef ELP3M
      case COULOMB_P3M:
	fprintf(stderr,"WARNING: Local stress tensor calculation cannot handle P3M electrostatics so it is left out\n");  
	break;
#endif
      case COULOMB_DH:
	for (i = 0; i < 3; i++)
	  eforce[i] = 0;
	add_dh_coulomb_pair_force(p1,p2,d,dist, eforce);
	for(i=0;i<3;i++)
	  force[i] += eforce[i];
	break;
      case COULOMB_RF:
	for (i = 0; i < 3; i++)
	  eforce[i] = 0;
	add_rf_coulomb_pair_force(p1,p2,d,dist, eforce);
	for(i=0;i<3;i++)
	    force[i] += eforce[i];
	break;
      case COULOMB_INTER_RF:
        // this is done elsewhere
	break;
      case COULOMB_MMM1D:
	fprintf(stderr,"WARNING: Local stress tensor calculation cannot handle MMM1D electrostatics so it is left out\n");  
      default:
	fprintf(stderr,"WARNING: Local stress tensor calculation does not recognise this electrostatic interaction\n");  
      }
    }
#endif /*ifdef ELECTROSTATICS */

#ifdef MAGNETOSTATICS
    if (coulomb.Dmethod != DIPOLAR_NONE) {
      switch (coulomb.Dmethod) {
#ifdef ELP3M
      case DIPOLAR_P3M:
    	fprintf(stderr,"WARNING: Local stress tensor calculation cannot handle P3M magnetostatics so it is left out\n");  
	break;
#endif
#ifdef DAWAANR
      case DIPOLAR_ALL_WITH_ALL_AND_NO_REPLICA:
    	fprintf(stderr,"WARNING: Local stress tensor calculation cannot handle DAWAANR magnetostatics so it is left out\n");  
	break;
#endif
#ifdef MAGNETIC_DIPOLAR_DIRECT_SUM
      case DIPOLAR_DS:
    	fprintf(stderr,"WARNING: Local stress tensor calculation cannot handle MAGNETIC DIPOLAR SUM magnetostatics so it is left out\n");  
	break;
#endif

      default:
	fprintf(stderr,"WARNING: Local stress tensor calculation does not recognise this magnetostatic interaction\n");  
      }
    }
#endif /*ifdef MAGNETOSTATICS */


  } /*if p1-> ... */
  return 0;
}

int local_stress_tensor_calc(DoubleList *TensorInBin, int bins[3], int periodic[3], double range_start[3], double range[3])
{
  /*calculates local stress tensors in cuboid bins
    uses Irving Kirkwood method
    we consider a cube of space starting with a corner at position range_start extending to 
      range_start + range
    if the variable periodic is set to 1 in dimension i then the cube is assumed to span the periodic box
    this cube is divided into bins[0] bins in the x direction bins[1] in the y direction, and bins[2] in the z direction
  */

  int i,j;                       /*counter for dimension */
  double binvolume;
  int c, np, n, bin;
  double centre[3];
  Cell *cell;
  Particle *p1, *p2, **pairs;
  Particle *particles;
  double force[3];
  int k,l;
  int type_num;
  Bonded_ia_parameters *iaparams;
  double dx[3];

  
  for (i = 0; i < 3; i ++) {
    if (periodic[i]) {
      range[i] = box_l[i];
      range_start[i] = 0;
    }
  }

  /* find centre of analyzed cube */
  for (i=0;i<3;i++) {
    centre[i] = range_start[i] + range[i]/2.0;
  }

  /* We consider all particles that are within a certain distance of the cube. The skin is used as this distance.  If the
     skin from on opposite sides of the box overlaps then we produce an error message.  To code dround this would be
     creating unnecessary work since I can't imagine when we might want that */

  if (skin < 0) {
    char *errtxt = runtime_error(128 + 3*TCL_INTEGER_SPACE);
    ERROR_SPRINTF(errtxt, "{analyze stess_profile: skin cannot be negative}");
    return 0;
  }
  for (i=0;i<3;i++) {
    if ((! periodic[i]) && (range[i] + 2*skin +2*max_cut > box_l[i])) {
      char *errtxt = runtime_error(128 + 3*TCL_INTEGER_SPACE);
      ERROR_SPRINTF(errtxt, "{analyze stress_profile: Analyzed box (%f) with skin+max_cut(%f) is larger than simulation box (%f).\n",range[i],skin+max_cut,box_l[i]);
      return 0;
    }
    range_start[i] = drem_down(range_start[i],box_l[i]);
  }
  PTENSOR_TRACE(fprintf(stderr,"%d: Running stress_profile\n",this_node));

  binvolume = range[0]*range[1]*range[2]/(double)bins[0]/(double)bins[1]/(double)bins[2];

  /* this next bit loops over all pair of particles, calculates the force between them, and distributes it amongst the tensors */

  // loop over all local cells
  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    particles   = cell->part;
    np  = cell->n;
    // loop over all particles in this cell
    for(i = 0; i < np; i++)  {
      p1 = &(particles[i]);
      whichbin(p1->r.p,bins,centre, range, &bin); 
      if (bin >= 0) {
	PTENSOR_TRACE(fprintf(stderr,"%d:Got particle number %d i is %d pos is %f %f %f \n",this_node,p1->p.identity,i,p1->r.p[0],p1->r.p[1],p1->r.p[2]));
	PTENSOR_TRACE(fprintf(stderr,"%d:Ideal gas component is {",this_node));
	for(k=0;k<3;k++) {
	  for(l=0;l<3;l++) {
	    TensorInBin[bin].e[k*3 + l] += (p1->m.v[k])*(p1->m.v[l])*PMASS(*p1)/time_step/time_step;
	    PTENSOR_TRACE(fprintf(stderr,"%f ",(p1->m.v[k])*(p1->m.v[l])*PMASS(*p1)/time_step/time_step));
	  }
	}

	PTENSOR_TRACE(fprintf(stderr,"}\n"));
      }
      
      /* bonded contributions */
      j = 0;
      while(j < p1->bl.n) {
	type_num = p1->bl.e[j++];
	iaparams = &bonded_ia_params[type_num];

	/* fetch particle 2 */
	p2 = local_particles[p1->bl.e[j++]];
	get_mi_vector(dx, p1->r.p, p2->r.p);
	calc_bonded_force(p1,p2,iaparams,&j,dx,force);
	PTENSOR_TRACE(fprintf(stderr,"%d: Bonded to particle %d with force %f %f %f\n",this_node,p2->p.identity,force[0],force[1],force[2]));
	if ((pow(force[0],2)+pow(force[1],2)+pow(force[2],2)) > 0) {
	  if (distribute_tensors(TensorInBin,force,bins,range_start,range,p1->r.p, p2->r.p) != 1) return 0;
	}
      }
    }

    // Loop cell neighbors
    for (n = 0; n < dd.cell_inter[c].n_neighbors; n++) {
      pairs = dd.cell_inter[c].nList[n].vList.pair;
      np    = dd.cell_inter[c].nList[n].vList.n;

      // verlet list loop //
      for(i=0; i<2*np; i+=2) {
	p1 = pairs[i];                    // pointer to particle 1
	p2 = pairs[i+1];                  // pointer to particle 2
	if ((incubewithskin(p1->r.p,centre,range)) && (incubewithskin(p2->r.p,centre,range))) {
	  get_nonbonded_interaction(p1,p2, force);
	  PTENSOR_TRACE(fprintf(stderr,"%d:Looking at pair %d %d force is %f %f %f\n",this_node,p1->p.identity, p2->p.identity,force[0],force[1], force[2]));
	  if ((pow(force[0],2)+pow(force[1],2)+pow(force[2],2)) > 0) {
	    if (distribute_tensors(TensorInBin,force,bins,range_start,range,p1->r.p, p2->r.p) != 1) return 0;
	  }
	} else {
	  // PTENSOR_TRACE(fprintf(stderr,"%d:Looking at pair %d %d not in cube with skin\n",this_node,p1->p.identity, p2->p.identity));
	}
      }
    }
  }

  for (i=0;i<bins[0]*bins[1]*bins[2];i++) {
    for (j=0;j<9;j++) {
	TensorInBin[i].e[j] /= binvolume;
    }
  }

  return 1;
}


/****************************************************************************************
 *                                 parser
 ****************************************************************************************/

static void print_detailed_pressure(Tcl_Interp *interp)
{
  char buffer[TCL_DOUBLE_SPACE + TCL_INTEGER_SPACE + 2];
  double value;
  int i, j;
  value = 0.0;

  value = total_pressure.data.e[0];
  for (i = 1; i < total_pressure.data.n; i++)
    value += total_pressure.data.e[i];

  Tcl_PrintDouble(interp, value, buffer);
  Tcl_AppendResult(interp, "{ pressure ", buffer, " } ", (char *)NULL);

  Tcl_PrintDouble(interp, total_pressure.data.e[0], buffer);
  Tcl_AppendResult(interp, "{ ideal ", buffer, " } ", (char *)NULL);

  for(i=0;i<n_bonded_ia;i++) {
    if (bonded_ia_params[i].type != BONDED_IA_NONE) {
      sprintf(buffer, "%d ", i);
      Tcl_AppendResult(interp, "{ ", buffer, (char *)NULL);
      Tcl_PrintDouble(interp, *obsstat_bonded(&total_pressure, i), buffer);
      Tcl_AppendResult(interp,
		       get_name_of_bonded_ia(bonded_ia_params[i].type),
		       " ", buffer, " } ", (char *) NULL);
    }
  }

  for (i = 0; i < n_particle_types; i++)
    for (j = i; j < n_particle_types; j++) {
      if (checkIfParticlesInteract(i, j)) {
	sprintf(buffer, "%d ", i);
	Tcl_AppendResult(interp, "{ ", buffer, (char *)NULL);
	sprintf(buffer, "%d ", j);
	Tcl_AppendResult(interp, " ", buffer, (char *)NULL);
	Tcl_PrintDouble(interp, *obsstat_nonbonded(&total_pressure, i, j), buffer);
	Tcl_AppendResult(interp, "nonbonded ", buffer, " } ", (char *)NULL);
      }
    }
  
/* In case we need intra- and inter- nonbonded (nb) contribution of total pressure  */
  value = 0.0;
  for (i = 0; i < n_particle_types; i++)
    for (j = i; j < n_particle_types; j++) {
      value += *obsstat_nonbonded_intra(&total_pressure_non_bonded, i, j);
    }
  Tcl_PrintDouble(interp, value, buffer);
  Tcl_AppendResult(interp, " { total_nb_intra ", (char *)NULL);
  Tcl_PrintDouble(interp, value, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  Tcl_AppendResult(interp, "} ", (char *)NULL);
 
  value = 0.0;
  for (i = 0; i < n_particle_types; i++)
    for (j = i; j < n_particle_types; j++) {
      value += *obsstat_nonbonded_inter(&total_pressure_non_bonded, i, j);
    }
  Tcl_PrintDouble(interp, value, buffer);
  Tcl_AppendResult(interp, " { total_nb_inter ", (char *)NULL);
  Tcl_PrintDouble(interp, value, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  Tcl_AppendResult(interp, "} ", (char *)NULL);
  
  for (i = 0; i < n_particle_types; i++)
    for (j = i; j < n_particle_types; j++) {
      if (checkIfParticlesInteract(i, j)) {
	sprintf(buffer, "%d ", i);
	Tcl_AppendResult(interp, "{ ", buffer, (char *)NULL);
	sprintf(buffer, "%d ", j);
	Tcl_AppendResult(interp, " ", buffer, (char *)NULL);
	Tcl_PrintDouble(interp, *obsstat_nonbonded_intra(&total_pressure_non_bonded, i, j), buffer);
	Tcl_AppendResult(interp, "nb_intra ", buffer, " } ", (char *)NULL);
      }
    }
  
  for (i = 0; i < n_particle_types; i++)
    for (j = i; j < n_particle_types; j++) {
      if (checkIfParticlesInteract(i, j)) {
	sprintf(buffer, "%d ", i);
	Tcl_AppendResult(interp, "{ ", buffer, (char *)NULL);
	sprintf(buffer, "%d ", j);
	Tcl_AppendResult(interp, " ", buffer, (char *)NULL);
	Tcl_PrintDouble(interp, *obsstat_nonbonded_inter(&total_pressure_non_bonded, i, j), buffer);
	Tcl_AppendResult(interp, "nb_inter ", buffer, " } ", (char *)NULL);
      }
    } 
  
#if  defined(ELECTROSTATICS) || defined (MAGNETOSTATICS)
  if( 
#ifdef ELECTROSTATICS
      coulomb.method != COULOMB_NONE
#else
      0
#endif
      ||
#ifdef MAGNETOSTATICS
      coulomb.Dmethod != DIPOLAR_NONE
#else
      0
#endif     
   ) {
    /* total Coulomb pressure */
    value = total_pressure.coulomb[0];
    for (i = 1; i < total_pressure.n_coulomb; i++)
      value += total_pressure.coulomb[i];
    for (i = 0; i < total_pressure.n_dipolar; i++)
      value += total_pressure.dipolar[i];
    Tcl_PrintDouble(interp, value, buffer);
#if  defined(ELECTROSTATICS) && defined (MAGNETOSTATICS)
    Tcl_AppendResult(interp, "{ coulomb+magdipoles ", buffer, (char *)NULL);
#else
#if defined(ELECTROSTATICS)
    Tcl_AppendResult(interp, "{ coulomb ", buffer, (char *)NULL);
#endif	
#if defined(MAGNETOSTATICS)
    Tcl_AppendResult(interp, "{ magdipoles ", buffer, (char *)NULL);
#endif		
#endif

    /* if it is split up, then print the split up parts */
    if (total_pressure.n_coulomb > 1) {
      for (i = 0; i < total_pressure.n_coulomb; i++) {
	Tcl_PrintDouble(interp, total_pressure.coulomb[i], buffer);
	Tcl_AppendResult(interp, " ", buffer, (char *)NULL);
      }
    }
    Tcl_AppendResult(interp, " }", (char *)NULL);
  }
#endif
}

/************************************************************/

int parse_and_print_pressure(Tcl_Interp *interp, int v_comp, int argc, char **argv)
{
  /* 'analyze pressure [{ bond <type_num> | nonbonded <type1> <type2> | coulomb | ideal | total }]' */
  char buffer[TCL_DOUBLE_SPACE + TCL_INTEGER_SPACE + 2];
  int i, j;
  double value, p_vel[3];
  value = 0.0;

  if (n_total_particles == 0) {
    Tcl_AppendResult(interp, "(no particles)",
		     (char *)NULL);
    return (TCL_OK);
  }

  /* if desired (v_comp==1) replace ideal component with instantaneous one */
  if (total_pressure.init_status != 1+v_comp ) {
    init_virials(&total_pressure);
    init_p_tensor(&total_p_tensor);
    
    init_virials_non_bonded(&total_pressure_non_bonded);
    init_p_tensor_non_bonded(&total_p_tensor_non_bonded);
    
    if(v_comp && (integ_switch == INTEG_METHOD_NPT_ISO) && !(nptiso.invalidate_p_vel)) {
      if (total_pressure.init_status == 0)
	master_pressure_calc(0);
      total_pressure.data.e[0] = 0.0;
      MPI_Reduce(nptiso.p_vel, p_vel, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      for(i=0; i<3; i++)
	if(nptiso.geometry & nptiso.nptgeom_dir[i])
	  total_pressure.data.e[0] += p_vel[i];
      total_pressure.data.e[0] /= (nptiso.dimension*nptiso.volume);
      total_pressure.init_status = 1+v_comp;   }
    else
      master_pressure_calc(v_comp);
  }

  if (argc == 0)
    print_detailed_pressure(interp);
  else {

    if      (ARG0_IS_S("ideal"))
      value = total_pressure.data.e[0];
    else if (ARG0_IS_S("bonded") ||
	     ARG0_IS_S("fene") ||
	     ARG0_IS_S("subt_lj_harm") ||
	     ARG0_IS_S("subt_lj_fene") ||
	     ARG0_IS_S("subt_lj") ||
	     ARG0_IS_S("harmonic")) {
      if(argc<2 || ! ARG1_IS_I(i)) {
	Tcl_ResetResult(interp);
	Tcl_AppendResult(interp, "wrong # or type of arguments for: analyze pressure bonded <type_num>",
			 (char *)NULL);
	return (TCL_ERROR);
      }
      if(i < 0 || i >= n_bonded_ia) {
	Tcl_AppendResult(interp, "bond type does not exist", (char *)NULL);
	return (TCL_ERROR);
      }
      value = *obsstat_bonded(&total_pressure, i);
    }
    else if (ARG0_IS_S("nonbonded") ||
	     ARG0_IS_S("lj") ||
	     ARG0_IS_S("buckingham") ||
             ARG0_IS_S("morse") ||
             ARG0_IS_S("soft-sphere") ||
	     ARG0_IS_S("lj-cos") ||
	     ARG0_IS_S("lj-cos2") ||
	     ARG0_IS_S("tabulated") ||
	     ARG0_IS_S("gb")) {
      if(argc<3 || ! ARG_IS_I(1, i) || ! ARG_IS_I(2, j)) {
	Tcl_ResetResult(interp);
	Tcl_AppendResult(interp, "wrong # or type of arguments for: analyze pressure nonbonded <type1> <type2>",
			 (char *)NULL);
	return (TCL_ERROR);
      }
      if(i < 0 || i >= n_particle_types || j < 0 || j >= n_particle_types) {
	Tcl_AppendResult(interp, "particle type does not exist", (char *)NULL);
	return (TCL_ERROR);
      }
      value = *obsstat_nonbonded(&total_pressure, i, j);
    }
    else if( ARG0_IS_S("tot_nb_intra") ||
	     ARG0_IS_S("tot_nonbonded_intra")) {
      value = 0.0;
      for (i = 0; i < n_particle_types; i++)
        for (j = i; j < n_particle_types; j++)
        value += *obsstat_nonbonded_intra(&total_pressure_non_bonded, i, j);
    }
    else if( ARG0_IS_S("tot_nb_inter") ||
	     ARG0_IS_S("tot_nonbonded_inter")) {
      value = 0.0;
      for (i = 0; i < n_particle_types; i++)
        for (j = i; j < n_particle_types; j++)
        value += *obsstat_nonbonded_inter(&total_pressure_non_bonded, i, j);
    }
    else if( ARG0_IS_S("nb_intra") ||
	     ARG0_IS_S("nonbonded_intra")) {
      if(argc<3 || ! ARG_IS_I(1, i) || ! ARG_IS_I(2, j)) {
	Tcl_ResetResult(interp);
	Tcl_AppendResult(interp, "wrong # or type of arguments for: analyze pressure nb_intra <type1> <type2>",
			 (char *)NULL);
	return (TCL_ERROR);
      }
      if(i < 0 || i >= n_particle_types || j < 0 || j >= n_particle_types) {
	Tcl_AppendResult(interp, "particle type does not exist", (char *)NULL);
	return (TCL_ERROR);
      }
      value = *obsstat_nonbonded_intra(&total_pressure_non_bonded, i, j);
    }   
    else if( ARG0_IS_S("nb_inter") ||
	     ARG0_IS_S("nonbonded_inter")) {
      if(argc<3 || ! ARG_IS_I(1, i) || ! ARG_IS_I(2, j)) {
	Tcl_ResetResult(interp);
	Tcl_AppendResult(interp, "wrong # or type of arguments for: analyze pressure nb_inter <type1> <type2>",
			 (char *)NULL);
	return (TCL_ERROR);
      }
      if(i < 0 || i >= n_particle_types || j < 0 || j >= n_particle_types) {
	Tcl_AppendResult(interp, "particle type does not exist", (char *)NULL);
	return (TCL_ERROR);
      }
      value = *obsstat_nonbonded_inter(&total_pressure_non_bonded, i, j);
    }
    else if( ARG0_IS_S("coulomb")) {
#ifdef ELECTROSTATICS
      value = 0;
      for (i = 0; i < total_pressure.n_coulomb; i++)
	value += total_pressure.coulomb[i];
#else
      Tcl_AppendResult(interp, "ELECTROSTATICS not compiled (see config.h)\n", (char *)NULL);
#endif
    }
    else if( ARG0_IS_S("dipolar")) {
#ifdef MAGNETOSTATICS
      value = 0;
      for (i = total_pressure.n_coulomb-1; i < total_pressure.n_coulomb; i++)  /*when DLC will be installed this has to be changed */
        value += total_pressure.coulomb[i];
#else
      Tcl_AppendResult(interp, "MAGNETOSTATICS not compiled (see config.h)\n", (char *)NULL);
#endif
    }
    else if (ARG0_IS_S("total")) {
      value = total_pressure.data.e[0];
      for (i = 1; i < total_pressure.data.n; i++) {
	value += total_pressure.data.e[i];
      }
    }
    else {
      Tcl_AppendResult(interp, "unknown feature of: analyze pressure",
		       (char *)NULL);
      return (TCL_ERROR);
    }
    Tcl_PrintDouble(interp, value, buffer);
    Tcl_AppendResult(interp, buffer, (char *)NULL);
  }

  return (TCL_OK);
}

int parse_and_print_p_IK1(Tcl_Interp *interp, int argc, char **argv)
{
  Tcl_AppendResult(interp, "analyze p_IK no longer exists", (char *)NULL);
  return (TCL_ERROR);
}
 
int parse_bins(Tcl_Interp *interp, int argc, char **argv)
{
  Tcl_AppendResult(interp, "analyze bins no longer exists", (char *)NULL);
  return (TCL_ERROR);
}

/****************************************************************************************
 *                                 parser
 ****************************************************************************************/
/* Parser routines used in the "analyze stress_tensor" command */ 
static void print_detailed_stress_tensor(Tcl_Interp *interp)
{
  char buffer[TCL_DOUBLE_SPACE + TCL_INTEGER_SPACE + 2];
  double value, tvalue[9];
  int i, j, k;
  value = 0.0;

  Tcl_AppendResult(interp, "{ pressure ", (char *)NULL);
  for(j=0; j<9; j++) {
    value = total_p_tensor.data.e[j];
    for (i = 1; i < total_p_tensor.data.n/9; i++) value += total_p_tensor.data.e[9*i + j];
    Tcl_PrintDouble(interp, value, buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  }
  Tcl_AppendResult(interp, "} ", (char *)NULL);

  Tcl_AppendResult(interp, "{ ideal ", (char *)NULL);
  for(j=0; j<9; j++) {
    Tcl_PrintDouble(interp, total_p_tensor.data.e[j], buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  }
  Tcl_AppendResult(interp, "} ", (char *)NULL);

  for(i=0;i<n_bonded_ia;i++) {
    if (bonded_ia_params[i].type != BONDED_IA_NONE) {
      sprintf(buffer, "%d ", i);
      Tcl_AppendResult(interp, "{ ", buffer, get_name_of_bonded_ia(bonded_ia_params[i].type)," ", (char *)NULL);
      for(j=0; j<9; j++) {
	Tcl_PrintDouble(interp, obsstat_bonded(&total_p_tensor, i)[j], buffer);
	Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
      }
      Tcl_AppendResult(interp, "} ", (char *)NULL);
    }
  }

  for (i = 0; i < n_particle_types; i++)
    for (j = i; j < n_particle_types; j++) {
      if (checkIfParticlesInteract(i, j)) {
	sprintf(buffer, "%d ", i);
	Tcl_AppendResult(interp, "{ ", buffer, (char *)NULL);
	sprintf(buffer, "%d ", j);
	Tcl_AppendResult(interp, " ", buffer, "nonbonded ", (char *)NULL);
	for(k=0; k<9; k++) {
	  Tcl_PrintDouble(interp, obsstat_nonbonded(&total_p_tensor, i, j)[k], buffer);
	  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
	}
	Tcl_AppendResult(interp, "} ", (char *)NULL);
      }
    }

  Tcl_AppendResult(interp, " { total_nb_intra ", (char *)NULL);
  for(k=0; k<9; k++) {
    tvalue[k] = 0.0;
    for (i = 0; i < n_particle_types; i++)
      for (j = i; j < n_particle_types; j++) {
        tvalue[k] += obsstat_nonbonded_intra(&total_p_tensor_non_bonded, i, j)[k];
      }
    Tcl_PrintDouble(interp, tvalue[k], buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  }
  Tcl_AppendResult(interp, "} ", (char *)NULL);
  
  Tcl_AppendResult(interp, " { total_nb_inter ", (char *)NULL);
  for(k=0; k<9; k++) {
    tvalue[k] = 0.0;
    for (i = 0; i < n_particle_types; i++)
      for (j = i; j < n_particle_types; j++) {
        tvalue[k] += obsstat_nonbonded_inter(&total_p_tensor_non_bonded, i, j)[k];
      }
    Tcl_PrintDouble(interp, tvalue[k], buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  }
  Tcl_AppendResult(interp, "} ", (char *)NULL);
 
  for (i = 0; i < n_particle_types; i++)
    for (j = i; j < n_particle_types; j++) {
      if (checkIfParticlesInteract(i, j)) {
	sprintf(buffer, "%d ", i);
	Tcl_AppendResult(interp, "{ ", buffer, (char *)NULL);
	sprintf(buffer, "%d ", j);
	Tcl_AppendResult(interp, " ", buffer, "nb_intra_tensor ", (char *)NULL);
	for(k=0; k<9; k++) {
	  Tcl_PrintDouble(interp, obsstat_nonbonded_intra(&total_p_tensor_non_bonded, i, j)[k], buffer);
	  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
	}
	Tcl_AppendResult(interp, "} ", (char *)NULL);
      }
    }
  
  for (i = 0; i < n_particle_types; i++)
    for (j = i; j < n_particle_types; j++) {
      if (checkIfParticlesInteract(i, j)) {
	sprintf(buffer, "%d ", i);
	Tcl_AppendResult(interp, "{ ", buffer, (char *)NULL);
	sprintf(buffer, "%d ", j);
	Tcl_AppendResult(interp, " ", buffer, "nb_inter_tensor ", (char *)NULL);
	for(k=0; k<9; k++) {
	  Tcl_PrintDouble(interp, obsstat_nonbonded_inter(&total_p_tensor_non_bonded, i, j)[k], buffer);
	  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
	}
	Tcl_AppendResult(interp, "} ", (char *)NULL);
      }
    }

#ifdef ELECTROSTATICS
  if(coulomb.method != COULOMB_NONE) {
    Tcl_AppendResult(interp, "{ coulomb ", (char *)NULL);
    for(j=0; j<9; j++) {
      Tcl_PrintDouble(interp, total_p_tensor.coulomb[j], buffer);
      Tcl_AppendResult(interp, buffer, (char *)NULL);
    }
    Tcl_AppendResult(interp, "} ", (char *)NULL);
  }
#endif

#ifdef MAGNETOSTATICS
  if(coulomb.Dmethod != DIPOLAR_NONE) {
    fprintf(stderr,"tensor magnetostatics, something should go here, file pressure.c ... \n");
  }  
#endif

  
}

/************************************************************/
int parse_and_print_stress_tensor(Tcl_Interp *interp, int v_comp, int argc, char **argv)
{
  /* 'analyze stress_tensor [{ bond <type_num> | nonbonded <type1> <type2> | coulomb | ideal | total }]' */
  char buffer[TCL_DOUBLE_SPACE + TCL_INTEGER_SPACE + 2];
  int i, j, k;
  double p_vel[3], tvalue[9];
  for(j=0; j<9; j++)  tvalue[j] = 0.0;

  if (n_total_particles == 0) {
    Tcl_AppendResult(interp, "(no particles)",
		     (char *)NULL);
    return (TCL_OK);
  }

  /* if desired (v_comp==1) replace ideal component with instantaneous one */
   if (total_pressure.init_status != 1+v_comp ) {
    init_virials(&total_pressure);
    init_p_tensor(&total_p_tensor);

    init_virials_non_bonded(&total_pressure_non_bonded);
    init_p_tensor_non_bonded(&total_p_tensor_non_bonded);

    if(v_comp && (integ_switch == INTEG_METHOD_NPT_ISO) && !(nptiso.invalidate_p_vel)) {
      if (total_pressure.init_status == 0)
	master_pressure_calc(0);
      p_tensor.data.e[0] = 0.0;
      MPI_Reduce(nptiso.p_vel, p_vel, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      for(i=0; i<3; i++)
	if(nptiso.geometry & nptiso.nptgeom_dir[i])
	  p_tensor.data.e[0] += p_vel[i];
      p_tensor.data.e[0] /= (nptiso.dimension*nptiso.volume);
      total_pressure.init_status = 1+v_comp;   }
    else
      master_pressure_calc(v_comp);
  }

  if (argc == 0)
    print_detailed_stress_tensor(interp);
  else {

    if      (ARG0_IS_S("ideal")) {
      for(j=0; j<9; j++)  tvalue[j] = total_p_tensor.data.e[j];
    }
    else if (ARG0_IS_S("bonded") ||
	     ARG0_IS_S("fene") ||
	     ARG0_IS_S("subt_lj_harm") ||
	     ARG0_IS_S("subt_lj_fene") ||
	     ARG0_IS_S("subt_lj") ||
	     ARG0_IS_S("harmonic")) {
      if(argc<2 || ! ARG1_IS_I(i)) {
	Tcl_ResetResult(interp);
	Tcl_AppendResult(interp, "wrong # or type of arguments for: analyze pressure bonded <type_num>",
			 (char *)NULL);
	return (TCL_ERROR);
      }
      if(i < 0 || i >= n_bonded_ia) {
	Tcl_AppendResult(interp, "bond type does not exist", (char *)NULL);
	return (TCL_ERROR);
      }
      for(k=0; k<9; k++) tvalue[k] = obsstat_bonded(&total_p_tensor, i)[k];
    }
    else if (ARG0_IS_S("nonbonded") ||
	     ARG0_IS_S("lj") ||
	     ARG0_IS_S("buckingham") ||
             ARG0_IS_S("morse") ||
             ARG0_IS_S("soft-sphere") ||
	     ARG0_IS_S("lj-cos") ||
	     ARG0_IS_S("lj-cos2") ||
	     ARG0_IS_S("tabulated") ||
	     ARG0_IS_S("gb")) {
      if(argc<3 || ! ARG_IS_I(1, i) || ! ARG_IS_I(2, j)) {
	Tcl_ResetResult(interp);
	Tcl_AppendResult(interp, "wrong # or type of arguments for: analyze pressure nonbonded <type1> <type2>",
			 (char *)NULL);
	return (TCL_ERROR);
      }
      if(i < 0 || i >= n_particle_types || j < 0 || j >= n_particle_types) {
	Tcl_AppendResult(interp, "particle type does not exist", (char *)NULL);
	return (TCL_ERROR);
      }
      for(k=0; k<9; k++) tvalue[k] = obsstat_nonbonded(&total_p_tensor, i, j)[k];
    }
    else if( ARG0_IS_S("tot_nb_intra")) {
      for(k=0; k<9; k++) {
        for (i = 0; i < n_particle_types; i++)
          for (j = i; j < n_particle_types; j++) {
            tvalue[k] += obsstat_nonbonded_intra(&total_p_tensor_non_bonded, i, j)[k];
          }
      }
    }
    else if( ARG0_IS_S("tot_nb_inter")) {
      for(k=0; k<9; k++) {
        for (i = 0; i < n_particle_types; i++)
          for (j = i; j < n_particle_types; j++) {
            tvalue[k] += obsstat_nonbonded_inter(&total_p_tensor_non_bonded, i, j)[k];
          }
      }
    }
    else if( ARG0_IS_S("nb_intra")) {
      if(argc<3 || ! ARG_IS_I(1, i) || ! ARG_IS_I(2, j)) {
	Tcl_ResetResult(interp);
	Tcl_AppendResult(interp, "wrong # or type of arguments for: analyze stress tensor nonbonded <type1> <type2>",
			 (char *)NULL);
	return (TCL_ERROR);
      }
      if(i < 0 || i >= n_particle_types || j < 0 || j >= n_particle_types) {
	Tcl_AppendResult(interp, "particle type does not exist", (char *)NULL);
	return (TCL_ERROR);
      }
      for(k=0; k<9; k++) tvalue[k] = obsstat_nonbonded_intra(&total_p_tensor_non_bonded, i, j)[k];
    }
    else if( ARG0_IS_S("nb_inter")) {
      if(argc<3 || ! ARG_IS_I(1, i) || ! ARG_IS_I(2, j)) {
	Tcl_ResetResult(interp);
	Tcl_AppendResult(interp, "wrong # or type of arguments for: analyze stress tensor nonbonded <type1> <type2>",
			 (char *)NULL);
	return (TCL_ERROR);
      }
      if(i < 0 || i >= n_particle_types || j < 0 || j >= n_particle_types) {
	Tcl_AppendResult(interp, "particle type does not exist", (char *)NULL);
	return (TCL_ERROR);
      }
      for(k=0; k<9; k++) tvalue[k] = obsstat_nonbonded_inter(&total_p_tensor_non_bonded, i, j)[k];
    }
    else if( ARG0_IS_S("coulomb")) {
#ifdef ELECTROSTATICS
      for(j=0; j<9; j++) tvalue[j] = total_p_tensor.coulomb[j];
#else
      Tcl_AppendResult(interp, "ELECTROSTATICS not compiled (see config.h)\n", (char *)NULL);
#endif
    }
    else if( ARG0_IS_S("dipolar")) {
#ifdef MAGNETOSTATICS
      /* for(j=0; j<9; j++) tvalue[j] = total_p_tensor.coulomb[j];*/
      fprintf(stderr," stress tensor, magnetostatics, something should go here, file pressure.c ");
#else
      Tcl_AppendResult(interp, "MAGNETOSTATICS not compiled (see config.h)\n", (char *)NULL);
#endif
    }
    else if (ARG0_IS_S("total")) {
      for(j=0; j<9; j++) {
        tvalue[j] = p_tensor.data.e[j];
        for (i = 1; i < p_tensor.data.n/9; i++) tvalue[j] += total_p_tensor.data.e[9*i + j];
     }
    }
    else {
      Tcl_AppendResult(interp, "unknown feature of: analyze pressure",
		       (char *)NULL);
      return (TCL_ERROR);
    }

    Tcl_AppendResult(interp, *argv, (char *)NULL);
    Tcl_AppendResult(interp, " ", (char *)NULL);
    for(j=0; j<9; j++) {
      Tcl_PrintDouble(interp, tvalue[j], buffer);
      Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
    }
  }

  return (TCL_OK);
}

int parse_local_stress_tensor(Tcl_Interp *interp, int argc, char **argv)
{
  char buffer[TCL_DOUBLE_SPACE];
  char* usage;
  int periodic[3];
  double range_start[3];
  double range[3];
  int bins[3];
  int i,j,k,l;
  DoubleList *TensorInBin;
  PTENSOR_TRACE(fprintf(stderr,"%d: Running parse_local_stress_tensor\n",this_node));
  /* 'analyze stress profile ' */
  if (argc != 12) {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "local_stress_tensor requires 12 inputs: x_periodic, y_periodic, z_periodic, x_range_start, y_range_start, z_range_start, x_range, y_range, z_range, x_bins, y_bins, z_bins", (char *)NULL);
    return(TCL_ERROR);
  }
  usage = "usage: analyse local_stress_tensor <x_periodic> <y_periodic> <z_periodic> <x_range_start> <y_range_start> <z_range_start> <x_range> <y_range> <z_range> <x_bins> <y_bins> <z_bins>";
 
  for (i=0;i<3;i++) {
    if ( !ARG0_IS_I(periodic[i]) ) {
      Tcl_ResetResult(interp);
      Tcl_AppendResult(interp,usage, (char *)NULL);
      return (TCL_ERROR);
    } else {
      argc -= 1;
      argv += 1;
    }
  }
  for (i=0;i<3;i++) {
    if ( !ARG0_IS_D(range_start[i]) ) {
      Tcl_ResetResult(interp);
      Tcl_AppendResult(interp,usage, (char *)NULL);
      return (TCL_ERROR);
    } else {
      argc -= 1;
      argv += 1;
    }
  }
  for (i=0;i<3;i++) {
    if ( !ARG0_IS_D(range[i]) ) {
      Tcl_ResetResult(interp);
      Tcl_AppendResult(interp,usage, (char *)NULL);
      return (TCL_ERROR);
    } else {
      argc -= 1;
      argv += 1;
    }
  }
  for (i=0;i<3;i++) {
    if ( !ARG0_IS_I(bins[i]) ) {
      Tcl_ResetResult(interp);
      Tcl_AppendResult(interp,usage, (char *)NULL);
      return (TCL_ERROR);
    } else {
      argc -= 1;
      argv += 1;
    }
  }

  /* Allocate a doublelist of bins to keep track of stress profile */
  TensorInBin = (DoubleList *)malloc(bins[0]*bins[1]*bins[2]*sizeof(DoubleList));
  if ( TensorInBin ) {
  /* Initialize the stress profile */
    for ( i = 0 ; i < bins[0]*bins[1]*bins[2]; i++ ) {
      init_doublelist(&TensorInBin[i]);
      alloc_doublelist(&TensorInBin[i],9);
      for ( j = 0 ; j < 9 ; j++ ) {
	TensorInBin[i].e[j] = 0.0;
      }
    }
  } else {
    Tcl_AppendResult(interp,"could not allocate memory for local_stress_tensor", (char *)NULL);
    return (TCL_ERROR);
  }

  mpi_local_stress_tensor(TensorInBin, bins, periodic,range_start, range);
  PTENSOR_TRACE(fprintf(stderr,"%d: parse_local_stress_tensor: finished mpi_local_stress_tensor \n",this_node));

  /* Write stress profile to Tcl export */
  Tcl_AppendResult(interp, "{ LocalStressTensor } ", (char *)NULL);
  for ( i = 0 ; i < bins[0] ; i++) {
    for ( j = 0 ; j < bins[1] ; j++) {
      for ( k = 0 ; k < bins[2] ; k++) {
	Tcl_AppendResult(interp, " { ", (char *)NULL);
	sprintf(buffer," { %d %d %d } ",i,j,k);
	Tcl_AppendResult(interp,buffer, (char *)NULL);
	Tcl_AppendResult(interp, " { ", (char *)NULL);
	for ( l = 0 ; l < 9 ; l++) {
	  Tcl_PrintDouble(interp,TensorInBin[i*bins[1]*bins[2]+j*bins[2]+k].e[l],buffer);
	  Tcl_AppendResult(interp, buffer, (char *)NULL);
	  Tcl_AppendResult(interp, " ", (char *)NULL);
	}
	Tcl_AppendResult(interp, " } ", (char *)NULL);
	Tcl_AppendResult(interp, " } ", (char *)NULL);
      }
    }
  }
  
  /* Free memory */
  for ( i = 0 ; i < bins[0]*bins[1]*bins[2] ; i++ ) {
    realloc_doublelist(&TensorInBin[i],0);
  }
  free(TensorInBin);
  return TCL_OK;
}
