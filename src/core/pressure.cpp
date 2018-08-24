/*
  Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
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
/** \file pressure.cpp
    Implementation of \ref pressure.hpp "pressure.h".
*/

#include "pressure_inline.hpp"
#include "cells.hpp"
#include "integrate.hpp"
#include "initialize.hpp"
#include "virtual_sites.hpp" 
#include "npt.hpp"
#include "p3m.hpp"
#include "p3m-dipolar.hpp"

#include "short_range_loop.hpp"

Observable_stat virials  = {0, {}, 0,0,0,0,0};
Observable_stat total_pressure = {0, {}, 0,0,0,0,0};
Observable_stat p_tensor = {0, {},0,0,0,0,0};
Observable_stat total_p_tensor = {0, {},0,0,0,0,0};

/* Observables used in the calculation of intra- and inter- molecular
   non-bonded contributions to pressure and to stress tensor */
Observable_stat_non_bonded virials_non_bonded  = {0, {}, 0,0,0};
Observable_stat_non_bonded total_pressure_non_bonded = {0, {}, 0,0,0};
Observable_stat_non_bonded p_tensor_non_bonded = {0, {},0,0,0};
Observable_stat_non_bonded total_p_tensor_non_bonded = {0, {},0,0,0};

nptiso_struct   nptiso   = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,{0.0,0.0,0.0},{0.0,0.0,0.0},1, 0 ,{NPTGEOM_XDIR, NPTGEOM_YDIR, NPTGEOM_ZDIR},0,0,0};

/************************************************************/
/* callbacks for setmd                                      */
/************************************************************/



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

/** Initializes stat to be used by \ref pressure_calc. */
void init_p_tensor(Observable_stat *stat);

/** Initializes stat_nb to be used by \ref pressure_calc. */
void init_p_tensor_non_bonded(Observable_stat_non_bonded *stat_nb);

/*********************************/
/* Scalar and Tensorial Pressure */
/*********************************/

void pressure_calc(double *result, double *result_t, double *result_nb, double *result_t_nb, int v_comp)
{
  int n, i;
  double volume = box_l[0]*box_l[1]*box_l[2];

  if (!interactions_sanity_checks())
    return;

  init_virials(&virials);

  init_p_tensor(&p_tensor);
  
  init_virials_non_bonded(&virials_non_bonded);

  init_p_tensor_non_bonded(&p_tensor_non_bonded);

  on_observable_calc();

  short_range_loop(
      [&v_comp](Particle &p) {
        add_kinetic_virials(&p, v_comp);
        add_bonded_virials(&p);
#ifdef BOND_ANGLE_OLD
        add_three_body_bonded_stress(&p);
#endif
#ifdef BOND_ANGLE
        add_three_body_bonded_stress(&p);
#endif
      },
      [](Particle &p1, Particle &p2, Distance &d) {
          add_non_bonded_pair_virials(&(p1), &(p2), d.vec21.data(), sqrt(d.dist2),
                                      d.dist2);
      });

/* rescale kinetic energy (=ideal contribution) */
  virials.data.e[0] /= (3.0*volume*time_step*time_step);

  calc_long_range_virials();

#ifdef VIRTUAL_SITES
  virtual_sites()->pressure_and_stress_tensor_contribution(virials.virtual_sites,p_tensor.virtual_sites);
#endif


  for (n = 1; n < virials.data.n; n++)
    virials.data.e[n] /= 3.0*volume;

    
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
  MPI_Reduce(virials.data.e, result, virials.data.n, MPI_DOUBLE, MPI_SUM, 0, comm_cart);
  MPI_Reduce(p_tensor.data.e, result_t, p_tensor.data.n, MPI_DOUBLE, MPI_SUM, 0, comm_cart);
  
  MPI_Reduce(virials_non_bonded.data_nb.e, result_nb, virials_non_bonded.data_nb.n, MPI_DOUBLE, MPI_SUM, 0, comm_cart);
  MPI_Reduce(p_tensor_non_bonded.data_nb.e, result_t_nb, p_tensor_non_bonded.data_nb.n, MPI_DOUBLE, MPI_SUM, 0, comm_cart); 
}

/************************************************************/

void calc_long_range_virials()
{
#ifdef ELECTROSTATICS
  /* calculate k-space part of electrostatic interaction. */
  switch (coulomb.method) {
#ifdef P3M
  case COULOMB_ELC_P3M:
    fprintf(stderr, "WARNING: pressure calculated, but ELC pressure not implemented\n");
    break;
  case COULOMB_P3M_GPU:
    fprintf(stderr, "WARNING: pressure calculated, but GPU P3M pressure not implemented\n");
    break;
  case COULOMB_P3M: {
    p3m_charge_assign();
    virials.coulomb[1] = p3m_calc_kspace_forces(0,1);
    p3m_charge_assign();
    p3m_calc_kspace_stress(p_tensor.coulomb + 9);
    break;
  }
#endif
  case COULOMB_MMM2D:
    fprintf(stderr, "WARNING: pressure calculated, but MMM2D pressure not implemented\n");
    break;
  case COULOMB_MMM1D:
  case COULOMB_MMM1D_GPU:
    fprintf(stderr, "WARNING: pressure calculated, but MMM1D pressure not implemented\n");
    break;
  default:
	  break;
  }
#endif /*ifdef ELECTROSTATICS */  
  
#ifdef DIPOLES
  /* calculate k-space part of magnetostatic interaction. */
  switch (coulomb.Dmethod) {
     case DIPOLAR_ALL_WITH_ALL_AND_NO_REPLICA:
    fprintf(stderr, "WARNING: pressure calculated, but DAWAANR pressure not implemented\n");
    break;
  case DIPOLAR_MDLC_DS:
    fprintf(stderr, "WARNING: pressure calculated, but DLC pressure not implemented\n");
    break;
  case DIPOLAR_DS:
    fprintf(stderr, "WARNING: pressure calculated, but  MAGNETIC DIRECT SUM pressure not implemented\n");
    break;

   
#ifdef DP3M
  case DIPOLAR_MDLC_P3M:
    fprintf(stderr, "WARNING: pressure calculated, but DLC pressure not implemented\n");
    break;
  case DIPOLAR_P3M: {
    int k;
    dp3m_dipole_assign();
    virials.dipolar[1] = dp3m_calc_kspace_forces(0,1);
     
    for(k=0;k<3;k++)
      p_tensor.coulomb[9+ k*3 + k] = virials.dipolar[1]/3.;
    fprintf(stderr, "WARNING: stress tensor calculated, but dipolar P3M stress tensor not implemented\n");
    fprintf(stderr, "WARNING: things have been added to the coulomb virial and p_tensor arrays !!!!!!!\n");
    
    break;
  }
#endif
  default:
      break;
 } 
#endif /*ifdef DIPOLES */
}

/* Initialize the virials used in the calculation of the scalar pressure */
/************************************************************/
void init_virials(Observable_stat *stat)
{
    // Determine number of contribution for different interaction types
    // bonded, nonbonded, coulomb, dipolar, rigid bodies
    int n_pre, n_non_bonded, n_coulomb, n_dipolar,n_vs;

  n_pre        = 1;
  n_non_bonded = (max_seen_particle_type*(max_seen_particle_type+1))/2;

  n_coulomb    = 0;
  n_dipolar    = 0;
  n_vs=0;

#ifdef ELECTROSTATICS
  switch (coulomb.method) {
  case COULOMB_NONE: n_coulomb = 0; break;
  case COULOMB_P3M_GPU:
  case COULOMB_P3M:  n_coulomb = 2; break;
  default: n_coulomb  = 1;
  }
#endif
#ifdef DIPOLES
  switch (coulomb.Dmethod) {
  case DIPOLAR_NONE:  n_dipolar = 0; break;
  case DIPOLAR_ALL_WITH_ALL_AND_NO_REPLICA:  n_dipolar = 0; break;
  case DIPOLAR_DS:  n_dipolar = 0; break;
  case DIPOLAR_P3M:   n_dipolar = 2; break;
  default:
      n_dipolar = 0;
      break;
  }
#endif
#ifdef VIRTUAL_SITES
  n_vs=virtual_sites()->n_pressure_contribs();
#endif


  // Allocate memory for the data
  obsstat_realloc_and_clear(stat, n_pre, bonded_ia_params.size(), n_non_bonded, n_coulomb, n_dipolar, n_vs, 1);
  stat->init_status = 0;
}

/************************************************************/
void init_virials_non_bonded(Observable_stat_non_bonded *stat_nb)
{
  int n_non_bonded;

  n_non_bonded = (max_seen_particle_type*(max_seen_particle_type+1))/2;

  obsstat_realloc_and_clear_non_bonded(stat_nb, n_non_bonded, 1);
  stat_nb->init_status_nb = 0;
}

/* Initialize the p_tensor */
/***************************/
void init_p_tensor(Observable_stat *stat)
{
    // Determine number of contribution for different interaction types
    // bonded, nonbonded, coulomb, dipolar, rigid bodies
    int n_pre, n_non_bonded, n_coulomb, n_dipolar,n_vs;


  n_pre        = 1;
  n_non_bonded = (max_seen_particle_type*(max_seen_particle_type+1))/2;

  n_coulomb = 0;
  n_dipolar = 0;
  n_vs=0;

#ifdef ELECTROSTATICS
  switch (coulomb.method) {
  case COULOMB_NONE: n_coulomb = 0; break;
  case COULOMB_P3M_GPU:
  case COULOMB_P3M:  n_coulomb = 2; break;
  default: n_coulomb  = 1;
  }
#endif
  
#ifdef DIPOLES
  switch (coulomb.Dmethod) {
  case DIPOLAR_NONE: n_dipolar = 0; break;
  case DIPOLAR_ALL_WITH_ALL_AND_NO_REPLICA:  n_dipolar = 0; break;
  case DIPOLAR_DS:  n_dipolar = 0; break;
  case DIPOLAR_P3M:  n_dipolar = 2; break;
  default: n_dipolar = 0;
  }
#endif
#ifdef VIRTUAL_SITES
  n_vs=virtual_sites()->n_pressure_contribs();
#endif

  obsstat_realloc_and_clear(stat, n_pre, bonded_ia_params.size(), n_non_bonded, n_coulomb, n_dipolar, n_vs, 9);
  stat->init_status = 0;
}

/***************************/
void init_p_tensor_non_bonded(Observable_stat_non_bonded *stat_nb)
{
  int n_nonbonded;
  n_nonbonded = (max_seen_particle_type*(max_seen_particle_type+1))/2;

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

namespace {
/** Calculates the remainder of a division */
 double drem_down(double a, double b) { return a - floor(a / b) * b; }

int getintersection(double pos1[3], double pos2[3],int given, int get, double value, double *answer, double box_size[3])
{
  /*pos1 and pos2 are two particle positions.                                                  */
  /*given and get are integers from 0 to 2. 0 = x direction. 1 = y direction. 2 = z direction  */
  /*there is a point on the line between the two particles p1 and p2 such that r[given]=value  */
  /*this procedure returns the value of r[get] at that point                                   */

  double p2r[3];
  int i;

  for (i=0;i<3;i++) {
    p2r[i] = drem_down((pos2[i]-pos1[i])+box_size[i]/2.0,box_size[i])-box_size[i]/2.0;
  }
  value = drem_down((value-pos1[given])+box_size[given]/2.0,box_size[given])-box_size[given]/2.0;
  //PTENSOR_TRACE(fprintf(stderr,"%d: getintersection: p1 is %f %f %f p2 is %f %f %f p2r is %f %f %f newvalue is %f\n",this_node,pos1[0],pos1[1],pos1[2],pos2[0],pos2[1],pos2[2],p2r[0],p2r[1],p2r[2],value););
  
  if ((value)*(p2r[given]) < -0.0001) {
      runtimeErrorMsg() <<"analyze stress_profile: getintersection: intersection is not between the two given particles - " << value << " is not between " << 0.0 << " and " << p2r[given] << " and box size is " << box_size[given] << ", given is " << given << "\n";
    return 0; 
  } else if (given == get) {
    *answer =  drem_down(value + pos1[given],box_size[given]);;
  } else if (0==p2r[given]) {
      runtimeErrorMsg() <<"analyze stress_profile: getintersection: intersection is a line, not a point - value is " << value << " same as " << 0.0 << " and " << p2r[given] << "\n";
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
  int totoccupiedybins;     /* total number of y-bins through which the line passes  */
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
    std::vector<int>starty(occupiedxbins+1);
    std::vector<int> occupiedybins(occupiedxbins);

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
    std::vector<int> occupiedzbins(totoccupiedybins);
    std::vector<int> startz(totoccupiedybins+1);
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
        runtimeErrorMsg() << this_node << ": analyze stress_profile: bug in distribute tensor code - calclength is " << calclength << " and length is " << length;
      return 0;
    }
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

int get_nonbonded_interaction(Particle *p1, Particle *p2, double *force, Distance &)
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
    dist2 = Utils::sqr(d[0]) + Utils::sqr(d[1]) + Utils::sqr(d[2]);
    dist  = sqrt(dist2);
    calc_non_bonded_pair_force(p1,p2,d,dist,dist2,force);
#ifdef ELECTROSTATICS
    if (coulomb.method != COULOMB_NONE) {
      switch (coulomb.method) {
#ifdef P3M
      case COULOMB_P3M_GPU:
	fprintf(stderr,"WARNING: Local stress tensor calculation cannot handle GPU P3M electrostatics so it is left out\n");  
	break;
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

#ifdef DIPOLES
    if (coulomb.Dmethod != DIPOLAR_NONE) {
      switch (coulomb.Dmethod) {
#ifdef DP3M
      case DIPOLAR_P3M:
    	fprintf(stderr,"WARNING: Local stress tensor calculation cannot handle P3M magnetostatics so it is left out\n");  
	break;
#endif
      case DIPOLAR_ALL_WITH_ALL_AND_NO_REPLICA:
    	fprintf(stderr,"WARNING: Local stress tensor calculation cannot handle DAWAANR magnetostatics so it is left out\n");  
	break;
      case DIPOLAR_DS:
    	fprintf(stderr,"WARNING: Local stress tensor calculation cannot handle MAGNETIC DIPOLAR SUM magnetostatics so it is left out\n");  
	break;

      default:
	fprintf(stderr,"WARNING: Local stress tensor calculation does not recognise this magnetostatic interaction\n");  
      }
    }
#endif /*ifdef DIPOLES */


  } /*if p1-> ... */
  return 0;
}

} /* namespace */

int local_stress_tensor_calc(DoubleList *TensorInBin, int bins[3],
                             int periodic[3], double range_start[3],
                             double range[3]) {
  /*calculates local stress tensors in cuboid bins
    uses Irving Kirkwood method
    we consider a cube of space starting with a corner at position range_start
    extending to
      range_start + range
    if the variable periodic is set to 1 in dimension i then the cube is assumed
    to span the periodic box
    this cube is divided into bins[0] bins in the x direction bins[1] in the y
    direction, and bins[2] in the z direction
  */

  double binvolume;
  double centre[3];

  for (int i = 0; i < 3; i++) {
    if (periodic[i]) {
      range[i] = box_l[i];
      range_start[i] = 0;
    }
  }

  /* find centre of analyzed cube */
  for (int i = 0; i < 3; i++) {
    centre[i] = range_start[i] + range[i] / 2.0;
  }

  /* We consider all particles that are within a certain distance of the cube.
     The skin is used as this distance.  If the
     skin from on opposite sides of the box overlaps then we produce an error
     message.  To code dround this would be
     creating unnecessary work since I can't imagine when we might want that */

  for (int i = 0; i < 3; i++) {
    if ((!periodic[i]) && (range[i] + 2 * skin + 2 * max_cut > box_l[i])) {
      runtimeErrorMsg() << "analyze stress_profile: Analyzed box (" << range[i]
                        << ") with skin+max_cut(" << skin + max_cut
                        << ") is larger than simulation box (" << box_l[i]
                        << ").\n";
      return 0;
    }
    range_start[i] = drem_down(range_start[i], box_l[i]);
  }
  PTENSOR_TRACE(fprintf(stderr, "%d: Running stress_profile\n", this_node));

  binvolume = range[0] * range[1] * range[2] / (double)bins[0] /
              (double)bins[1] / (double)bins[2];

  auto add_ideal = [&](Particle &p) {
    int bin;
    whichbin(p.r.p.data(), bins, centre, range, &bin);
    if (bin >= 0) {
      for (int k = 0; k < 3; k++) {
        for (int l = 0; l < 3; l++) {
          TensorInBin[bin].e[k * 3 + l] +=
              (p.m.v[k]) * (p.m.v[l]) * p.p.mass / time_step / time_step;
        }
      }
    }
  };

  auto add_bonded =
      [&](Particle &p) {
        int j = 0;
        while (j < p.bl.n) {
          auto type_num = p.bl.e[j++];
          auto iaparams = &bonded_ia_params[type_num];

          /* fetch particle 2 */
          auto p2 = local_particles[p.bl.e[j++]];
          double dx[3];
          get_mi_vector(dx, p.r.p, p2->r.p);
          std::array<double,3> force;
          calc_bonded_force(&p, p2, iaparams, &j, dx, force.data());
          PTENSOR_TRACE(
              fprintf(stderr, "%d: Bonded to particle %d with force %f %f %f\n",
                      this_node, p2->p.identity, force[0], force[1], force[2]));
          if ((pow(force[0], 2) + pow(force[1], 2) + pow(force[2], 2)) > 0) {
            if (distribute_tensors(TensorInBin, force.data(), bins, range_start, range,
                                   p.r.p.data(), p2->r.p.data()) != 1)
              return 0;
          }
        }
        return 0;
  };

  auto add_single_particle_contribution = [&add_ideal, &add_bonded](Particle &p) {
    add_ideal(p);
    add_bonded(p);
  };

  /* this next bit loops over all pair of particles, calculates the force
   * between them, and distributes it amongst the tensors */
  short_range_loop(
      add_single_particle_contribution,
      [&centre, &range, &TensorInBin, &range_start,
       &bins](Particle &p1, Particle &p2, Distance &d) {
        if ((incubewithskin(p1.r.p.data(), centre, range)) &&
            (incubewithskin(p2.r.p.data(), centre, range))) {
          double force[3];
          get_nonbonded_interaction(&p1, &p2, force, d);
          if ((pow(force[0], 2) + pow(force[1], 2) + pow(force[2], 2)) > 0) {
            if (distribute_tensors(TensorInBin, force, bins, range_start, range,
                                   p1.r.p.data(), p2.r.p.data()) != 1)
              return 0;
          }
        }
        return 0;
      });

  for (int i = 0; i < bins[0] * bins[1] * bins[2]; i++) {
    for (int j = 0; j < 9; j++) {
      TensorInBin[i].e[j] /= binvolume;
    }
  }

  return 1;
}

/************************************************************/
int observable_compute_stress_tensor(int v_comp, double *A)
{
  int i, j;
  double value;
  double p_vel[3];

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

  for(j=0; j<9; j++) {
    value = total_p_tensor.data.e[j];
    for (i = 1; i < total_p_tensor.data.n/9; i++) value += total_p_tensor.data.e[9*i + j];
    A[j]=value;
  }
  return 0;
}
void update_stress_tensor (int v_comp) {
	int i;
	double p_vel[3];
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
}

int analyze_local_stress_tensor(int* periodic, double* range_start, double* range, int* bins, DoubleList* TensorInBin)
{
	PTENSOR_TRACE(fprintf(stderr,"%d: Running analyze_local_stress_tensor\n",this_node));

	mpi_local_stress_tensor(TensorInBin, bins, periodic,range_start, range);
	PTENSOR_TRACE(fprintf(stderr,"%d: analyze_local_stress_tensor: finished mpi_local_stress_tensor \n",this_node));

	return ES_OK;
}
