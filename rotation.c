// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2009; all rights reserved unless otherwise stated.
/** \file rotation.c  Molecular dynamics integrator for rotational motion.
 *
 *  A velocity Verlet <a HREF="http://ciks.cbt.nist.gov/~garbocz/dpd1/dpd.html">algorithm</a>
 *  using quaternions is implemented to tackle rotational motion. A random torque and a friction
 *  term are added to provide the constant NVT conditions. Due to this feature all particles are
 *  treated as 3D objects with 3 translational and 3 rotational degrees of freedom if ROTATION
 *  flag is set in \ref config.h "config.h".
*/

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "utils.h"
#include "integrate.h"
#include "interaction_data.h"
#include "particle_data.h"
#include "communication.h"
#include "grid.h"
#include "cells.h"
#include "verlet.h"
#include "rotation.h"
#include "ghosts.h"
#include "forces.h"
#include "p3m.h"
#include "thermostat.h"
#include "initialize.h"

/****************************************************
 *                     DEFINES
 ***************************************************/
/**************** local variables  *******************/

#ifdef ROTATION

/** moment of inertia. Currently we define the inertia tensor here to be constant.
    If it is not spherical the angular velocities have to be refined several times
    in the \ref convert_torqes_propagate_omega. Also the kinetic energy in file
    \ref statistics.c is calculated assuming that I[0] =  I[1] =  I[2] = 1  */
static double I[3] = { 1, 1, 1};

/** \name Privat Functions */
/************************************************************/
/*@{*/

/** define rotation matrix A for a given particle */
static void define_rotation_matrix(Particle *p, double A[9]);

/** define first and second time derivatives of a quaternion, as well
    as the angular acceleration. */
static void define_Qdd(Particle *p, double Qd[4], double Qdd[4], double S[3], double Wd[3]);

/*@}*/

/** convert quaternions to the director */
void convert_quat_to_quatu(double quat[4], double quatu[3])
{
  /* director */
  quatu[0] = 2*(quat[1]*quat[3] + quat[0]*quat[2]);
  quatu[1] = 2*(quat[2]*quat[3] - quat[0]*quat[1]);
  quatu[2] =   (quat[0]*quat[0] - quat[1]*quat[1] -
		quat[2]*quat[2] + quat[3]*quat[3]); 
}

#ifdef DIPOLES

/** convert a dipole moment to quaternions and dipolar strength  */
int convert_dip_to_quat(double dip[3], double quat[4], double *dipm)
{
  double dip_xy, dm;
  double theta2, phi2;

  // Calculate magnitude of dipole moment
  dm = sqrt(dip[0]*dip[0] + dip[1]*dip[1] + dip[2]*dip[2]);
  *dipm = dm;

  // The dipole vector needs to be != 0 to be converted into a quaternion
  if (dm < ROUND_ERROR_PREC) {
    return 1;
  }
  else {
    // Calculate angles 
    dip_xy = sqrt(dip[0]*dip[0] + dip[1]*dip[1]);
    // If dipole points along z axis:
    if (dip_xy == 0){
      // We need to distinguish between (0,0,d_z) and (0,0,d_z)
      if (dip[2]>0)
       theta2 = 0;
      else
       theta2 = PI/2.;
      phi2 = 0;
    }
    else {
      // Here, we take care of all other directions
      //Here we suppose that theta2 = 0.5*theta and phi2 = 0.5*(phi - PI/2),
      //where theta and phi - angles are in spherical coordinates
      theta2 = 0.5*acos(dip[2]/dm);
      if (dip[1] < 0) phi2 = -0.5*acos(dip[0]/dip_xy) - PI*0.25;
      else phi2 = 0.5*acos(dip[0]/dip_xy) - PI*0.25;
    }

    // Calculate the quaternion from the angles
    quat[0] =  cos(theta2) * cos(phi2);
    quat[1] = -sin(theta2) * cos(phi2);
    quat[2] = -sin(theta2) * sin(phi2);
    quat[3] =  cos(theta2) * sin(phi2);
  }
  return 0;
}

/** convert quaternion director to the dipole moment */
void convert_quatu_to_dip(double quatu[3], double dipm, double dip[3])
{
  /* dipole moment */
  dip[0] = quatu[0]*dipm;
  dip[1] = quatu[1]*dipm;
  dip[2] = quatu[2]*dipm;
}

#endif

/** Here we use quaternions to calculate the rotation matrix which
    will be used then to transform torques from the laboratory to
    the body-fixed frames */  
void define_rotation_matrix(Particle *p, double A[9])
{
  A[0 + 3*0] = p->r.quat[0]*p->r.quat[0] + p->r.quat[1]*p->r.quat[1] -
               p->r.quat[2]*p->r.quat[2] - p->r.quat[3]*p->r.quat[3] ;
  A[1 + 3*1] = p->r.quat[0]*p->r.quat[0] - p->r.quat[1]*p->r.quat[1] +
               p->r.quat[2]*p->r.quat[2] - p->r.quat[3]*p->r.quat[3] ;
  A[2 + 3*2] = p->r.quat[0]*p->r.quat[0] - p->r.quat[1]*p->r.quat[1] -
               p->r.quat[2]*p->r.quat[2] + p->r.quat[3]*p->r.quat[3] ;

  A[0 + 3*1] = 2*(p->r.quat[1]*p->r.quat[2] + p->r.quat[0]*p->r.quat[3]);
  A[0 + 3*2] = 2*(p->r.quat[1]*p->r.quat[3] - p->r.quat[0]*p->r.quat[2]);
  A[1 + 3*0] = 2*(p->r.quat[1]*p->r.quat[2] - p->r.quat[0]*p->r.quat[3]);

  A[1 + 3*2] = 2*(p->r.quat[2]*p->r.quat[3] + p->r.quat[0]*p->r.quat[1]);
  A[2 + 3*0] = 2*(p->r.quat[1]*p->r.quat[3] + p->r.quat[0]*p->r.quat[2]);
  A[2 + 3*1] = 2*(p->r.quat[2]*p->r.quat[3] - p->r.quat[0]*p->r.quat[1]);
}

/** calculate the second derivative of the quaternion of a given particle
    as well as Wd vector which is the angular acceleration of this particle */
void define_Qdd(Particle *p, double Qd[4], double Qdd[4], double S[3], double Wd[3])
{
  double S1;

  /* calculate the first derivative of the quaternion */
  Qd[0] = 0.5 * ( -p->r.quat[1] * p->m.omega[0] -
                   p->r.quat[2] * p->m.omega[1] -
		   p->r.quat[3] * p->m.omega[2] );

  Qd[1] = 0.5 * (  p->r.quat[0] * p->m.omega[0] -
                   p->r.quat[3] * p->m.omega[1] +
		   p->r.quat[2] * p->m.omega[2] );

  Qd[2] = 0.5 * (  p->r.quat[3] * p->m.omega[0] +
                   p->r.quat[0] * p->m.omega[1] -
		   p->r.quat[1] * p->m.omega[2] );

  Qd[3] = 0.5 * ( -p->r.quat[2] * p->m.omega[0] +
                   p->r.quat[1] * p->m.omega[1] +
		   p->r.quat[0] * p->m.omega[2] );

  /* calculate the second derivative of the quaternion */
  
#ifdef ROTATIONAL_INERTIA
  Wd[0] =  (p->f.torque[0] + p->m.omega[1]*p->m.omega[2]*(I[1]-I[2]))/I[0]/p->p.rinertia[0];
  Wd[1] =  (p->f.torque[1] + p->m.omega[2]*p->m.omega[0]*(I[2]-I[0]))/I[1]/p->p.rinertia[1];
  Wd[2] =  (p->f.torque[2] + p->m.omega[0]*p->m.omega[1]*(I[0]-I[1]))/I[2]/p->p.rinertia[2];
#else
  Wd[0] =  (p->f.torque[0] + p->m.omega[1]*p->m.omega[2]*(I[1]-I[2]))/I[0];
  Wd[1] =  (p->f.torque[1] + p->m.omega[2]*p->m.omega[0]*(I[2]-I[0]))/I[1];
  Wd[2] =  (p->f.torque[2] + p->m.omega[0]*p->m.omega[1]*(I[0]-I[1]))/I[2];
#endif

  S1 = Qd[0]*Qd[0] + Qd[1]*Qd[1] + Qd[2]*Qd[2] + Qd[3]*Qd[3];

  Qdd[0] = 0.5 * ( -p->r.quat[1] * Wd[0] -
		    p->r.quat[2] * Wd[1] -
	  	    p->r.quat[3] * Wd[2] ) - p->r.quat[0] * S1;

  Qdd[1] = 0.5 * (  p->r.quat[0] * Wd[0] -
		    p->r.quat[3] * Wd[1] +
		    p->r.quat[2] * Wd[2] ) - p->r.quat[1] * S1;
  
  Qdd[2] = 0.5 * (  p->r.quat[3] * Wd[0] +
		    p->r.quat[0] * Wd[1] -
		    p->r.quat[1] * Wd[2] ) - p->r.quat[2] * S1;
  
  Qdd[3] = 0.5 * ( -p->r.quat[2] * Wd[0] +
                    p->r.quat[1] * Wd[1] +
		    p->r.quat[0] * Wd[2] ) - p->r.quat[3] * S1;

  S[0] = S1;
  S[1] = Qd[0]*Qdd[0]  + Qd[1]*Qdd[1]  + Qd[2]*Qdd[2]  + Qd[3]*Qdd[3];
  S[2] = Qdd[0]*Qdd[0] + Qdd[1]*Qdd[1] + Qdd[2]*Qdd[2] + Qdd[3]*Qdd[3];
}

/** propagate angular velocities and quaternions \todo implement for
       fixed_coord_flag */
void propagate_omega_quat()
{
  Particle *p;
  Cell *cell;
  int c,i,j, np;
  double lambda, dt2, dt4, dtdt, dtdt2;

  dt2 = time_step*0.5;
  dt4 = time_step*0.25;
  dtdt = time_step*time_step;
  dtdt2 = dtdt*0.5;

  INTEG_TRACE(fprintf(stderr,"%d: propagate_omega_quat:\n",this_node));

  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    p  = cell->part;
    np = cell->n;
    for(i = 0; i < np; i++) {
	  double Qd[4], Qdd[4], S[3], Wd[3];
	  define_Qdd(&p[i], Qd, Qdd, S, Wd);
	  
	  lambda = 1 - S[0]*dtdt2 - sqrt(1 - dtdt*(S[0] + time_step*(S[1] + dt4*(S[2]-S[0]*S[0]))));

	  for(j=0; j < 3; j++){
	    p[i].m.omega[j]+= dt2*Wd[j];
	  }
	  ONEPART_TRACE(if(p[i].p.identity==check_id) fprintf(stderr,"%d: OPT: PV_1 v_new = (%.3e,%.3e,%.3e)\n",this_node,p[i].m.v[0],p[i].m.v[1],p[i].m.v[2]));

	  p[i].r.quat[0]+= time_step*(Qd[0] + dt2*Qdd[0]) - lambda*p[i].r.quat[0];
	  p[i].r.quat[1]+= time_step*(Qd[1] + dt2*Qdd[1]) - lambda*p[i].r.quat[1];
	  p[i].r.quat[2]+= time_step*(Qd[2] + dt2*Qdd[2]) - lambda*p[i].r.quat[2];
	  p[i].r.quat[3]+= time_step*(Qd[3] + dt2*Qdd[3]) - lambda*p[i].r.quat[3];

	  convert_quat_to_quatu(p[i].r.quat, p[i].r.quatu);
#ifdef DIPOLES
	  convert_quatu_to_dip(p[i].r.quatu, p[i].p.dipm, p[i].r.dip);
#endif
	

      ONEPART_TRACE(if(p[i].p.identity==check_id) fprintf(stderr,"%d: OPT: PPOS p = (%.3f,%.3f,%.3f)\n",this_node,p[i].r.p[0],p[i].r.p[1],p[i].r.p[2]));
    }
  }
}

/** convert the torques to the body-fixed frames and propagate angular velocities */
void convert_torqes_propagate_omega()
{
  Particle *p;
  Cell *cell;
  int c,i, np, times;
  double dt2, tx, ty, tz;

  dt2 = time_step*0.5;
  INTEG_TRACE(fprintf(stderr,"%d: convert_torqes_propagate_omega:\n",this_node));
  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    p  = cell->part;
    np = cell->n;
    for(i = 0; i < np; i++) {
      double A[9];
      define_rotation_matrix(&p[i], A);

      tx = A[0 + 3*0]*p[i].f.torque[0] + A[0 + 3*1]*p[i].f.torque[1] + A[0 + 3*2]*p[i].f.torque[2];
      ty = A[1 + 3*0]*p[i].f.torque[0] + A[1 + 3*1]*p[i].f.torque[1] + A[1 + 3*2]*p[i].f.torque[2];
      tz = A[2 + 3*0]*p[i].f.torque[0] + A[2 + 3*1]*p[i].f.torque[1] + A[2 + 3*2]*p[i].f.torque[2];

      if ( thermo_switch & THERMO_LANGEVIN ) {
	friction_thermo_langevin_rotation(&p[i]);

	p[i].f.torque[0]+= tx;
	p[i].f.torque[1]+= ty;
	p[i].f.torque[2]+= tz;
      } else {
	p[i].f.torque[0] = tx;
	p[i].f.torque[1] = ty;
	p[i].f.torque[2] = tz;
      }
    
      ONEPART_TRACE(if(p[i].p.identity==check_id) fprintf(stderr,"%d: OPT: SCAL f = (%.3e,%.3e,%.3e) v_old = (%.3e,%.3e,%.3e)\n",this_node,p[i].f.f[0],p[i].f.f[1],p[i].f.f[2],p[i].m.v[0],p[i].m.v[1],p[i].m.v[2]));
    
#ifdef ROTATIONAL_INERTIA
	  p[i].m.omega[0]+= dt2*p[i].f.torque[0]/p[i].p.rinertia[0]/I[0];
	  p[i].m.omega[1]+= dt2*p[i].f.torque[1]/p[i].p.rinertia[1]/I[1];
	  p[i].m.omega[2]+= dt2*p[i].f.torque[2]/p[i].p.rinertia[2]/I[2];
#else
	  p[i].m.omega[0]+= dt2*p[i].f.torque[0]/I[0];
	  p[i].m.omega[1]+= dt2*p[i].f.torque[1]/I[1];
	  p[i].m.omega[2]+= dt2*p[i].f.torque[2]/I[2];
#endif
	  /* if the tensor of inertia is isotrpic, the following refinement is not needed.
	     Otherwise repeat this loop 2-3 times depending on the required accuracy */
	  for(times=0;times<=0;times++) { 
	    double Wd[3];

#ifdef ROTATIONAL_INERTIA
	    Wd[0] = (p[i].m.omega[1]*p[i].m.omega[2]*(I[1]-I[2]))/I[0]/p[i].p.rinertia[0];
	    Wd[1] = (p[i].m.omega[2]*p[i].m.omega[0]*(I[2]-I[0]))/I[1]/p[i].p.rinertia[1]; 
	    Wd[2] = (p[i].m.omega[0]*p[i].m.omega[1]*(I[0]-I[1]))/I[2]/p[i].p.rinertia[2];
#else
	    Wd[0] = (p[i].m.omega[1]*p[i].m.omega[2]*(I[1]-I[2]))/I[0];
	    Wd[1] = (p[i].m.omega[2]*p[i].m.omega[0]*(I[2]-I[0]))/I[1]; 
	    Wd[2] = (p[i].m.omega[0]*p[i].m.omega[1]*(I[0]-I[1]))/I[2];
#endif
 
	    p[i].m.omega[0]+= dt2*Wd[0];
	    p[i].m.omega[1]+= dt2*Wd[1];
	    p[i].m.omega[2]+= dt2*Wd[2];
	  }
	
      
      ONEPART_TRACE(if(p[i].p.identity==check_id) fprintf(stderr,"%d: OPT: PV_2 v_new = (%.3e,%.3e,%.3e)\n",this_node,p[i].m.v[0],p[i].m.v[1],p[i].m.v[2]));
      
    }
  }
}

/** convert the torques to the body-fixed frames before the integration loop */
void convert_initial_torques()
{
  Particle *p;
  Cell *cell;
  int c,i, np;
  double tx, ty, tz;

  INTEG_TRACE(fprintf(stderr,"%d: convert_initial_torqes:\n",this_node));
  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    p  = cell->part;
    np = cell->n;
    for(i = 0; i < np; i++) {
      double A[9];
      define_rotation_matrix(&p[i], A);

      tx = A[0 + 3*0]*p[i].f.torque[0] + A[0 + 3*1]*p[i].f.torque[1] + A[0 + 3*2]*p[i].f.torque[2];
      ty = A[1 + 3*0]*p[i].f.torque[0] + A[1 + 3*1]*p[i].f.torque[1] + A[1 + 3*2]*p[i].f.torque[2];
      tz = A[2 + 3*0]*p[i].f.torque[0] + A[2 + 3*1]*p[i].f.torque[1] + A[2 + 3*2]*p[i].f.torque[2];

      if ( thermo_switch & THERMO_LANGEVIN ) {
	friction_thermo_langevin_rotation(&p[i]);
	p[i].f.torque[0]+= tx;
	p[i].f.torque[1]+= ty;
	p[i].f.torque[2]+= tz;
      } else {
	p[i].f.torque[0] = tx;
	p[i].f.torque[1] = ty;
	p[i].f.torque[2] = tz;
      }

      ONEPART_TRACE(if(p[i].p.identity==check_id) fprintf(stderr,"%d: OPT: SCAL f = (%.3e,%.3e,%.3e) v_old = (%.3e,%.3e,%.3e)\n",this_node,p[i].f.f[0],p[i].f.f[1],p[i].f.f[2],p[i].m.v[0],p[i].m.v[1],p[i].m.v[2]));
    }
  }
}

/** convert torques from the body-fixed frames to space-fixed coordinates */
void convert_torques_body_to_space(Particle *p, double *torque)
{
  double A[9];
  define_rotation_matrix(p, A);

  torque[0] = A[0 + 3*0]*p->f.torque[0] + A[1 + 3*0]*p->f.torque[1] + A[2 + 3*0]*p->f.torque[2];
  torque[1] = A[0 + 3*1]*p->f.torque[0] + A[1 + 3*1]*p->f.torque[1] + A[2 + 3*1]*p->f.torque[2];
  torque[2] = A[0 + 3*2]*p->f.torque[0] + A[1 + 3*2]*p->f.torque[1] + A[2 + 3*2]*p->f.torque[2];

}

#endif
