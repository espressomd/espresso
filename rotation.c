/** \file rotation.c  Molecular dynamics integrator for rotational motion.
 *
 *  For more information about the integrator 
 *  see \ref rotation.h "rotation.h".
*/

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "integrate.h"
#include "interaction_data.h"
#include "particle_data.h"
#include "communication.h"
#include "grid.h"
#include "cells.h"
#include "verlet.h"
#include "rotation.h"
#include "ghosts.h"
#include "debug.h"
#include "p3m.h"
#include "utils.h"
#include "thermostat.h"
#include "initialize.h"
#include "forces.h"

/****************************************************
 *                     DEFINES
 ***************************************************/
/**************** local variables  *******************/

#ifdef ROTATION
/** rotation matrix */
static double A[3][3];

/** moment of inertia (at the moment it is defined in \ref init_torques "init_torques()") */
static double I[3];

/** first and second time derivative of a quaternion */
static double Qd[4], Qdd[4];

/** Qd squared */
static double S1;

/** angular accelaration */
static double Wd[3];

#endif

/** \name Privat Functions */
/************************************************************/
/*@{*/

/** define rotation matrix A for a given particle */
void define_rotation_matrix(Particle *p);

/** define the first time derivative of a quaternion */
void define_Qd(Particle *p);

/** define the second time derivative of a quaternion */
void define_Qdd(Particle *p);

/*@}*/

/** set torques to zero for all particles */
void init_torques()
{
#ifdef ROTATION
 Particle *p;
   int np, m, n, o, i;
/* currently we define the inertia tensor here. If it is not spherical
 the angular velocities have to be refined several times in the
 convert_torqes_propagate_omega(). See the bottom of this file.  */
 
 I[0] =  I[1] =  I[2] = 1;

  CELLS_LOOP(m, n, o) {
    p  = CELL_PTR(m, n, o)->pList.part;
    np = CELL_PTR(m, n, o)->pList.n;
    /* all the particles */
    {
      for (i = 0; i < np; i++) {

	p[i].torque[0] = 0;
	p[i].torque[1] = 0;
	p[i].torque[2] = 0;
      }
    }
  }
#endif
}

#ifdef ROTATION
void define_rotation_matrix(Particle *p)
{   
/** Here we use quaternions to calculate the rotation matrix which
    will be used then to transform torques from the laboratory to
    the body-fixed frames */
    
   A[0][0] = p->r.quat[0]*p->r.quat[0] + p->r.quat[1]*p->r.quat[1] -
             p->r.quat[2]*p->r.quat[2] - p->r.quat[3]*p->r.quat[3] ;
             
   A[1][1] = p->r.quat[0]*p->r.quat[0] - p->r.quat[1]*p->r.quat[1] +
             p->r.quat[2]*p->r.quat[2] - p->r.quat[3]*p->r.quat[3] ;
                          
   A[2][2] = p->r.quat[0]*p->r.quat[0] - p->r.quat[1]*p->r.quat[1] -
             p->r.quat[2]*p->r.quat[2] + p->r.quat[3]*p->r.quat[3] ;
	     	     
   A[0][1] = 2*(p->r.quat[1]*p->r.quat[2] + p->r.quat[0]*p->r.quat[3]);
	     
   A[0][2] = 2*(p->r.quat[1]*p->r.quat[3] - p->r.quat[0]*p->r.quat[2]);
        
   A[1][0] = 2*(p->r.quat[1]*p->r.quat[2] - p->r.quat[0]*p->r.quat[3]);
	     
      
   A[1][2] = 2*(p->r.quat[2]*p->r.quat[3] + p->r.quat[0]*p->r.quat[1]);
	     
   A[2][0] = 2*(p->r.quat[1]*p->r.quat[3] + p->r.quat[0]*p->r.quat[2]);
	     
   A[2][1] = 2*(p->r.quat[2]*p->r.quat[3] - p->r.quat[0]*p->r.quat[1]);
}

void define_Qd(Particle *p)
{
/** calculate the first derivative of the quaternion of a given particle */

  Qd[0] = 0.5 * ( -p->r.quat[1] * p->omega[0] -
                   p->r.quat[2] * p->omega[1] -
		   p->r.quat[3] * p->omega[2] );
		   
  Qd[1] = 0.5 * (  p->r.quat[0] * p->omega[0] -
                   p->r.quat[3] * p->omega[1] +
		   p->r.quat[2] * p->omega[2] );
		   
  Qd[2] = 0.5 * (  p->r.quat[3] * p->omega[0] +
                   p->r.quat[0] * p->omega[1] -
		   p->r.quat[1] * p->omega[2] );
		   
  Qd[3] = 0.5 * ( -p->r.quat[2] * p->omega[0] +
                   p->r.quat[1] * p->omega[1] +
		   p->r.quat[0] * p->omega[2] );
}

void define_Qdd(Particle *p)
{
/** calculate the second derivative of the quaternion of a given particle
    as well as Wd vector wich is the angular acceleration of this particle */
    
  S1 = Qd[0]*Qd[0] + Qd[1]*Qd[1] + Qd[2]*Qd[2] + Qd[3]*Qd[3] ;
   
  Wd[0] =  (p->torque[0] + p->omega[1]*p->omega[2]*(I[1]-I[2]))/I[0];  
  Wd[1] =  (p->torque[1] + p->omega[2]*p->omega[0]*(I[2]-I[0]))/I[1]; 
  Wd[2] =  (p->torque[2] + p->omega[0]*p->omega[1]*(I[0]-I[1]))/I[2];
        
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
}

/** propagate angular velocities and quaternions */
void propagate_omega_quat() 
{
  Particle *p;
  int m,n,o,i, np;
double lambda, dt2, S2, S3;

 dt2 = time_step*0.5;

  INTEG_TRACE(fprintf(stderr,"%d: propagate_omega_quat:\n",this_node));

  INNER_CELLS_LOOP(m, n, o) {
    p  = CELL_PTR(m, n, o)->pList.part;
    np = CELL_PTR(m, n, o)->pList.n;
    for(i = 0; i < np; i++) {
#ifdef EXTERNAL_FORCES
      if(p[i].ext_flag != PARTICLE_FIXED) 
#endif
	{ define_Qd(&p[i]); 
          define_Qdd(&p[i]);
	  
S2 = Qd[0]*Qdd[0]  + Qd[1]*Qdd[1]  + Qd[2]*Qdd[2]  + Qd[3]*Qdd[3];
S3 = Qdd[0]*Qdd[0] + Qdd[1]*Qdd[1] + Qdd[2]*Qdd[2] + Qdd[3]*Qdd[3];

lambda = ( 1 - S1*time_step*dt2 - sqrt(1 - S1*time_step*time_step - (S2 + 0.25*time_step*(S3-S1*S1))*time_step*time_step*time_step) );
	  
	p[i].omega[0]+= dt2*Wd[0];
	p[i].omega[1]+= dt2*Wd[1];
	p[i].omega[2]+= dt2*Wd[2];
	
      ONEPART_TRACE(if(p[i].r.identity==check_id) fprintf(stderr,"%d: OPT: PV_1 v_new = (%.3e,%.3e,%.3e)\n",this_node,p[i].v[0],p[i].v[1],p[i].v[2]));
	  
	  p[i].r.quat[0]+=  time_step*(Qd[0] + dt2*Qdd[0]) - lambda*p[i].r.quat[0]; 
	  p[i].r.quat[1]+=  time_step*(Qd[1] + dt2*Qdd[1]) - lambda*p[i].r.quat[1]; 
	  p[i].r.quat[2]+=  time_step*(Qd[2] + dt2*Qdd[2]) - lambda*p[i].r.quat[2]; 
	  p[i].r.quat[3]+=  time_step*(Qd[3] + dt2*Qdd[3]) - lambda*p[i].r.quat[3];  
	  
	}
	
      ONEPART_TRACE(if(p[i].r.identity==check_id) fprintf(stderr,"%d: OPT: PPOS p = (%.3f,%.3f,%.3f)\n",this_node,p[i].r.p[0],p[i].r.p[1],p[i].r.p[2]));
    }
  }
}

/** convert the toques to the body-fixed frames and propagate angular velocities */
void convert_torqes_propagate_omega()
{
  Particle *p;
  int m,n,o,i, np, times;
  double dt2, tx, ty, tz;
  
   dt2 = time_step*0.5;
  INTEG_TRACE(fprintf(stderr,"%d: convert_torqes_propagate_omega:\n",this_node));
  INNER_CELLS_LOOP(m, n, o) {
    p  = CELL_PTR(m, n, o)->pList.part;
    np = CELL_PTR(m, n, o)->pList.n;
    for(i = 0; i < np; i++) {define_rotation_matrix(&p[i]);
           
tx = A[0][0]*p[i].torque[0] + A[0][1]*p[i].torque[1] + A[0][2]*p[i].torque[2]; 
ty = A[1][0]*p[i].torque[0] + A[1][1]*p[i].torque[1] + A[1][2]*p[i].torque[2];
tz = A[2][0]*p[i].torque[0] + A[2][1]*p[i].torque[1] + A[2][2]*p[i].torque[2];

        friction_thermo_rotation(&p[i]);

	p[i].torque[0]+= tx;
	p[i].torque[1]+= ty;
	p[i].torque[2]+= tz;

      ONEPART_TRACE(if(p[i].r.identity==check_id) fprintf(stderr,"%d: OPT: SCAL f = (%.3e,%.3e,%.3e) v_old = (%.3e,%.3e,%.3e)\n",this_node,p[i].f[0],p[i].f[1],p[i].f[2],p[i].v[0],p[i].v[1],p[i].v[2]));

#ifdef EXTERNAL_FORCES
      if(p[i].ext_flag != PARTICLE_FIXED) 
#endif
	{ p[i].omega[0]+= dt2*p[i].torque[0]/I[0];
	  p[i].omega[1]+= dt2*p[i].torque[1]/I[1];
	  p[i].omega[2]+= dt2*p[i].torque[2]/I[2];
	  
/* if the tensor of inertia is isotrpic, the following refinement is not needed.
   Otherwise repeat this loop 2-3 times depending on the required accuracy */
for(times=0;times<=0;times++) { 
		   		 
  Wd[0] = (p[i].omega[1]*p[i].omega[2]*(I[1]-I[2]))/I[0];  
  Wd[1] = (p[i].omega[2]*p[i].omega[0]*(I[2]-I[0]))/I[1]; 
  Wd[2] = (p[i].omega[0]*p[i].omega[1]*(I[0]-I[1]))/I[2];
 
  p[i].omega[0]+= dt2*Wd[0];
  p[i].omega[1]+= dt2*Wd[1];
  p[i].omega[2]+= dt2*Wd[2];
		 }
	}

      ONEPART_TRACE(if(p[i].r.identity==check_id) fprintf(stderr,"%d: OPT: PV_2 v_new = (%.3e,%.3e,%.3e)\n",this_node,p[i].v[0],p[i].v[1],p[i].v[2]));
      
    }
  }
}

/** convert the toques to the body-fixed frames before the integration loop */
void convert_initial_torques()
{
  Particle *p;
  int m,n,o,i, np;
  double tx, ty, tz;
  
  INTEG_TRACE(fprintf(stderr,"%d: convert_initial_torqes:\n",this_node));
  INNER_CELLS_LOOP(m, n, o) {
    p  = CELL_PTR(m, n, o)->pList.part;
    np = CELL_PTR(m, n, o)->pList.n;
    for(i = 0; i < np; i++) {define_rotation_matrix(&p[i]);
           
tx = A[0][0]*p[i].torque[0] + A[0][1]*p[i].torque[1] + A[0][2]*p[i].torque[2]; 
ty = A[1][0]*p[i].torque[0] + A[1][1]*p[i].torque[1] + A[1][2]*p[i].torque[2];
tz = A[2][0]*p[i].torque[0] + A[2][1]*p[i].torque[1] + A[2][2]*p[i].torque[2];

	p[i].torque[0] = tx;
	p[i].torque[1] = ty;
	p[i].torque[2] = tz;

      ONEPART_TRACE(if(p[i].r.identity==check_id) fprintf(stderr,"%d: OPT: SCAL f = (%.3e,%.3e,%.3e) v_old = (%.3e,%.3e,%.3e)\n",this_node,p[i].f[0],p[i].f[1],p[i].f[2],p[i].v[0],p[i].v[1],p[i].v[2]));      
    }
  }  
}
#endif
