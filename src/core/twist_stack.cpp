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

#include "twist_stack.hpp"
#include "communication.hpp"

#ifdef TWIST_STACK

// extrema for cos(theta), used for the force calculations that involve angles
#define COS_MAX (0.99999999)
#define COS_MIN (-0.99999999)

#define SQR(x) ((x)*(x))

/* @TODO: Remove code duplication, use vector class */

/* calcualte dot product of two 3D vectors */
inline double dot(double *x1, double *x2)
{
  return x1[0]*x2[0] + x1[1]*x2[1] + x1[2]*x2[2];
}

/* calcualte cross product of two 3D vectors, x1 and x2. Store answer in x3 */
inline void cross(double *x1, double *x2, double *x3)
{
  x3[0] = x1[1]*x2[2] - x1[2]*x2[1];
  x3[1] = x1[2]*x2[0] - x1[0]*x2[2];
  x3[2] = x1[0]*x2[1] - x1[1]*x2[0];
}

inline double norm2(double *x) {
  return SQR(x[0]) + SQR(x[1]) + SQR(x[2]);
}

inline double norm(double *x) {
  return sqrt(norm2(x));
}

inline void normalize(double *x) {
  const double n = norm(x);
  if(n > 0.0) {
    x[0] /= n;
    x[1] /= n;
    x[2] /= n;
  }
}

inline double cos_angle(double *x1, double *x2) {
  return dot(x1, x2) / (norm(x1) * norm(x2) );
}

inline double angle(double *x1, double *x2) {
  return acos(cos_angle(x1, x2));
}

#ifdef TWIST_STACK_DEBUG
#define EX(A) A[0], A[1], A[2]

#define PV(A) printf(#A " = (%lf %lf %lf)\n",EX(A))
#define PS(A) printf(#A " = %lf\n", (double) A)
#else
#define PV(A) 
#define PS(A)
#endif /* TWIST_STACK_DEBUG */

int calc_twist_stack_energy(Particle *si1, Particle *bi1, Particle *bi2, Particle *si2,
				      Particle *sj1, Particle *bj1, Particle *bj2, Particle *sj2,
				      Bonded_ia_parameters *iaparams, double *_energy) {

  /* Base-Base and Sugar-Sugar vectors */
  double rcci[3], rccj[3];
  /* Sugar-Base Vectors */
  double rcb1[3], rcb2[3];
  double rcb1j[3], rcb2j[3];
  /* Mean basepair distance */
  double r;
  /* Base normals */
  double n1[3], n2[3], n1_l, n2_l;
  double n1j[3], n2j[3], n1j_l, n2j_l;

  double vec1[3];
  double vec2[3];
  double vec3[3];
  double vec4[3];

  double ani[3], anj[3];
  double ani_l, anj_l;
   
  get_mi_vector(vec1, sj1->r.p, si1->r.p);
  get_mi_vector(vec2, bj1->r.p, si1->r.p);
  get_mi_vector(vec3, sj2->r.p, si2->r.p);
  get_mi_vector(vec4, bj2->r.p, si2->r.p);

  get_mi_vector(rcci, si2->r.p, si1->r.p);
  get_mi_vector(rcb1, bi1->r.p, si1->r.p);
  get_mi_vector(rcb2, bi2->r.p, si2->r.p);
  
  get_mi_vector(rccj, sj2->r.p, sj1->r.p);
  get_mi_vector(rcb1j, bj1->r.p, sj1->r.p);
  get_mi_vector(rcb2j, bj2->r.p, sj2->r.p);

  cross(rcci, rcb1, n1);
  cross(rcci, rcb2, n2);

  cross(rccj, rcb1j, n1j);
  cross(rccj, rcb2j, n2j);

  n1_l = norm(n1);
  n2_l = norm(n2);
  n1j_l = norm(n1j);
  n2j_l = norm(n2j);
 
  n1_l = ( n1_l == 0 ) ? 1 : n1_l;
  n2_l = ( n2_l == 0 ) ? 1 : n2_l;

  n1j_l = ( n1j_l == 0 ) ? 1 : n1j_l;
  n2j_l = ( n2j_l == 0 ) ? 1 : n2j_l;

  for(int i = 0; i < 3; i++) {    
    n1[i] /= n1_l;
    n2[i] /= n2_l;
    n1j[i] /= n1j_l;
    n2j[i] /= n2j_l;
    ani[i] = n1[i] + n2[i];
    anj[i] = n1j[i] + n2j[i];
  }

  ani_l = norm(ani);
  anj_l = norm(anj);

  /* Prevent division by 0 in degenerate case */
  ani_l = (ani_l == 0) ? 1 : ani_l;
  anj_l = (anj_l == 0) ? 1 : anj_l;

  for(int i = 0; i < 3; i++) {    
    ani[i] /= ani_l;
    anj[i] /= anj_l;
  }

  r = 0.25*(dot(vec1, n1) + dot(vec2,n1) + dot(vec3,n2) + dot(vec4,n2));

  double pot_stack;

  const double ir = 1. / r;
  const double ir2 =SQR(ir);        
  const double ir5 = ir2*ir2*ir;
  const double ir6 = ir5*ir;

  const double rm = iaparams->p.twist_stack.rm; 
  const double epsilon = iaparams->p.twist_stack.epsilon;
  const double *a = iaparams->p.twist_stack.a;
  const double *b = iaparams->p.twist_stack.b;

  const double rm2 = rm*rm;
  const double rm5 = rm2*rm2*rm;
  const double rm6 = rm*rm5;
  const double eps5rm6 = 5.*epsilon*rm6;
  const double eps6rm5 = 6.*epsilon*rm5;

  pot_stack = eps5rm6*ir6 - eps6rm5*ir5;  

  /* Parallel projection of rccj */
  const double rccj_parallel = dot(rccj, n1);

  /* Projection of rccj to plane of bp i */
  double rccj_p[3];

  rccj_p[0] = rccj[0] - rccj_parallel*n1[0];
  rccj_p[1] = rccj[1] - rccj_parallel*n1[1];
  rccj_p[2] = rccj[2] - rccj_parallel*n1[2];
  
  const double rcci_l = norm(rcci);
  const double rccj_l = norm(rccj);

  const double rccj_p_l2 = SQR(rccj_l) - SQR(rccj_parallel);
  const double rccj_p_l = sqrt(rccj_p_l2);

  double cos1 = dot(rccj_p, rcci)/(rccj_p_l * rcci_l);
  
  cos1 = (cos1 > COS_MAX) ? COS_MAX : cos1;
  cos1 = (cos1 < COS_MIN) ? COS_MIN : cos1;

  cross(rcci, rccj_p, vec1);

  const double sin1 = (dot(vec1, n1) < 0.) ? -sqrt(1.-SQR(cos1)) : sqrt(1.-SQR(cos1));

  // Evaluation of cos(n*theta) by Chebyshev polynomials
  const double cos2 = 2.*cos1*cos1 - 1.;
  const double cos3 = 2.*cos2*cos1 - cos1;
  const double cos4 = 2.*cos3*cos1 - cos2;
  const double cos5 = 2.*cos4*cos1 - cos3;
  const double cos6 = 2.*cos5*cos1 - cos4;
  const double cos7 = 2.*cos6*cos1 - cos5;

  // Evaluation of sin(n*theta) by Chebyshev polynomials
  const double sin2 = 2.*sin1*cos1;
  const double sin3 = 2.*sin2*cos1 - sin1;
  const double sin4 = 2.*sin3*cos1 - sin2;
  const double sin5 = 2.*sin4*cos1 - sin3;
  const double sin6 = 2.*sin5*cos1 - sin4;
  const double sin7 = 2.*sin6*cos1 - sin5;

  const double pot_twist_ref = iaparams->p.twist_stack.ref_pot;

  const double pot_twist = a[0] + a[1]*cos1 + a[2]*cos2 + a[3]*cos3 + a[4]*cos4
    +a[5]*cos5 + a[6]*cos6 + a[7]*cos7 + b[0]*sin1 + b[1]*sin2 + b[2]*sin3 + b[3]*sin4 + b[4]*sin5 + b[5]*sin6 + b[6] *sin7;


  cos1 = dot(ani, anj);
  double tau_tilt;

  if(cos1 < 0) {
    tau_tilt = 0.0;
  } else {
    tau_tilt = cos1*cos1;
  }

  double tau_twist = pot_twist/pot_twist_ref;

  *_energy += pot_stack*tau_twist*tau_tilt;

  return 0;
}

int calc_twist_stack_force(Particle *si1, Particle *bi1, Particle *bi2, Particle *si2,
				      Particle *sj1, Particle *bj1, Particle *bj2, Particle *sj2,
				      Bonded_ia_parameters *iaparams,
				      double f_si1[3], double f_bi1[3], double f_bi2[3], double f_si2[3],
				      double f_sj1[3], double f_bj1[3], double f_bj2[3], double f_sj2[3]) {

  /* Base-Base and Sugar-Sugar vectors */
  double rcci[3], rccj[3];
  /* Sugar-Base Vectors */
  double rcb1[3], rcb2[3];
  double rcb1j[3], rcb2j[3];
  /* Mean basepair distance */
  double r;
  /* Base normals */
  double n1[3], n2[3], n1_l, n2_l;
  double n1j[3], n2j[3];
  double f_tilt_si1[3], f_tilt_si2[3], f_tilt_sj1[3], f_tilt_sj2[3];
  double f_tilt_bi1[3], f_tilt_bi2[3], f_tilt_bj1[3], f_tilt_bj2[3];
  double f_twist_si1[3], f_twist_si2[3], f_twist_sj1[3], f_twist_sj2[3];
  double f_twist_bi1[3];
  double f_stack_si1[3], f_stack_si2[3], f_stack_sj1[3], f_stack_sj2[3];
  double f_stack_bi1[3], f_stack_bi2[3], f_stack_bj1[3], f_stack_bj2[3];

  double vec1[3], u1[3], dot01, dot11;
  double vec2[3], u2[3], dot02, dot12;
  double vec3[3], u3[3], dot03, dot13;
  double vec4[3], u4[3], dot04, dot14;

  double ani[3], anj[3];
  double ani_l, anj_l;
   
  get_mi_vector(vec1, sj1->r.p, si1->r.p);
  get_mi_vector(vec2, bj1->r.p, si1->r.p);
  get_mi_vector(vec3, sj2->r.p, si2->r.p);
  get_mi_vector(vec4, bj2->r.p, si2->r.p);

  get_mi_vector(rcci, si2->r.p, si1->r.p);
  get_mi_vector(rcb1, bi1->r.p, si1->r.p);
  get_mi_vector(rcb2, bi2->r.p, si2->r.p);
  
  get_mi_vector(rccj, sj2->r.p, sj1->r.p);
  get_mi_vector(rcb1j, bj1->r.p, sj1->r.p);
  get_mi_vector(rcb2j, bj2->r.p, sj2->r.p);

  cross(rcci, rcb1, n1);
  cross(rcci, rcb2, n2);

  cross(rccj, rcb1j, n1j);
  cross(rccj, rcb2j, n2j);

  n1_l = norm(n1);
  n2_l = norm(n2);
 
  n1_l = ( n1_l == 0 ) ? 1 : n1_l;
  n2_l = ( n2_l == 0 ) ? 1 : n2_l;

  for(int i = 0; i < 3; i++) {    
    ani[i] = n1[i] + n2[i];
    anj[i] = n1j[i] + n2j[i];
    n1[i] /= n1_l;
    n2[i] /= n2_l;
  }

  ani_l = norm(ani);
  anj_l = norm(anj);

  /* Prevent division by 0 in degenerate case */
  ani_l = (ani_l == 0) ? 1 : ani_l;
  anj_l = (anj_l == 0) ? 1 : anj_l;

  for(int i = 0; i < 3; i++) {    
    ani[i] /= ani_l;
    anj[i] /= anj_l;
  }

  r = 0.25*(dot(vec1, n1) + dot(vec2,n1) + dot(vec3,n2) + dot(vec4,n2));

  double f_r;
  double pot_stack;

  const double ir = 1. / r;
  const double ir2 =SQR(ir);        
  const double ir5 = ir2*ir2*ir;
  const double ir6 = ir5*ir;

  const double rm = iaparams->p.twist_stack.rm; 
  const double epsilon = iaparams->p.twist_stack.epsilon;
  const double *a = iaparams->p.twist_stack.a;
  const double *b = iaparams->p.twist_stack.b;

  const double rm2 = rm*rm;
  const double rm5 = rm2*rm2*rm;
  const double rm6 = rm*rm5;
  const double eps5rm6 = 5.*epsilon*rm6;
  const double eps6rm5 = 6.*epsilon*rm5;
  const double eps30rm6 = 6.*eps5rm6;
  const double eps30rm5 = 5.*eps6rm5;

  pot_stack = eps5rm6*ir6 - eps6rm5*ir5;  
  f_r = 0.25*ir*(eps30rm6*ir6 - eps30rm5*ir5);

  cross(n1, vec1, u1);
  cross(n1, vec2, u2);
  cross(n2, vec3, u3);
  cross(n2, vec4, u4);
  dot01 = dot(u1, rcci)/n1_l;
  dot02 = dot(u2, rcci)/n1_l;
  dot03 = dot(u3, rcci)/n2_l;
  dot04 = dot(u4, rcci)/n2_l;
  dot11 = dot(u1, rcb1)/n1_l;
  dot12 = dot(u2, rcb1)/n1_l;
  dot13 = dot(u3, rcb2)/n2_l;
  dot14 = dot(u4, rcb2)/n2_l;

  double mag1, mag2;

  PV(n1);
  PV(n2);
  PS(f_r);

  for(int k = 0; k < 3; k++) {
    mag1 = f_r*n1[k];
    mag2 = f_r*n2[k];
    PS(k);
    PS(mag1);
    PS(mag2);
    PS(n1[k]);
    PS(n2[k]);
    f_stack_bi1[k] = (dot01+dot02)*mag1;
    f_stack_bi2[k] = (dot03+dot04)*mag2;
    f_stack_si1[k] = -(2.+dot01-dot11+dot02-dot12)*mag1 + (dot13+dot14)*mag2;
    f_stack_si2[k] = -(2.+dot03+dot13+dot04+dot14)*mag2 - (dot11+dot12)*mag1;
    f_stack_bj1[k] = mag1;
    f_stack_bj2[k] = mag2;
    f_stack_sj1[k] = mag1;
    f_stack_sj2[k] = mag2;    
  }

  /* Parallel projection of rccj */
  const double rccj_parallel = dot(rccj, n1);

  /* Projection of rccj to plane of bp i */
  double rccj_p[3];

  rccj_p[0] = rccj[0] - rccj_parallel*n1[0];
  rccj_p[1] = rccj[1] - rccj_parallel*n1[1];
  rccj_p[2] = rccj[2] - rccj_parallel*n1[2];
  
  const double rcci_l = norm(rcci);
  const double rccj_l = norm(rccj);

  const double rccj_p_l2 = SQR(rccj_l) - SQR(rccj_parallel);
  const double rccj_p_l = sqrt(rccj_p_l2);

  double cos1 = dot(rccj_p, rcci)/(rccj_p_l * rcci_l);

  // static FILE *fd = fopen("cos1.dat", "w");

  // fprintf(fd, "%e\n", cos1);
  // fflush(fd);
  
  cos1 = (cos1 > COS_MAX) ? COS_MAX : cos1;
  cos1 = (cos1 < COS_MIN) ? COS_MIN : cos1;

#ifdef TWIST_STACK_DEBUG
  double vec1o[3];
  vec1o[0] = vec1[0];
  vec1o[1] = vec1[1];
  vec1o[2] = vec1[2];
#endif


  cross(rcci, rccj_p, vec1);

  const double sin1 = (dot(vec1, n1) < 0.) ? -sqrt(1.-SQR(cos1)) : sqrt(1.-SQR(cos1));

  // Evaluation of cos(n*theta) by Chebyshev polynomials
  const double cos2 = 2.*cos1*cos1 - 1.;
  const double cos3 = 2.*cos2*cos1 - cos1;
  const double cos4 = 2.*cos3*cos1 - cos2;
  const double cos5 = 2.*cos4*cos1 - cos3;
  const double cos6 = 2.*cos5*cos1 - cos4;
  const double cos7 = 2.*cos6*cos1 - cos5;

  // Evaluation of sin(n*theta) by Chebyshev polynomials
  const double sin2 = 2.*sin1*cos1;
  const double sin3 = 2.*sin2*cos1 - sin1;
  const double sin4 = 2.*sin3*cos1 - sin2;
  const double sin5 = 2.*sin4*cos1 - sin3;
  const double sin6 = 2.*sin5*cos1 - sin4;
  const double sin7 = 2.*sin6*cos1 - sin5;

  const double pot_twist_ref = iaparams->p.twist_stack.ref_pot;

  const double pot_twist = a[0] + a[1]*cos1 + a[2]*cos2 + a[3]*cos3 + a[4]*cos4
    +a[5]*cos5 + a[6]*cos6 + a[7]*cos7 + b[0]*sin1 + b[1]*sin2 + b[2]*sin3 + b[3]*sin4 + b[4]*sin5 + b[5]*sin6 + b[6] *sin7;

  double fmag = a[1]*sin1 - b[0] * cos1 + 2.* (a[2]*sin2 - b[1] *cos2) + 3.*(a[3]*sin3 - b[2]*cos3) + 4. * (a[4]*sin4 - b[3]*cos4) + 5.*(a[5]*sin5 - b[4]*cos5) + 6.*(a[6]*sin6 - b[5]*cos6) + 7. * (a[7]*sin7 - b[6]*cos7);

  fmag = -fmag/sin1;

  cross(n1, rccj, u1);
  dot01 = dot(u1, rcci)/n1_l;
  dot11 = dot(u1, rcb1)/n1_l;
  double factor1 = fmag/(rcci_l*rccj_p_l);
  double factor2 = fmag*cos1/SQR(rcci_l);
  double factor3 = fmag*cos1/rccj_p_l2;
  const double factor4 = factor3*rccj_parallel;

  double mag0; 
  double mag3; 
  double mag4; 

  for(int k = 0; k < 3; k++) {
    mag0 = factor1*rcci[k];
    mag1 = factor1*rccj[k];
    mag2 = factor2*rcci[k];
    mag3 = factor3*rccj[k];
    mag4 = factor4*n1[k];

    f_twist_bi1[k] = dot01*mag4;
    f_twist_si1[k] = -mag1 + mag2 + (dot11-dot01)*mag4;
    f_twist_si2[k] = mag1 - mag2 - dot11*mag4;
    f_twist_sj1[k] = -mag0 + mag3 - mag4;
    f_twist_sj2[k] = -f_twist_sj1[k];
  }

  cos1 = dot(ani, anj);
  double tau_tilt;
  double f_tilt;

  double ui[3], uj[3], veci[3], vecj[3];

  if(cos1 < 0) {
    tau_tilt = 0.0;
    for(int i = 0; i < 3; i++) {
      f_tilt_si1[i] = f_tilt_si2[i] = f_tilt_sj1[i] = f_tilt_sj2[i] = 0;
      f_tilt_bi1[i] = f_tilt_bi2[i] = f_tilt_bj1[i] = f_tilt_bj2[i] = 0;
    }
  } else {
    tau_tilt = cos1*cos1;
    f_tilt = -2.*pot_twist_ref*cos1;

    get_mi_vector(veci, bi1->r.p, si2->r.p);
    get_mi_vector(vecj, bj1->r.p, sj2->r.p);

    for(int i = 0; i < 3; i++) {
      ui[i] = f_tilt*(anj[i] - cos1*ani[i])/ani_l;
      uj[i] = f_tilt*(ani[i] - cos1*anj[i])/anj_l;
      veci[i] += rcb2[i];
      vecj[i] += rcb2j[i];
    }

    cross(ui, rcci, f_tilt_bi1);
    cross(ui, veci, f_tilt_si1);
    cross(uj, rccj, f_tilt_bj1);
    cross(uj, vecj, f_tilt_sj1);
    
    for(int k = 0; k < 3; k++) {
      f_tilt_bi2[k] = f_tilt_bi1[k];
      f_tilt_si2[k] = -f_tilt_si1[k] - 2.*f_tilt_bi1[k];
      f_tilt_bj2[k] = f_tilt_bj1[k];
      f_tilt_sj2[k] = -f_tilt_sj1[k] - 2.*f_tilt_bj1[k];
    }
  }

  const double tau_stack = pot_stack/pot_twist_ref;
  const double tau_twist = pot_twist/pot_twist_ref;

  // printf("tau_stack %e, tau_twist %e, tau_tilt %e cos1 %e, \n",
  // 	 tau_stack, tau_twist, tau_tilt, cos1);
  // printf("pot_twist_ref %e\n", pot_twist_ref);

  factor1 = tau_twist*tau_tilt;
  factor2 = tau_stack*tau_tilt;
  factor3 = tau_stack*tau_twist;


#ifdef TWIST_STACK_DEBUG
  fflush(NULL);
  int big_force = 0;
#endif

  for(int k = 0; k < 3; k++) {
    f_bi1[k] = factor1*f_stack_bi1[k] + factor2*f_twist_bi1[k] + factor3*f_tilt_bi1[k];
    f_bi2[k] = factor1*f_stack_bi2[k]                          + factor3*f_tilt_bi2[k];
    f_si1[k] = factor1*f_stack_si1[k] + factor2*f_twist_si1[k] + factor3*f_tilt_si1[k];
    f_si2[k] = factor1*f_stack_si2[k] + factor2*f_twist_si2[k] + factor3*f_tilt_si2[k];

    f_bj1[k] = factor1*f_stack_bj1[k]                          + factor3*f_tilt_bj1[k];
    f_bj2[k] = factor1*f_stack_bj2[k]                          + factor3*f_tilt_bj2[k];
    f_sj1[k] = factor1*f_stack_sj1[k] + factor2*f_twist_sj1[k] + factor3*f_tilt_sj1[k];
    f_sj2[k] = factor1*f_stack_sj2[k] + factor2*f_twist_sj2[k] + factor3*f_tilt_sj2[k];

#ifdef TWIST_STACK_DEBUG
    if((f_bi1[k] >= 100.) || (f_bi2[k] >= 100.) || (f_si1[k] >= 100.) || (f_si2[k] >= 100.) || (f_bj1[k] >= 100.) || (f_bj2[k] >= 100.) || (f_sj1[k] >= 100.) || (f_sj2[k] >= 100.)) {
      big_force = 1;
    }
#endif

 
  }	  

#ifdef TWIST_STACK_DEBUG
  if(big_force) {
    puts("Big Force Twist/Stack.");
    PS(si1->p.identity);

    PS(sim_time);
    PS(time_step);
    PS(sim_time/time_step);

    PV(si1->r.p);
    PV(si2->r.p);
    PV(bi1->r.p);
    PV(bi2->r.p);

    PV(sj1->r.p);
    PV(sj2->r.p);
    PV(bj1->r.p);
    PV(bj2->r.p);

    PS(rcci_l);
    PS(rccj_l);

    PV(vec1o);
    PV(vec2);
    PV(vec3);
    PV(vec4);

    PV(n1);
    PV(n2);
    PV(n1j);
    PV(n2j);
    PV(ani);
    PV(anj);

    PS(dot(vec1o,n1));
    PS(dot(vec2,n1));
    PS(dot(vec3,n2));
    PS(dot(vec4,n2));

    PS(r);
    PS(cos1);
    PS(epsilon);

    PS(tau_twist);
    PS(tau_tilt);
    PS(tau_stack);

    PS(factor1);
    PS(factor2);
    PS(factor3);

    PS(fmag);
    PS(rccj_p_l2);
    PS(pot_stack);
    PS(pot_twist_ref);

    PV(factor1*f_stack_si1);
    PV(factor1*f_stack_si2);
    PV(factor1*f_stack_bi1);
    PV(factor1*f_stack_bi2);
    PV(factor1*f_stack_sj1);
    PV(factor1*f_stack_sj2);
    PV(factor1*f_stack_bj1);
    PV(factor1*f_stack_bj2);

    PV(factor2*f_twist_si1);
    PV(factor2*f_twist_si2);
    PV(factor2*f_twist_bi1);
    PV(factor2*f_twist_sj1);
    PV(factor2*f_twist_sj2);

    PV(factor3*f_tilt_si1);
    PV(factor3*f_tilt_si2);
    PV(factor3*f_tilt_bi1);
    PV(factor3*f_tilt_bi2);
    PV(factor3*f_tilt_sj1);
    PV(factor3*f_tilt_sj2);
    PV(factor3*f_tilt_bj1);
    PV(factor3*f_tilt_bj2);


    PV(f_bi1);
    PV(f_bi2);
    PV(f_si1);
    PV(f_si2);

    PV(f_bj1);
    PV(f_bj2);
    PV(f_sj1);
    PV(f_sj2);

  }
#endif

  return 0;
}

int twist_stack_set_params(int bond_type, DoubleList *params) {
  if(bond_type < 0)
    return ES_ERROR;

  make_bond_type_exist(bond_type);

  bonded_ia_params[bond_type].p.twist_stack.rm = params->e[0];
  bonded_ia_params[bond_type].p.twist_stack.epsilon = params->e[1];

  for(int i = 0; i < 8; i++)
    bonded_ia_params[bond_type].p.twist_stack.a[i] = params->e[2+i];

  for(int i = 0; i < 7; i++)
    bonded_ia_params[bond_type].p.twist_stack.b[i] = params->e[10+i];
  
  const double dt = PI*36./180.;
  const double *a = bonded_ia_params[bond_type].p.twist_stack.a;
  const double *b = bonded_ia_params[bond_type].p.twist_stack.b;

  bonded_ia_params[bond_type].p.twist_stack.ref_pot = a[0] + a[1] * cos(dt) + a[2] * cos(2*dt) + a[3] * cos(3*dt) + a[4] * cos(4*dt) + a[5] * cos(5*dt) + a[6] * cos(6*dt) + a[7] * cos(7*dt) + b[0] * sin(dt) + b[1] * sin(2*dt) + b[2] * sin(3*dt) + b[3] * sin(4*dt) + b[5] * sin(5*dt) + b[6] * sin(6*dt);
  bonded_ia_params[bond_type].type = BONDED_IA_CG_DNA_STACKING;
  bonded_ia_params[bond_type].num = 7;

  mpi_bcast_ia_params(bond_type, -1);

  return ES_OK;
}

#endif /* TWIST_STACK */
