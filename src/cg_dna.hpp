/*
  Copyright (C) 2010,2011,2012,2013,2014 The ESPResSo project
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

#ifndef CG_DNA_HPP
#define CG_DNA_HPP

#include "utils.hpp"
#include "interaction_data.hpp"
#include "particle_data.hpp"
#include "grid.hpp"

#ifdef CG_DNA

// extrema for cos(theta), used for the force calculations that involve angles
#define COS_MAX (0.99999999)
#define COS_MIN (-0.99999999)

int cg_dna_basepair_set_params(int bond_type, DoubleList *params);

int cg_dna_stacking_set_params(int bond_type, DoubleList *params);

#define SQR(x) ((x)*(x))

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

inline double cos_angle(double *x1, double *x2) {
  return dot(x1, x2) / (norm(x1) * norm(x2) );
}

inline double angle(double *x1, double *x2) {
  return acos(cos_angle(x1, x2));
}

/* get_mi_vector(c, a, b ) (a,b) -> a -x */

#define EX(A) A[0], A[1], A[2]
#define PV(A) printf(#A " = (%lf %lf %lf)\n",EX(A))
#define PS(A) printf(#A " = %lf\n", (double) A)

inline int calc_cg_dna_stacking_force(Particle *si1, Particle *bi1, Particle *si2, Particle *bi2, 
				      Particle *sj1, Particle *bj1, Particle *sj2, Particle *bj2, 
				      Bonded_ia_parameters *iaparams,
				      double force1[3], double force2[3], double force3[3], double force4[3],
				      double force5to8[12]) {

  /* Base-Base and Sugar-Sugar vectors */
  double rcci[3], rccj[3];
  /* Sugar-Base Vectors */
  double rcb1[3], rcb2[3], rcb1_l, rcb2_l;
  double rcb1j[3], rcb2j[3];
  /* Mean basepair distance */
  double r;
  /* Base normals */
  double n1[3], n2[3], n1_l, n2_l;
  double n1j[3], n2j[3], n1j_l, n2j_l;
  double f_tilt_si1[3], f_tilt_si2[3], f_tilt_sj1[3], f_tilt_sj2[3];
  double f_tilt_bi1[3], f_tilt_bi2[3], f_tilt_bj1[3], f_tilt_bj2[3];
  double f_twist_si1[3], f_twist_si2[3], f_twist_sj1[3], f_twist_sj2[3];
  double f_twist_bi1[3], f_twist_bi2[3], f_twist_bj1[3], f_twist_bj2[3];
  double f_stack_si1[3], f_stack_si2[3], f_stack_sj1[3], f_stack_sj2[3];
  double f_stack_bi1[3], f_stack_bi2[3], f_stack_bj1[3], f_stack_bj2[3];

  double vec1[3], u1[3], dot01, dot11;
  double vec2[3], u2[3], dot02, dot12;
  double vec3[3], u3[3], dot03, dot13;
  double vec4[3], u4[3], dot04, dot14;

  double ani[3], anj[3];
  double ani_l, anj_l;

  // x[i][iRs][k]  = x[i][iSugar2][k] - x[i][iSugar1][k];  // sugar-sugar vector
  // x[i][iRb][k]  = x[i][iBase2][k]  - x[i][iBase1][k];   // base-base vector
  // x[i][iRb1][k] = x[i][iBase1][k]  - x[i][iSugar1][k];  // base-sugar vector 1
  // x[i][iRb2][k] = x[i][iBase2][k]  - x[i][iSugar2][k];  // base-sugar vector 2

  // get_mi_vector(c, a, b) -> c = a - b

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

  // PV(ani);
  // PV(anj); 
  // PV(vec1);
  // PV(vec2);
  // PV(vec3);
  // PV(vec4);
  // PV(n1);
  // PV(n2);

  r = 0.25*(dot(vec1, n1) + dot(vec2,n1) + dot(vec3,n2) + dot(vec4,n2));

  // PS(dot(vec1, n1));
  // PS(dot(vec2, n1));
  // PS(dot(vec3, n2));
  // PS(dot(vec4, n2));

  // PS((dot(vec1, n1) + dot(vec2,n1) + dot(vec3,n2) + dot(vec4,n2)));

  // PS(r);

  /* In the degenrate case take the average euclidian distance between the sugars. */
  if( r == 0) {
    r = 0.5 * (norm(vec1) +norm(vec3));
  }
  
  // PS(r);

  double f_r;
  double pot_stack;

  const double ir = 1. / r;
  const double ir2 =SQR(ir);        
  const double ir5 = ir2*ir2*ir;
  const double ir6 = ir5*ir;

  const double rm = iaparams->p.cg_dna_stacking.rm; 
  const double epsilon = iaparams->p.cg_dna_stacking.epsilon;
  const double *a = iaparams->p.cg_dna_stacking.a;
  const double *b = iaparams->p.cg_dna_stacking.b;

  const double rm2 = rm*rm;
  const double rm5 = rm2*rm2*rm;
  const double rm6 = rm*rm5;
  const double eps5rm6 = 5.*epsilon*rm6;
  const double eps6rm5 = 6.*epsilon*rm5;
  const double eps30rm6 = 6.*eps5rm6;
  const double eps30rm5 = 5.*eps6rm5;

  pot_stack = eps5rm6*ir6 - eps6rm5*ir5;  
  f_r = 0.25*ir*(eps30rm6*ir6 - eps30rm5*ir5);

  // PS(pot_stack);
  // PS(f_r);

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

  // PS(dot01); 
  // PS(dot02); 
  // PS(dot03); 
  // PS(dot04); 

  // PS(dot11); 
  // PS(dot12); 
  // PS(dot13); 
  // PS(dot14); 

  // PV(u1);
  // PV(u2);
  // PV(u3);
  // PV(u4);

  double mag1, mag2;

  for(int k = 0; k < 3; k++) {
    mag1 = f_r*n1[k];
    mag2 = f_r*n2[k];
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

  const double pot_twist = 0.5*a[0] + a[1]*cos1 + a[2]*cos2 + a[3]*cos3 + a[4]*cos4
    +a[5]*cos5 + a[6]*cos6 + a[7]*cos7 + b[0]*sin1 + b[1]*sin2 + b[2]*sin3 + b[3]*sin4 + b[4]*sin5 + b[5]*sin6 + b[6] *sin7;

  double fmag = a[1]*sin1 - b[0] * cos1 + 2.* (a[2]*sin2 - b[1] *cos2) + 3.*(a[3]*sin3 - b[2]*cos3) + 4. * (a[4]*sin4 - b[3]*cos4) + 5.*(a[5]*sin5 - b[4]*cos5) + 6.*(a[6]*sin6 - b[5]*cos6) + 7. * (a[7]*sin7 - b[6]*cos7);

  fmag = -fmag/sin1;
 
  // PS(fmag);

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
    f_tilt = -2.*epsilon*cos1;
    
    // PS(tau_tilt);
    // PS(cos1);

    for(int i = 0; i < 3; i++) {
      ui[i] = f_tilt*(anj[i] - cos1*ani[i])/ani_l;
      uj[i] = f_tilt*(ani[i] - cos1*anj[i])/anj_l;
      veci[i] = bi1->r.p[i] - si2->r.p[i] + rcb2[i];
      vecj[i] = bj1->r.p[i] - sj2->r.p[i] + rcb2j[i];
    }

    // PV(ui);
    // PV(rcci);

    cross(ui, rcci, f_tilt_bi1);
    cross(ui, veci, f_tilt_si1);
    cross(uj, rccj, f_tilt_bj1);
    cross(uj, vecj, f_tilt_sj1);
    
  // printf("f_tilt_bi1 (%e %e %e)\n", EX(f_tilt_bi1));

    for(int k = 0; k < 3; k++) {
      f_tilt_bi2[k] = f_tilt_bi1[k];
      f_tilt_sj2[k] = -f_tilt_si1[k] - 2.*f_tilt_bi1[k];
      f_tilt_bj2[k] = f_tilt_bj1[k];
      f_tilt_sj2[k] = -f_tilt_sj1[k] - 2.*f_tilt_bj1[k];
    }
  }

  double tau_stack = pot_stack/epsilon;
  double tau_twist = pot_twist/epsilon;

  factor1 = tau_twist*tau_tilt;
  factor2 = tau_stack*tau_tilt;
  factor3 = tau_stack*tau_twist;

  // PS(tau_twist);
  // PS(tau_tilt);
  // PS(tau_stack);

  // PS(factor1);
  // PS(factor2);
  // PS(factor3);
  // PS(factor4);

  // PV(factor1*f_stack_si1);
  // PV(factor1*f_stack_si2);
  // PV(factor1*f_stack_bi1);
  // PV(factor1*f_stack_bi2);
  // PV(factor1*f_stack_sj1);
  // PV(factor1*f_stack_sj2);
  // PV(factor1*f_stack_bj1);
  // PV(factor1*f_stack_bj2);

  // PV(factor2*f_twist_si1);
  // PV(factor2*f_twist_si2);
  // PV(factor2*f_twist_bi1);
  // PV(factor2*f_twist_bi2);
  // PV(factor2*f_twist_sj1);
  // PV(factor2*f_twist_sj2);
  // PV(factor2*f_twist_bj1);
  // PV(factor2*f_twist_bj2);

  // PV(factor3*f_tilt_si1);
  // PV(factor3*f_tilt_si2);
  // PV(factor3*f_tilt_bi1);
  // PV(factor3*f_tilt_bi2);
  // PV(factor3*f_tilt_sj1);
  // PV(factor3*f_tilt_sj2);
  // PV(factor3*f_tilt_bj1);
  // PV(factor3*f_tilt_bj2);


  for(int k = 0; k < 3; k++) {
    force2[k] = factor1*f_stack_bi1[k] + factor2*f_twist_bi1[k] + factor3*f_tilt_bi1[k];
    force4[k] = factor1*f_stack_bi2[k]                          + factor3*f_tilt_bi2[k];
    force1[k] = factor1*f_stack_si1[k] + factor2*f_twist_si1[k] + factor3*f_tilt_si1[k];
    force3[k] = factor1*f_stack_si2[k] + factor2*f_twist_si2[k] + factor3*f_tilt_si2[k];
    force5to8[3 + k] = factor1*f_stack_bj1[k] + factor2*f_twist_bj1[k] + factor3*f_tilt_bj1[k];
    force5to8[9 + k] = factor1*f_stack_bj2[k]                          + factor3*f_tilt_bj2[k];
    force5to8[0 + k] = factor1*f_stack_sj1[k] + factor2*f_twist_sj1[k] + factor3*f_tilt_sj1[k];
    force5to8[6 + k] = factor1*f_stack_sj2[k] + factor2*f_twist_sj2[k] + factor3*f_tilt_sj2[k];
  }	  

  return 0;
}

inline int calc_cg_dna_basepair_force(Particle *p1, Particle *p2, Particle *p3, Particle *p4, Bonded_ia_parameters *iaparams, double force1[3], double force2[3], double force3[3], double force4[3]) {

  /* Base-Base and Sugar-Sugar vectors */
  double rhb[3], rcc[3];
  /* Sugar-Base Vectors */
  double rcb1[3], rcb2[3], rcb1_l, rcb2_l;
  /* Normal vectors of the base pairs */
  double n1[3], n2[3];
  /* Dihedral of the base pair */
  double thetad;
  /* Base angles */
  double psi1, psi2;
  double dpsi1, dpsi2;
  /* gamma1 = cos(psi1) */
  /* length of rhb */
  double gamma1, gamma2;
  double rhb_l;
  /* magnitude of the _r contribution to the force */
  double f_r;
  double f_d;
  double f_f1, f_f2;
  double f_sb1, f_sb2;

  /* helper variables */
  double ra, c, tau_r, tau_d, tau_flip, tau_rd;
  const Bonded_ia_parameters params = *iaparams;
  const double E0 = params.p.cg_dna_basepair.E0;

  get_mi_vector(rhb, p4->r.p, p2->r.p);
  get_mi_vector(rcc, p3->r.p, p1->r.p);
  get_mi_vector(rcb1, p2->r.p, p1->r.p);
  get_mi_vector(rcb2, p4->r.p, p3->r.p);

  rcb1_l = norm(rcb1);
  rcb2_l = norm(rcb2);

  gamma1 = cos_angle(rcc, rcb1);
  gamma2 = cos_angle(rcc, rcb2);

  cross(rcc, rcb1, n1);
  cross(rcc, rcb2, n2);

  double n1_l2 = norm2(n1);
  double n2_l2 = norm2(n2);

  if((n1_l2 == 0) || (n2_l2 == 0)) {
    // puts("(n1_l2 == 0) || (n2_l2 == 0)");
    thetad = 0;
    n1_l2 = 1;
    n2_l2 = 1;
  } else {
    thetad = (gamma1 != 1. && gamma2 != 1.) ? angle(n1, n2) : 0;
  }

  // printf("cos(<(n1,n2)) = %lf, cos(<(rcc,rcb1)) = %lf\n", cos_angle(n1,n2), gamma1);
  // PRINT_VECTOR(n1);
  // PRINT_VECTOR(n2);

  if((gamma1 > COS_MAX) || (gamma1 < COS_MIN)) {
    gamma1 = (gamma1 > COS_MAX)? COS_MAX : COS_MIN;
  }
  if((gamma2 > COS_MAX) || (gamma2 < COS_MIN))
    gamma2 = (gamma2 > COS_MAX)? COS_MAX : COS_MIN;

  psi1 = gamma1 >= 1. ? 0. : (gamma1 <= -1. ? M_PI : acos(gamma1));
  psi2 = gamma2 >= 1. ? 0. : (gamma2 <= -1. ? M_PI : acos(gamma2));  

  rhb_l = norm(rhb);

  // printf("p1 %lf %lf %lf, p2 %lf %lf %lf, p3 %lf %lf %lf, p4 %lf %lf %lf\n",
  // 	 p1->r.p[0], p1->r.p[1], p1->r.p[2], p2->r.p[0], p2->r.p[1], p2->r.p[2], 
  // 	 p3->r.p[0], p3->r.p[1], p3->r.p[2], p4->r.p[0], p4->r.p[1], p4->r.p[2]);

  // PRINT_VECTOR(rhb);
  // PRINT_VECTOR(rcc);
  // PRINT_VECTOR(rcb1);
  // PRINT_VECTOR(rcb2);
  // PRINT_VECTOR(n1);
  // PRINT_VECTOR(n2);

  /* Sugar base interaction strand 1 */
  
  const double r0sb = params.p.cg_dna_basepair.r0sb;
  const double alphasb = params.p.cg_dna_basepair.alphasb;
  const double f2 = params.p.cg_dna_basepair.f2;
  const double f3 = params.p.cg_dna_basepair.f3;
  const double E0sb = params.p.cg_dna_basepair.E0sb;
  const double c0sb = (1. - 2.*f2)*E0sb*alphasb;
  const double c1sb = (f2-3.*f3)*E0sb*alphasb;
  const double c2sb = f3*E0sb*alphasb;

  double temp;

  temp = (rcb1_l - r0sb)*alphasb;
  f_sb1 = exp(-temp)*temp*(c0sb+c1sb*temp+c2sb*temp*temp)/rcb1_l;

  temp = (rcb2_l - r0sb)*alphasb;
  f_sb2 = exp(-temp)*temp*(c0sb+c1sb*temp+c2sb*temp*temp)/rcb2_l;

  /* Radial part */

  ra = (rhb_l - params.p.cg_dna_basepair.r0)*params.p.cg_dna_basepair.alpha;
  c = exp(-ra);
  tau_r = c *(1.+ra);
  f_r = E0*params.p.cg_dna_basepair.alpha*c*ra/rhb_l;

  /* Dihedral part */

  tau_d = exp(params.p.cg_dna_basepair.kd*(thetad - 1));
  f_d = -E0*params.p.cg_dna_basepair.kd*tau_d;  

  /* Flip part */

  dpsi1 = psi1 - params.p.cg_dna_basepair.psi10;
  dpsi2 = psi2 - params.p.cg_dna_basepair.psi20;

  double s1sqr = (SQR(params.p.cg_dna_basepair.sigma1));
  double s2sqr = (SQR(params.p.cg_dna_basepair.sigma2));

  // printf("gamma1 %lf gamma2 %lf\n", gamma1, gamma2);

  f_f1 = (gamma1 < 1.) ? -dpsi1 / sqrt(1. - SQR(gamma1)) * params.p.cg_dna_basepair.E0/s1sqr : 1.0;
  f_f2 = (gamma2 < 1.) ? -dpsi2 / sqrt(1. - SQR(gamma2)) * params.p.cg_dna_basepair.E0/s2sqr : 1.0;

  // printf("f_f1 %lf, f_f2 %lf\n", f_f1, f_f2);
  
  tau_rd = tau_r * tau_d;

  if(dpsi1 > 0 && dpsi2 > 0) {
    tau_flip = exp(-(SQR(dpsi1)/(2*s1sqr)+SQR(dpsi2)/(2*s2sqr)));
    f_f1 *= tau_flip * tau_rd;
    f_f2 *= tau_flip * tau_rd;
  } else if (dpsi1 > 0.) {
    tau_flip = exp(-(SQR(dpsi1)/(2*s1sqr)));
    f_f1 *= tau_flip * tau_rd;
  }
  else if (dpsi2 > 0.) {
    tau_flip = exp(-(SQR(dpsi2)/(2*s2sqr)));
    f_f2 *= tau_flip * tau_rd;
  } else {
    tau_flip = 1.;
  }

  /* Dihedral force */
  double vec[3];
  cross(n1, n2, vec);
  // puts("vec");
  // PRINT_VECTOR(vec);
  const double dot1 = f_d * dot(vec, rhb);
  const double dot2 = f_d * dot(vec, rcb1);
  const double dot3 = f_d * dot(vec, rcb2);
  const double dot4 = dot2 - dot1;
  const double dot5 = dot1 + dot3;

  const double rcc_l = norm(rcc);

  // printf("rcb1_l %lf, rcb2_l %lf, rcc_l %lf\n", rcb1_l, rcb2_l, rcc_l);
  // printf("f_f1 %lf, f_f2 %lf\n", f_f1, f_f2);

  const double factor1 = f_f1/(rcc_l * rcb1_l);
  const double factor2 = f_f1*gamma1/SQR(rcb1_l);
  const double factor3 = f_f1*gamma1/SQR(rcc_l);

  const double factor4 = f_f1/(rcc_l * rcb2_l);
  const double factor5 = f_f2*gamma2/SQR(rcb2_l);
  const double factor6 = f_f2*gamma2/SQR(rcc_l);

  // printf("factorx (%lf, %lf, %lf, %lf, %lf, %lf)\n",
  // 	 factor1,factor2,factor3,
  // 	 factor4,factor5,factor6);

  double fBase1, fSugar2, fBase2, fSugar1;
  double fr;

  for(int i = 0; i < 3; i++) {
    fr = f_r * rhb[i];
    n1[i] /= n1_l2;
    n2[i] /= n2_l2;

    fBase1  =  factor1*rcc[i]  - factor2 * rcb1[i];
    fSugar2 =  factor1*rcb1[i] - factor3 * rcc[i];
    fBase2  = -factor4*rcc[i]  - factor5 * rcb2[i];
    fSugar1 =  factor4*rcb2[i] + factor6 * rcc[i];

    // printf("fBase1 %lf, fSugar2 %lf, fBase2 %lf, fSugar1 %lf\n", fBase1, fSugar2, fBase2, fSugar1);

    force1[i] = dot4*n1[i] - dot3*n2[i] + fSugar1 - fBase1 - fSugar2 - f_sb1 *rcb1[i];
    force2[i] = -fr + dot1*n1[i] + fBase1 + f_sb1 *rcb1[i];
    force3[i] = dot5*n2[i] - dot2*n1[i] + fSugar2 - fBase2 - fSugar1 - f_sb2 *rcb2[i];
    force4[i] = fr - dot1*n2[i] + fBase2 + f_sb2 *rcb2[i];
  }

  // printf("rhb_l %lf thetad %lf f_r %lf, f_d %lf, tau_d = %lf\n", rhb_l, thetad, f_r, f_d, tau_d);

  // PRINT_VECTOR(force1);
  // PRINT_VECTOR(force2);
  // PRINT_VECTOR(force3);
  // PRINT_VECTOR(force4);
 
  return 0;
}


#endif

#endif
