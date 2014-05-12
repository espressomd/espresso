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

// extrema for cos(theta), used for the force calculations that involve angles
#define COS_MAX (0.99999999)
#define COS_MIN (-0.99999999)

int cg_dna_basepair_set_params(int bond_type, DoubleList *params);

#define PRINT_VECTOR(A) printf("%lf %lf %lf\n", A[0], A[1], A[2])

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

inline int calc_cg_dna_basepair_force(Particle *p1, Particle *p2, Particle *p3, Particle *p4, Bonded_ia_parameters *iaparams, double force1[3], double force2[3], double force3[3], double force4[3]) {

  /* Base-Base and Sugar-Sugar vectors */
  double rhb[3], rcc[3];
  /* Sugar-Base Vectors */
  double rcb1[3], rcb2[3];
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

  /* helper variables */
  double ra, c, tau_r, tau_d, tau_flip, tau_rd;
  const Bonded_ia_parameters params = *iaparams;
  const double E0 = params.p.cg_dna_basepair.E0;

  get_mi_vector(rhb, p4->r.p, p2->r.p);
  get_mi_vector(rcc, p3->r.p, p1->r.p);
  get_mi_vector(rcb1, p2->r.p, p1->r.p);
  get_mi_vector(rcb2, p4->r.p, p3->r.p);

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

  psi1 = gamma1 >= 1. ? 0. : (gamma1 <= -1. ? M_PI : acos(gamma1));
  psi2 = gamma2 >= 1. ? 0. : (gamma2 <= -1. ? M_PI : acos(gamma2));  

  if((psi1 > COS_MAX) || (psi1 < COS_MIN) || (psi2 > COS_MAX) || (psi2 < COS_MIN)) {
    printf("Warning basepair angle out of bounds.\n");
  }

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

  f_f1 = (gamma1 < 1.) ? -dpsi1 / sqrt(1. - SQR(gamma1)) * params.p.cg_dna_basepair.E01/s1sqr : 1.0;
  f_f2 = (gamma2 < 1.) ? -dpsi2 / sqrt(1. - SQR(gamma2)) * params.p.cg_dna_basepair.E02/s2sqr : 1.0;

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

  const double rcb1_l = norm(rcb1);
  const double rcb2_l = norm(rcb2);
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

    force1[i] = dot4*n1[i] - dot3*n2[i] + fSugar1 - fBase1 - fSugar2;
    force2[i] = -fr + dot1*n1[i] + fBase1;
    force3[i] = dot5*n2[i] - dot2*n1[i] + fSugar2 - fBase2 - fSugar1;
    force4[i] = fr - dot1*n2[i] + fBase2;
  }

  // printf("rhb_l %lf thetad %lf f_r %lf, f_d %lf, tau_d = %lf\n", rhb_l, thetad, f_r, f_d, tau_d);

  // PRINT_VECTOR(force1);
  // PRINT_VECTOR(force2);
  // PRINT_VECTOR(force3);
  // PRINT_VECTOR(force4);
 
  return 0;
}


#endif
