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

#include "hydrogen_bond.hpp"
#include "communication.hpp"

#ifdef CG_DNA_DEBUG
#include "integrate.hpp"
#endif

#ifdef HYDROGEN_BOND

// extrema for cos(theta), used for the force calculations that involve angles
#define COS_MAX (0.99999999)
#define COS_MIN (-0.99999999)

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
  return Utils::sqr(x[0]) + Utils::sqr(x[1]) + Utils::sqr(x[2]);
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

#define EX(A) A[0], A[1], A[2]

#ifdef CG_DNA_DEBUG
#define PV(A) printf(#A " = (%lf %lf %lf)\n",EX(A))
#define PS(A) printf(#A " = %lf\n", (double) A)
#else
#define PV(A) 
#define PS(A)
#endif

int calc_hydrogen_bond_force(Particle *s1, Particle *b1, Particle *b2, Particle *s2, Bonded_ia_parameters *iaparams, double f_s1[3], double f_b1[3], double f_b2[3], double f_s2[3]) {

  /* Base-Base and Sugar-Sugar vectors */
  double rhb[3], rcc[3];
  /* Sugar-Base Vectors */
  double rcb1[3], rcb2[3], rcb1_l, rcb2_l;
  /* Normal vectors of the base pairs */
  double n1[3], n2[3];
  /* Dihedral of the base pair */
  double gammad;
  /* Base angles */
  double psi1, psi2;
  double dpsi1, dpsi2;
  /* gamma1 = cos(psi1) */
  /* length of rhb */
  double gamma1, gamma2;
  double rhb_l, rcc_l;
  /* magnitude of the _r contribution to the force */
  double f_r;
  double f_d;
  double f_f1, f_f2;
  double f_sb1, f_sb2;

  /* helper variables */
  double ra, temp, tau_r, tau_d, tau_flip, tau_rd;
  const Bonded_ia_parameters params = *iaparams;
  const double E0 = params.p.hydrogen_bond.E0;
  
#ifdef CG_DNA_DEBUG
  //  puts("calc_hydrogen_bond_force():");
#endif

  /* Calculate geometry variables */

  get_mi_vector(rcc, s2->r.p, s1->r.p);
  get_mi_vector(rhb, b2->r.p, b1->r.p);
  get_mi_vector(rcb1, b1->r.p, s1->r.p);
  get_mi_vector(rcb2, b2->r.p, s2->r.p);

  rcb1_l = norm(rcb1);
  rcb2_l = norm(rcb2);
  rcc_l = norm(rcc);

  cross(rcc, rcb1, n1);
  cross(rcc, rcb2, n2);

  const double n1_l = norm(n1);
  const double n2_l = norm(n2);

  rhb_l = norm(rhb);

  /* Sugar base interaction */
  
  const double r0sb = params.p.hydrogen_bond.r0sb;
  const double alphasb = params.p.hydrogen_bond.alphasb;
  const double f2 = params.p.hydrogen_bond.f2;
  const double f3 = params.p.hydrogen_bond.f3;
  const double E0sb = params.p.hydrogen_bond.E0sb;
  const double c0sb = (1. - 2.*f2)*E0sb*alphasb;
  const double c1sb = (f2-3.*f3)*E0sb*alphasb;
  const double c2sb = f3*E0sb*alphasb;

  ra = (rcb1_l - r0sb)*alphasb;
  f_sb1 = exp(-ra)*ra*(c0sb+c1sb*ra+c2sb*ra*ra)/rcb1_l;

  ra = (rcb2_l - r0sb)*alphasb;
  f_sb2 = exp(-ra)*ra*(c0sb+c1sb*ra+c2sb*ra*ra)/rcb2_l;

  /* Hydrogen bond interaction */

  /* Radial part */

  ra = (rhb_l - params.p.hydrogen_bond.r0)*params.p.hydrogen_bond.alpha;
  temp = exp(-ra);
  tau_r = temp *(1.+ra);
  f_r = E0*params.p.hydrogen_bond.alpha*temp*ra/rhb_l;

  /* Dihedral part */

  gammad = dot(n1, n2)/(n1_l*n2_l);
  tau_d = exp(params.p.hydrogen_bond.kd*(gammad - 1));
  f_d = -E0*params.p.hydrogen_bond.kd*tau_d/(n1_l*n2_l);  

  /* Flip part */

  gamma1 =  dot(rcc, rcb1)/(rcc_l*rcb1_l);
  gamma2 = -dot(rcc, rcb2)/(rcc_l*rcb2_l);

  /* Avoid illdefined values */

  if((gamma1 > COS_MAX) || (gamma1 < COS_MIN)) {
    gamma1 = (gamma1 > COS_MAX)? COS_MAX : COS_MIN;
  }
  if((gamma2 > COS_MAX) || (gamma2 < COS_MIN))
    gamma2 = (gamma2 > COS_MAX)? COS_MAX : COS_MIN;

  psi1 = gamma1 >= 1. ? 0. : (gamma1 <= -1. ? M_PI : acos(gamma1));
  psi2 = gamma2 >= 1. ? 0. : (gamma2 <= -1. ? M_PI : acos(gamma2));  

  dpsi1 = psi1 - params.p.hydrogen_bond.psi10;
  dpsi2 = psi2 - params.p.hydrogen_bond.psi20;

  const double sigma1sqr = (Utils::sqr(params.p.hydrogen_bond.sigma1));
  const double sigma2sqr = (Utils::sqr(params.p.hydrogen_bond.sigma2));

  f_f1 = -dpsi1 / sqrt(1. - Utils::sqr(gamma1)) * params.p.hydrogen_bond.E0/sigma1sqr;
  f_f2 = -dpsi2 / sqrt(1. - Utils::sqr(gamma2)) * params.p.hydrogen_bond.E0/sigma2sqr;
  
  tau_rd = tau_r * tau_d;

  if(dpsi1 > 0 && dpsi2 > 0) {
    tau_flip = exp(-(Utils::sqr(dpsi1)/(2.*sigma1sqr)+
		     Utils::sqr(dpsi2)/(2.*sigma2sqr)));
    f_f1 *= tau_flip * tau_rd;
    f_f2 *= tau_flip * tau_rd;
  } else if (dpsi1 > 0.) {
    tau_flip = exp(-(Utils::sqr(dpsi1)/(2.*sigma1sqr)));
    f_f1 *= tau_flip * tau_rd;
  } else if (dpsi2 > 0.) {
    tau_flip = exp(-(Utils::sqr(dpsi2)/(2.*sigma2sqr)));
    f_f2 *= tau_flip * tau_rd;
  } else {
    tau_flip = 1.;
  }

  /* Angle at which the constraint sets in */
  const double psi_cutoff = PI*140./180.;
  /* Spring constant for the angle constraint */
  const double k_constraint = 50.;

  if(psi1 > psi_cutoff) {
    f_f1 += k_constraint*(psi1-psi_cutoff)/sqrt(1. - Utils::sqr(gamma1));
  }
  if(psi2 > psi_cutoff) {
    f_f2 += k_constraint*(psi2-psi_cutoff)/sqrt(1. - Utils::sqr(gamma2));
  }

  f_r *= tau_d*tau_flip;  
  f_d *= tau_r*tau_flip;

  /* Dihedral force */
  double vec[3];
  cross(n1, n2, vec);
  
  const double dot1 = f_d * dot(vec, rcc);
  const double dot2 = f_d * dot(vec, rcb1);
  const double dot3 = f_d * dot(vec, rcb2);
  const double dot4 = dot2 - dot1;
  const double dot5 = dot1 + dot3;

  const double factor1 = f_f1/(rcc_l * rcb1_l);
  const double factor2 = f_f1*gamma1/Utils::sqr(rcb1_l);
  const double factor3 = f_f1*gamma1/Utils::sqr(rcc_l);

  const double factor4 = f_f2/(rcc_l * rcb2_l);
  const double factor5 = f_f2*gamma2/Utils::sqr(rcb2_l);
  const double factor6 = f_f2*gamma2/Utils::sqr(rcc_l);


#ifdef CG_DNA_DEBUG
  int big_force = 0;
#endif

  for(int i = 0; i < 3; i++) {
    double fr = f_r * rhb[i];
    double n1n = n1[i]/n1_l;
    double n2n = n2[i]/n2_l;
   
    double fBase1  =  factor1*rcc[i]  - factor2 * rcb1[i];
    double fSugar2 =  factor1*rcb1[i] - factor3 * rcc[i];
    double fBase2  = -factor4*rcc[i]  - factor5 * rcb2[i];
    double fSugar1 =  factor4*rcb2[i] + factor6 * rcc[i];

    f_b1[i] = -fr + dot1*n1n + fBase1 + f_sb1 *rcb1[i];
    f_b2[i] =  fr - dot1*n2n + fBase2 + f_sb2 *rcb2[i];

    f_s1[i] = dot4*n1n - dot3*n2n + fSugar1 - fBase1 - fSugar2 - f_sb1 *rcb1[i];
    f_s2[i] = dot5*n2n - dot2*n1n + fSugar2 - fBase2 - fSugar1 - f_sb2 *rcb2[i];

#ifdef CG_DNA_DEBUG
    if((f_b1[i] >= 100.) || (f_b2[i] >= 100.) || (f_s1[i] >= 100.) || (f_s2[i] >= 100.)) 
      big_force = 1;
#endif
  }

#ifdef CG_DNA_DEBUG
  if(big_force) {  
    puts("Big Force Basepair.");
    PS(s1->p.identity);
    PS(b1->p.identity);
    PS(b2->p.identity);
    PS(s2->p.identity);
    PV(s1->r.p);
    PV(b1->r.p);
    PV(b2->r.p);
    PV(s2->r.p);

    PV(rhb);
    PV(rcc);
    PV(rcb1);
    PV(rcb2);

    PS(rhb_l);
    PS(rcc_l);
    PS(rcb1_l);
    PS(rcb2_l);

    PS(gamma1);
    PS(gamma2);
    PS(dot(n1,n2));
    PV(n1);
    PV(n2);

    PV(dot1*n1);
    PV(dot2*n2);
    
    PV(f_sb1 *rcb1);
    PV(f_sb2 *rcb2);
    PV(f_r*rhb);

    PS(f_r);
    PS(f_r*rhb_l);
    PS(f_sb1);
    PS(f_sb2);

    PS(tau_d);
    PS(tau_r);
    PS(tau_flip);

    PS(dot1);
    PS(dot2);
    PS(dot3);
    PS(dot4);
    PS(dot5);

    PS(factor1);
    PS(factor2);
    PS(factor3);
    PS(factor4);
    PS(factor5);
    PS(factor6);

    PS(fSugar1);
    PS(fBase1);
    PS(fBase2);
    PS(fSugar2);

    PV(f_s1);
    PV(f_b1);
    PV(f_b2);
    PV(f_s2);
  }
#endif
 
  return 0;
}

int calc_hydrogen_bond_energy(Particle *s1, Particle *b1, Particle *b2, Particle *s2, Bonded_ia_parameters *iaparams, double *_energy) {

  /* Base-Base and Sugar-Sugar vectors */
  double rhb[3], rcc[3];
  /* Sugar-Base Vectors */
  double rcb1[3], rcb2[3], rcb1_l, rcb2_l;
  /* Normal vectors of the base pairs */
  double n1[3], n2[3];
  /* Dihedral of the base pair */
  double gammad;
  /* Base angles */
  double psi1, psi2;
  double dpsi1, dpsi2;
  /* gamma1 = cos(psi1) */
  /* length of rhb */
  double gamma1, gamma2;
  double rhb_l, rcc_l;

  double potential = 0.0;

  /* helper variables */
  double ra, tau_r, tau_d, tau_flip;
  const Bonded_ia_parameters params = *iaparams;
  const double E0 = params.p.hydrogen_bond.E0;
  
#ifdef CG_DNA_DEBUG
  //  puts("calc_hydrogen_bond_force():");
#endif

  /* Calculate geometry variables */

  get_mi_vector(rcc, s2->r.p, s1->r.p);
  get_mi_vector(rhb, b2->r.p, b1->r.p);
  get_mi_vector(rcb1, b1->r.p, s1->r.p);
  get_mi_vector(rcb2, b2->r.p, s2->r.p);

  rcb1_l = norm(rcb1);
  rcb2_l = norm(rcb2);
  rcc_l = norm(rcc);

  cross(rcc, rcb1, n1);
  cross(rcc, rcb2, n2);

  normalize(n1);
  normalize(n2);

  rhb_l = norm(rhb);

  /* Sugar base interaction */
  
  const double r0sb = params.p.hydrogen_bond.r0sb;
  const double alphasb = params.p.hydrogen_bond.alphasb;
  const double f2 = params.p.hydrogen_bond.f2;
  const double f3 = params.p.hydrogen_bond.f3;
  const double E0sb = params.p.hydrogen_bond.E0sb;

  ra = (rcb1_l - r0sb)*alphasb;
  potential += E0sb * exp(-ra)*(1.+ra+f2*ra*ra+f3*ra*ra*ra);

  ra = (rcb2_l - r0sb)*alphasb;
  potential += E0sb * exp(-ra)*(1.+ra+f2*ra*ra+f3*ra*ra*ra);
 

  /* Hydrogen bond interaction */

  /* Radial part */

  ra = (rhb_l - params.p.hydrogen_bond.r0)*params.p.hydrogen_bond.alpha;
  tau_r = exp(-ra)*(1.+ra);

  /* Dihedral part */

  gammad = dot(n1, n2);
  tau_d = exp(params.p.hydrogen_bond.kd*(gammad - 1));

  /* Flip part */

  gamma1 =  dot(rcc, rcb1)/(rcc_l*rcb1_l);
  gamma2 = -dot(rcc, rcb2)/(rcc_l*rcb2_l);

  /* Avoid illdefined values */

  if((gamma1 > COS_MAX) || (gamma1 < COS_MIN)) {
    gamma1 = (gamma1 > COS_MAX)? COS_MAX : COS_MIN;
  }
  if((gamma2 > COS_MAX) || (gamma2 < COS_MIN))
    gamma2 = (gamma2 > COS_MAX)? COS_MAX : COS_MIN;

  psi1 = gamma1 >= 1. ? 0. : (gamma1 <= -1. ? M_PI : acos(gamma1));
  psi2 = gamma2 >= 1. ? 0. : (gamma2 <= -1. ? M_PI : acos(gamma2));  

  dpsi1 = psi1 - params.p.hydrogen_bond.psi10;
  dpsi2 = psi2 - params.p.hydrogen_bond.psi20;

  const double sigma1sqr = (Utils::sqr(params.p.hydrogen_bond.sigma1));
  const double sigma2sqr = (Utils::sqr(params.p.hydrogen_bond.sigma2));
  
  if(dpsi1 > 0 && dpsi2 > 0) {
    tau_flip = exp(-(Utils::sqr(dpsi1)/(2.*sigma1sqr)+Utils::sqr(dpsi2)/(2.*sigma2sqr)));
  } else if (dpsi1 > 0.) {
    tau_flip = exp(-(Utils::sqr(dpsi1)/(2.*sigma1sqr)));
    potential += Utils::sqr(dpsi2)*(-0.5*E0/sigma2sqr);
  } else if (dpsi2 > 0.) {
    tau_flip = exp(-(Utils::sqr(dpsi2)/(2.*sigma2sqr)));
    potential += Utils::sqr(dpsi1)*(-0.5*E0/sigma1sqr);
  } else {
    tau_flip = 1.;
    potential += Utils::sqr(dpsi2)*(-0.5*E0/sigma2sqr) + Utils::sqr(dpsi1)*(-0.5*E0/sigma1sqr);
  }

  /* Angle at which the constraint sets in */
  const double psi_cutoff = PI*140./180.;
  /* Spring constant for the angle constraint */
  const double k_constraint = 50.;

  if(psi1 > psi_cutoff) {
    potential += 0.5*k_constraint*Utils::sqr(psi1-psi_cutoff);
  }
  if(psi2 > psi_cutoff) {
    potential += 0.5*k_constraint*Utils::sqr(psi2-psi_cutoff);
  }
 
  potential += tau_r * tau_flip * tau_d * E0;

  *_energy = potential;

  return 0;
}

int hydrogen_bond_set_params(int bond_type, DoubleList *params) {
  if(bond_type < 0)
    return ES_ERROR;

  make_bond_type_exist(bond_type);

  bonded_ia_params[bond_type].p.hydrogen_bond.r0 = params->e[0];
  bonded_ia_params[bond_type].p.hydrogen_bond.alpha = params->e[1];
  bonded_ia_params[bond_type].p.hydrogen_bond.E0 = params->e[2];
  bonded_ia_params[bond_type].p.hydrogen_bond.kd = params->e[3];
  bonded_ia_params[bond_type].p.hydrogen_bond.sigma1 = params->e[4];
  bonded_ia_params[bond_type].p.hydrogen_bond.sigma2 = params->e[5];
  bonded_ia_params[bond_type].p.hydrogen_bond.psi10 = params->e[6];
  bonded_ia_params[bond_type].p.hydrogen_bond.psi20 = params->e[7];
  bonded_ia_params[bond_type].p.hydrogen_bond.E0sb = params->e[8];
  bonded_ia_params[bond_type].p.hydrogen_bond.r0sb = params->e[9];
  bonded_ia_params[bond_type].p.hydrogen_bond.alphasb = params->e[10];
  bonded_ia_params[bond_type].p.hydrogen_bond.f2 = params->e[11];
  bonded_ia_params[bond_type].p.hydrogen_bond.f3 = params->e[12];  

  bonded_ia_params[bond_type].type = BONDED_IA_CG_DNA_BASEPAIR;
  bonded_ia_params[bond_type].num = 3;

  mpi_bcast_ia_params(bond_type, -1);

  return ES_OK;
}

#endif

