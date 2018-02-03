#include "HydrogenBond.hpp"
#include "grid.hpp" // get_mi_vector

/*****************************************************/
// Help Functions
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

/****************************************************/
// class implementation

//force calculation
int Bond::HydrogenBond::calc_bonded_four_particle_force(Particle *p1, Particle *p2, Particle *p3, 
				     Particle *p4, double force[3], double force2[3], 
				     double force3[3], double force4[3]) const {
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
  const double E0 = m_E0;
  
#ifdef CG_DNA_DEBUG
  //  puts("calc_hydrogen_bond_force():");
#endif

  /* Calculate geometry variables */

  get_mi_vector(rcc, p4->r.p, p1->r.p);
  get_mi_vector(rhb, p3->r.p, p2->r.p);
  get_mi_vector(rcb1, p2->r.p, p1->r.p);
  get_mi_vector(rcb2, p3->r.p, p4->r.p);

  rcb1_l = norm(rcb1);
  rcb2_l = norm(rcb2);
  rcc_l = norm(rcc);

  cross(rcc, rcb1, n1);
  cross(rcc, rcb2, n2);

  const double n1_l = norm(n1);
  const double n2_l = norm(n2);

  rhb_l = norm(rhb);

  /* Sugar base interaction */
  
  const double r0sb = m_r0sb;
  const double alphasb = m_alphasb;
  const double f2 = m_f2;
  const double f3 = m_f3;
  const double E0sb = m_E0sb;
  const double c0sb = (1. - 2.*f2)*E0sb*alphasb;
  const double c1sb = (f2-3.*f3)*E0sb*alphasb;
  const double c2sb = f3*E0sb*alphasb;

  ra = (rcb1_l - r0sb)*alphasb;
  f_sb1 = exp(-ra)*ra*(c0sb+c1sb*ra+c2sb*ra*ra)/rcb1_l;

  ra = (rcb2_l - r0sb)*alphasb;
  f_sb2 = exp(-ra)*ra*(c0sb+c1sb*ra+c2sb*ra*ra)/rcb2_l;

  /* Hydrogen bond interaction */

  /* Radial part */

  ra = (rhb_l - m_r0)*m_alpha;
  temp = exp(-ra);
  tau_r = temp *(1.+ra);
  f_r = E0*m_alpha*temp*ra/rhb_l;

  /* Dihedral part */

  gammad = dot(n1, n2)/(n1_l*n2_l);
  tau_d = exp(m_kd*(gammad - 1));
  f_d = -E0*m_kd*tau_d/(n1_l*n2_l);  

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

  dpsi1 = psi1 - m_psi10;
  dpsi2 = psi2 - m_psi20;

  const double sigma1sqr = (Utils::sqr(m_sigma1));
  const double sigma2sqr = (Utils::sqr(m_sigma2));

  f_f1 = -dpsi1 / sqrt(1. - Utils::sqr(gamma1)) * m_E0/sigma1sqr;
  f_f2 = -dpsi2 / sqrt(1. - Utils::sqr(gamma2)) * m_E0/sigma2sqr;
  
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

  double fBase1, fSugar2, fBase2, fSugar1;
  double fr, n1n, n2n;

#ifdef CG_DNA_DEBUG
  int big_force = 0;
#endif

  for(int i = 0; i < 3; i++) {
    fr = f_r * rhb[i];
    n1n = n1[i]/n1_l;
    n2n = n2[i]/n2_l;
   
    fBase1  =  factor1*rcc[i]  - factor2 * rcb1[i];
    fSugar2 =  factor1*rcb1[i] - factor3 * rcc[i];
    fBase2  = -factor4*rcc[i]  - factor5 * rcb2[i];
    fSugar1 =  factor4*rcb2[i] + factor6 * rcc[i];

    force2[i] = -fr + dot1*n1n + fBase1 + f_sb1 *rcb1[i];
    force3[i] =  fr - dot1*n2n + fBase2 + f_sb2 *rcb2[i];

    force[i] = dot4*n1n - dot3*n2n + fSugar1 - fBase1 - fSugar2 - f_sb1 *rcb1[i];
    force4[i] = dot5*n2n - dot2*n1n + fSugar2 - fBase2 - fSugar1 - f_sb2 *rcb2[i];

#ifdef CG_DNA_DEBUG
    if((force2[i] >= 100.) || (force3[i] >= 100.) || (force[i] >= 100.) || (force4[i] >= 100.)) 
      big_force = 1;
#endif
  }

#ifdef CG_DNA_DEBUG
  if(big_force) {  
    puts("Big Force Basepair.");
    PS(p1->p.identity);
    PS(p2->p.identity);
    PS(p3->p.identity);
    PS(p4->p.identity);
    PV(p1->r.p);
    PV(p2->r.p);
    PV(p3->r.p);
    PV(p4->r.p);

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

    PV(force);
    PV(force2);
    PV(force3);
    PV(force4);
  }
#endif

    return 0;
}


//energy calculation
int Bond::HydrogenBond::calc_bonded_four_particle_energy(Particle *p1, Particle *p2, Particle *p3, 
				    Particle *p4, double *_energy) const {

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
  const double E0 = m_E0;
  
#ifdef CG_DNA_DEBUG
  //  puts("calc_hydrogen_bond_force():");
#endif

  /* Calculate geometry variables */

  get_mi_vector(rcc, p4->r.p, p1->r.p);
  get_mi_vector(rhb, p3->r.p, p2->r.p);
  get_mi_vector(rcb1, p2->r.p, p1->r.p);
  get_mi_vector(rcb2, p3->r.p, p4->r.p);

  rcb1_l = norm(rcb1);
  rcb2_l = norm(rcb2);
  rcc_l = norm(rcc);

  cross(rcc, rcb1, n1);
  cross(rcc, rcb2, n2);

  normalize(n1);
  normalize(n2);

  rhb_l = norm(rhb);

  /* Sugar base interaction */
  
  const double r0sb = m_r0sb;
  const double alphasb = m_alphasb;
  const double f2 = m_f2;
  const double f3 = m_f3;
  const double E0sb = m_E0sb;

  ra = (rcb1_l - r0sb)*alphasb;
  potential += E0sb * exp(-ra)*(1.+ra+f2*ra*ra+f3*ra*ra*ra);

  ra = (rcb2_l - r0sb)*alphasb;
  potential += E0sb * exp(-ra)*(1.+ra+f2*ra*ra+f3*ra*ra*ra);
 

  /* Hydrogen bond interaction */

  /* Radial part */

  ra = (rhb_l - m_r0)*m_alpha;
  tau_r = exp(-ra)*(1.+ra);

  /* Dihedral part */

  gammad = dot(n1, n2);
  tau_d = exp(m_kd*(gammad - 1));

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

  dpsi1 = psi1 - m_psi10;
  dpsi2 = psi2 - m_psi20;

  const double sigma1sqr = (Utils::sqr(m_sigma1));
  const double sigma2sqr = (Utils::sqr(m_sigma2));
  
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
