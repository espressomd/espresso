#include"bond.hpp"
//headers from fene.hpp
#include "utils.hpp"
//#include "interaction_data.hpp"
#include "particle_data.hpp"
#include "bonded_interaction.hpp"
#include "random.hpp"
#include "errorhandling.hpp"
#include "debug.hpp"
/*
Member functions and construcor of concrete classes in bond.hpp
 */

//---BOND---
//abstract class
BondedInteraction Bond::get_interaction_code(){
  return interaction_type_code;
}

//---NO BOND---
NO_BOND::NO_BOND(){
  interaction_type_code = BONDED_IA_NONE;
}

int NO_BOND::add_bonded_force(Particle *p1, Particle *p2, double dx[3], double force[3]){
  return 1;
}

int NO_BOND::add_bonded_energy(Particle *p1, Particle *p2, double dx[3], double *_energy){
  return 1;
}


//---FENE---
//constructor needs parameters of fene bond
FENE::FENE(double r0_input, double drmax_input, double drmax2_input, double drmax2i_input, double k_input){
  r0 = r0_input;
  drmax = drmax_input;
  drmax2 = drmax2_input;
  drmax2i = drmax2i_input;
  k = k_input;
  interaction_type_code = BONDED_IA_FENE;
}
//calculating the fene bond force: virtual function
int FENE::add_bonded_force(Particle *p1, Particle *p2, double dx[3], double force[3]){
  int i;
 
  const double len2 = sqrlen(dx);
  const double len = sqrt(len2);
  const double dr = r0;

  if (dr >= drmax) return 1;

  double fac = k * dr / ((1.0 - dr*dr*drmax2i));
  if (fabs(dr) > ROUND_ERROR_PREC) {
     if(len > ROUND_ERROR_PREC) {  /* Regular case */
	fac /= len ; 
     } else { /* dx[] == 0: the force is undefined. Let's use a random direction */
        for(int i = 0;i < 3;i++) dx[i] = d_random()-0.5;
        fac /= sqrt(sqrlen(dx));
     }
  } else { 
    fac = 0.0;
  }
  FENE_TRACE(if(fac > 50) fprintf(stderr,"WARNING: FENE force factor between Pair (%d,%d) large: %f at distance %f\n", p1->p.identity,p2->p.identity,fac,sqrt(len2)) );
  
  for(i=0;i<3;i++)
    force[i] = fac*dx[i];

  ONEPART_TRACE(if(p1->p.identity==check_id) fprintf(stderr,"%d: OPT: FENE f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p1->f.f[0],p1->f.f[1],p1->f.f[2],p2->p.identity,sqrt(len2),fac));
  ONEPART_TRACE(if(p2->p.identity==check_id) fprintf(stderr,"%d: OPT: FENE f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p2->f.f[0],p2->f.f[1],p2->f.f[2],p1->p.identity,sqrt(len2),fac));

  return 0;
}

//calculating FENE Energy: virtual function
int FENE::add_bonded_energy(Particle *p1, Particle *p2, double dx[3], double *_energy){
  /* compute bond stretching (r-r0) */
  double dr = sqrt(sqrlen(dx))-r0;

  /* check bond stretching */
  if(dr >= drmax) {
    runtimeErrorMsg() <<"FENE bond broken between particles "<< p1->p.identity << " and " << p2->p.identity;
    return 1;
  }

  double energy = -0.5*k*drmax2;
  energy *= log((1.0 - dr*dr*drmax2i));
  *_energy = energy;
  return 0;

}

//---HARMONIC DUMBBELL
// consturctor
HARMONIC_DUMBBELL::HARMONIC_DUMBBELL(double k1_i, double k_2_i, double r_i, double r_cut_i){
  k1 = k1_i;
  k2 = k_2_i;
  r = r_i;
  r_cut = r_cut_i;
  interaction_type_code = BONDED_IA_HARMONIC_DUMBBELL;
}

int HARMONIC_DUMBBELL::add_bonded_force(Particle *p1, Particle *p2, double dx[3], double force[3]){
  #ifdef ROTATION
  int i;
  double fac;
  double dist2 = sqrlen(dx);
  double dist = sqrt(dist2);
  double dr;

  if ((r_cut > 0.0) &&
      (dist > r_cut)) 
    return 1;

  dr = dist - r;
  fac = -k1 * dr;
  if (fabs(dr) > ROUND_ERROR_PREC) {
     if (dist > ROUND_ERROR_PREC)  /* Regular case */
        fac /= dist;
     else { /* dx[] == 0: the force is undefined. Let's use a random direction */
        for(i=0;i<3;i++)
	  dx[i] = d_random()-0.5;
        fac /= sqrt(sqrlen(dx));
     }
  } else { 
     fac = 0;
  }
  
  for (int i=0; i<3; i++)
    force[i] = fac*dx[i];

  double dhat[3];
  dhat[0] = dx[0]/dist;
  dhat[1] = dx[1]/dist;
  dhat[2] = dx[2]/dist;

  double da[3];
  da[0] = dhat[1]*p1->r.quatu[2] - dhat[2]*p1->r.quatu[1];
  da[1] = dhat[2]*p1->r.quatu[0] - dhat[0]*p1->r.quatu[2];
  da[2] = dhat[0]*p1->r.quatu[1] - dhat[1]*p1->r.quatu[0];

  p1->f.torque[0] += k2 * da[0];
  p1->f.torque[1] += k2 * da[1];
  p1->f.torque[2] += k2 * da[2];

  ONEPART_TRACE(if(p1->p.identity==check_id) fprintf(stderr,"%d: OPT: HARMONIC f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p1->f.f[0],p1->f.f[1],p1->f.f[2],p2->p.identity,dist2,fac));
  ONEPART_TRACE(if(p2->p.identity==check_id) fprintf(stderr,"%d: OPT: HARMONIC f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p2->f.f[0],p2->f.f[1],p2->f.f[2],p1->p.identity,dist2,fac));
  #endif
  return 0;
}

int HARMONIC_DUMBBELL::add_bonded_energy(Particle *p1, Particle *p2, double dx[3], double *_energy){

  #ifdef ROTATION
  double dist2 = sqrlen(dx);
  double dist = sqrt(dist2);

  if ((r_cut > 0.0) && 
      (dist > r_cut)) 
    return 1;

  double dhat[3];
  dhat[0] = dx[0]/dist;
  dhat[1] = dx[1]/dist;
  dhat[2] = dx[2]/dist;

  double da[3];
  da[0] = dhat[1]*p1->r.quatu[2] - dhat[2]*p1->r.quatu[1];
  da[1] = dhat[2]*p1->r.quatu[0] - dhat[0]*p1->r.quatu[2];
  da[2] = dhat[0]*p1->r.quatu[1] - dhat[1]*p1->r.quatu[0];

  double torque[3];
  torque[0] = k2 * da[0];
  torque[1] = k2 * da[1];
  torque[2] = k2 * da[2];

  double diff[3];
  diff[0] = dhat[0] - p1->r.quatu[0];
  diff[1] = dhat[1] - p1->r.quatu[1];
  diff[2] = dhat[2] - p1->r.quatu[2];

  *_energy = 0.5*k1*SQR(dist - r)
           + 0.5*k2*(torque[0]*diff[0] + torque[1]*diff[1] + torque[2]*diff[2]);
  #endif
  return 0;

}

//---HARMONIC BOND---
// constuctor
HARMONIC::HARMONIC(double k_i, double r_i, double r_cut_i){
  k = k_i;
  r = r_i;
  r_cut = r_cut_i;
  interaction_type_code = BONDED_IA_HARMONIC;
} 

int HARMONIC::add_bonded_force(Particle *p1, Particle *p2, double dx[3], double force[3]){
  int i;
  double fac;
  double dist2 = sqrlen(dx);
  double dist = sqrt(dist2);
  double dr;

  if ((r_cut > 0.0) &&
      (dist > r_cut)) 
    return 1;

  dr = dist - r;
  fac = -k * dr;
  if (fabs(dr) > ROUND_ERROR_PREC) {
     if(dist>ROUND_ERROR_PREC) {  /* Regular case */
        fac /= dist;
     } else { /* dx[] == 0: the force is undefined. Let's use a random direction */
        for(i=0;i<3;i++) dx[i] = d_random()-0.5;
        fac /= sqrt(sqrlen(dx));
     }
  } else { 
     fac=0;
  }
  
  for(i=0;i<3;i++)
    force[i] = fac*dx[i];
  ONEPART_TRACE(if(p1->p.identity==check_id) fprintf(stderr,"%d: OPT: HARMONIC f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p1->f.f[0],p1->f.f[1],p1->f.f[2],p2->p.identity,dist2,fac));
  ONEPART_TRACE(if(p2->p.identity==check_id) fprintf(stderr,"%d: OPT: HARMONIC f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p2->f.f[0],p2->f.f[1],p2->f.f[2],p1->p.identity,dist2,fac));

#ifdef CONFIGTEMP
  extern double configtemp[2];
  int numfac = 0;
  if (p1->p.configtemp) numfac+=1;
  if (p2->p.configtemp) numfac+=1;
  configtemp[0] += numfac*SQR(k * dr);
  configtemp[1] -= numfac*k*(3-2.*r/dist);
#endif

  return 0;
}

int HARMONIC::add_bonded_energy(Particle *p1, Particle *p2, double dx[3], double *_energy){
  double dist2 = sqrlen(dx);
  double dist = sqrt(dist2);

  if ((r_cut > 0.0) && 
      (dist > r_cut)) 
    return 1;

  *_energy = 0.5*k*SQR(dist - r);
  return 0;
}

//---QUARTIC BOND---
//constructor
QUARTIC::QUARTIC(double k0_i, double k_1_i, double r_i, double r_cut_i){
  k0 = k0_i;
  k1 = k_1_i;
  r = r_i;
  r_cut = r_cut_i;
  interaction_type_code = BONDED_IA_QUARTIC;
}

int QUARTIC::add_bonded_force(Particle *p1, Particle *p2, double dx[3], double force[3]){
  int i;
  double fac;
  double dist2 = sqrlen(dx);
  double dist = sqrt(dist2);
  double dr;

  // printf("Quartic dist2 %e, dist %e\n", dist2, dist);

  if ((r_cut > 0.0) &&
      (dist > r_cut)) 
    return 1;

  dr = dist - r;

  fac = (k0 * dr + k1 * dr * dr * dr)/dist;
  
  for(i=0;i<3;i++)
    force[i] = -fac*dx[i];

  //  printf("Quartic (%d-%d), dist %e, dx %e %e %e, dr %e, f %e %e %e\n", p1->p.identity, p2->p.identity, dist, dx[0], dx[1], dx[2], dr, force[0], force[1], force[2]);

  ONEPART_TRACE(if(p1->p.identity==check_id) fprintf(stderr,"%d: OPT: QUARTIC f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p1->f.f[0],p1->f.f[1],p1->f.f[2],p2->p.identity,dist2,fac));
  ONEPART_TRACE(if(p2->p.identity==check_id) fprintf(stderr,"%d: OPT: QUARTIC f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p2->f.f[0],p2->f.f[1],p2->f.f[2],p1->p.identity,dist2,fac));

  return 0;
}

int QUARTIC::add_bonded_energy(Particle *p1, Particle *p2, double dx[3], double *_energy){
  double dist2 = sqrlen(dx);
  double dist = sqrt(dist2);

  if ((r_cut > 0.0) && 
      (dist > r_cut)) 
    return 1;
 
  double dr2 = SQR(dist -r);

  *_energy = 0.5*k0*dr2 + 0.25 * k1 * SQR(dr2);
  return 0;
}

//---BONDED_COULOMB---
//constructor
BONDED_COULOMB::BONDED_COULOMB(double prefactor_i){
  prefactor = prefactor_i;
  interaction_type_code = BONDED_IA_BONDED_COULOMB;
}

int BONDED_COULOMB::add_bonded_force(Particle *p1, Particle *p2, double dx[3], double force[3]){
#ifdef ELECTROSTATICS
  int i;
  double fac;
  double dist2 = sqrlen(dx);
  double dist = sqrt(dist2);

  fac = prefactor * p1->p.q * p2->p.q / (dist*dist2);

  for(i=0;i<3;i++)
    force[i] = fac*dx[i];
  ONEPART_TRACE(if(p1->p.identity==check_id) fprintf(stderr,"%d: OPT: BONDED_COULOMB f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p1->f.f[0],p1->f.f[1],p1->f.f[2],p2->p.identity,dist2,fac));
  ONEPART_TRACE(if(p2->p.identity==check_id) fprintf(stderr,"%d: OPT: BONDED_COULOMB f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p2->f.f[0],p2->f.f[1],p2->f.f[2],p1->p.identity,dist2,fac));
#endif
  return 0;
}

int BONDED_COULOMB::add_bonded_energy(Particle *p1, Particle *p2, double dx[3], double *_energy){
#ifdef ELECTROSTATICS
  double dist = sqrt(sqrlen(dx));
  *_energy = prefactor * p1->p.q * p2->p.q / dist;
#endif
  return 0;
}
