#include"Harmonic.hpp"
#include "debug.hpp"
#include "core/random.hpp"
#include "utils.hpp"
#include "interaction_data.hpp" //Harmonic_bond_parameters

//---HARMONIC BOND---
int Bond::Harmonic::calc_bonded_pair_force(Particle *p1, Particle *p2, double dx[3], double force[3]) const {
  int i;
  double fac;
  double dist2 = sqrlen(dx);
  double dist = sqrt(dist2);
  double dr;

  if ((m_r_cut > 0.0) &&
      (dist > m_r_cut)) 
    return 1;

  dr = dist - m_r;
  fac = -m_k * dr;
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
  return 0;
}

int Bond::Harmonic::calc_bonded_pair_energy(Particle *p1, Particle *p2, double dx[3], double *_energy) const {
  double dist2 = sqrlen(dx);
  double dist = sqrt(dist2);

  if ((m_r_cut > 0.0) && (dist > m_r_cut)){
    runtimeErrorMsg() <<"Harmonic bond broken between particles "<< p1->p.identity <<
      " and " << p2->p.identity;
    return 1;
  };

  *_energy = 0.5*m_k*Utils::sqr(dist - m_r);
  return 0;
}

boost::any Bond::Harmonic::get_bond_parameters_from_bond() const
{

  Harmonic_bond_parameters params = {m_k,  m_r,  m_r_cut};
  return boost::any(params);
  
}
