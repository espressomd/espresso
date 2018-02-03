#include"Fene.hpp"
#include "core/random.hpp"
#include "core/errorhandling.hpp"
#include "debug.hpp"

//---FENE---
//calculating the fene bond force: virtual function
int Bond::Fene::calc_bonded_pair_force(Particle *p1, Particle *p2, double dx[3], double force[3]) const {

  int i;
 
  const double len2 = sqrlen(dx);
  const double len = sqrt(len2);
  const double dr = len - m_r0;

  if (dr >= m_drmax){
    return 1;
  };
  
  double fac = -m_k * dr / ((1.0 - dr*dr*m_drmax2i));
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
int Bond::Fene::calc_bonded_pair_energy(Particle *p1, Particle *p2, double dx[3], double *_energy) const {
  /* compute bond stretching (r-r0) */
  double dr = sqrt(sqrlen(dx))-m_r0;

  /* check bond stretching */
  if(dr >= m_drmax) {
    runtimeErrorMsg() <<"FENE bond broken between particles "<< p1->p.identity << " and " << p2->p.identity;
    return 1;
  };

  double energy = -0.5*m_k*m_drmax2;
  energy *= log((1.0 - dr*dr*m_drmax2i));
  *_energy = energy;
  return 0;

}
