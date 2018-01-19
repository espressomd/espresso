#include"HarmonicDumbbell.hpp"
#include "debug.hpp"
#include "core/random.hpp"

//---HARMONIC DUMBBELL
int Bond::HarmonicDumbbell::calc_bonded_pair_force(Particle *p1, Particle *p2, double dx[3], double force[3]) const {
  #ifdef ROTATION
  int i;
  double fac;
  double dist2 = sqrlen(dx);
  double dist = sqrt(dist2);
  double dr;

  if ((m_r_cut > 0.0) &&
      (dist > m_r_cut)) 
    return 1;

  dr = dist - m_r;
  fac = -m_k1 * dr;
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

  p1->f.torque[0] += m_k2 * da[0];
  p1->f.torque[1] += m_k2 * da[1];
  p1->f.torque[2] += m_k2 * da[2];

  ONEPART_TRACE(if(p1->p.identity==check_id) fprintf(stderr,"%d: OPT: HARMONIC f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p1->f.f[0],p1->f.f[1],p1->f.f[2],p2->p.identity,dist2,fac));
  ONEPART_TRACE(if(p2->p.identity==check_id) fprintf(stderr,"%d: OPT: HARMONIC f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p2->f.f[0],p2->f.f[1],p2->f.f[2],p1->p.identity,dist2,fac));
  #endif
  return 0;
}

int Bond::HarmonicDumbbell::calc_bonded_pair_energy(Particle *p1, Particle *p2, double dx[3], double *_energy) const {

  #ifdef ROTATION
  double dist2 = sqrlen(dx);
  double dist = sqrt(dist2);

  if ((m_r_cut > 0.0) && 
      (dist > m_r_cut)) 
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
  torque[0] = m_k2 * da[0];
  torque[1] = m_k2 * da[1];
  torque[2] = m_k2 * da[2];

  double diff[3];
  diff[0] = dhat[0] - p1->r.quatu[0];
  diff[1] = dhat[1] - p1->r.quatu[1];
  diff[2] = dhat[2] - p1->r.quatu[2];

  *_energy = 0.5*m_k1*SQR(dist - m_r)
           + 0.5*m_k2*(torque[0]*diff[0] + torque[1]*diff[1] + torque[2]*diff[2]);
  #endif
  return 0;

}
