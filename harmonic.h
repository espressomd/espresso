#ifndef HARMONIC_H
#define HARMONIC_H
/** \file harmonic.h
 *  Routines to calculate the HARMONIC Energy or/and HARMONIC force 
 *  for a particle pair.
 *  \ref forces.c
 *
 *  <b>Responsible:</b>
 *  <a href="mailto:limbach@mpip-mainz.mpg.de">Hanjo</a>
*/

/************************************************************/

/** Computes the HARMONIC pair force and adds this
    force to the particle forces (see \ref #inter). 
    @param p1        Pointer to first particle.
    @param p2        Pointer to second/middle particle.
    @param type_num  bond type number of the harmonic interaction (see \ref #inter).
*/
MDINLINE void add_harmonic_pair_force(Particle *p1, Particle *p2, int type_num)
{
  int i;
  double dx[3], dist=0.0, dist2=0.0, fac;
  for(i=0;i<3;i++) {
    dx[i] = p1->r.p[i] - p2->r.p[i];
    dx[i] -= dround(dx[i]/box_l[i])*box_l[i];
    dist2 += SQR(dx[i]);
  }
  dist=sqrt(dist2);

  fac = bonded_ia_params[type_num].p.harmonic.k;
  fac *= (dist-bonded_ia_params[type_num].p.harmonic.r);
  fac /= dist;


  for(i=0;i<3;i++) {
    p1->f[i] -= fac*dx[i];
    p2->f[i] += fac*dx[i];
  }

  ONEPART_TRACE(if(p1->r.identity==check_id) fprintf(stderr,"%d: OPT: HARMONIC f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p1->f[0],p1->f[1],p1->f[2],p2->r.identity,sqrt(dist2),fac));
  ONEPART_TRACE(if(p2->r.identity==check_id) fprintf(stderr,"%d: OPT: HARMONIC f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p2->f[0],p2->f[1],p2->f[2],p1->r.identity,sqrt(dist2),fac));

}

MDINLINE double harmonic_pair_energy(Particle *p1, Particle *p2, int type_num)
{
  int i;
  double dx[3], dist2=0.0, dist=0.0, energy;
  for(i=0;i<3;i++) {
    dx[i] = p1->r.p[i] - p2->r.p[i];
    dx[i] -= dround(dx[i]/box_l[i])*box_l[i];
    dist2 += SQR(dx[i]);
  }
  dist=sqrt(dist2);
  
  energy = 0.5*bonded_ia_params[type_num].p.harmonic.k;
  energy *= SQR(dist-bonded_ia_params[type_num].p.harmonic.r);
  
  return energy;
}

#endif
