// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2003; all rights reserved unless otherwise stated.
#ifndef HARMONIC_H
#define HARMONIC_H
/** \file harmonic.h
 *  Routines to calculate the HARMONIC Energy or/and HARMONIC force 
 *  for a particle pair.
 *  \ref forces.c
 *
 *  <b>Responsible:</b>
 *  <a href="mailto:sayar@mpip-mainz.mpg.de">Mehmet</a>
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
  double dx[3], dist, dist2, fac;
  get_mi_vector(dx, p1->r.p, p2->r.p);
  dist2=sqrlen(dx);
  dist=sqrt(dist2);

  fac = bonded_ia_params[type_num].p.harmonic.k;
  fac *= (dist-bonded_ia_params[type_num].p.harmonic.r);
  fac /= dist;


  for(i=0;i<3;i++) {
    p1->f.f[i] -= fac*dx[i];
    p2->f.f[i] += fac*dx[i];
  }

  ONEPART_TRACE(if(p1->p.identity==check_id) fprintf(stderr,"%d: OPT: HARMONIC f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p1->f.f[0],p1->f.f[1],p1->f.f[2],p2->p.identity,sqrt(dist2),fac));
  ONEPART_TRACE(if(p2->p.identity==check_id) fprintf(stderr,"%d: OPT: HARMONIC f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p2->f.f[0],p2->f.f[1],p2->f.f[2],p1->p.identity,sqrt(dist2),fac));

}

MDINLINE double harmonic_pair_energy(Particle *p1, Particle *p2, int type_num)
{
  double dx[3], dist2=0.0, dist=0.0, energy;
  get_mi_vector(dx, p1->r.p, p2->r.p);
  dist2=sqrlen(dx);
  dist=sqrt(dist2);
  
  energy = 0.5*bonded_ia_params[type_num].p.harmonic.k;
  energy *= SQR(dist-bonded_ia_params[type_num].p.harmonic.r);
  
  return energy;
}

#endif
