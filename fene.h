// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2003; all rights reserved unless otherwise stated.
#ifndef FENE_H
#define FENE_H
/** \file fene.h
 *  Routines to calculate the FENE Energy or/and FENE force 
 *  for a particle pair.
 *  \ref forces.c
 *
 *  <b>Responsible:</b>
 *  <a href="mailto:limbach@mpip-mainz.mpg.de">Hanjo</a>
*/

/************************************************************/

/** Computes the FENE pair force and adds this
    force to the particle forces (see \ref #inter). 
    @param p1        Pointer to first particle.
    @param p2        Pointer to second/middle particle.
    @param type_num  bond type number of the fene interaction (see \ref #inter).
*/
MDINLINE void add_fene_pair_force(Particle *p1, Particle *p2, int type_num)
{
  int i;
  double dx[3], dist2, fac;
  get_mi_vector(dx, p1->r.p, p2->r.p);
  dist2=sqrlen(dx);

  if(dist2 >= bonded_ia_params[type_num].p.fene.r2) {
    fprintf(stderr,"%d: add_fene_pair_force: ERROR: FENE Bond between Pair (%d,%d) broken: dist=%f\n",this_node,
	    p1->r.identity,p2->r.identity,sqrt(dist2)); 
    errexit();
  }
  fac = bonded_ia_params[type_num].p.fene.k;
  fac /= (1.0 - dist2/bonded_ia_params[type_num].p.fene.r2);

  FENE_TRACE(if(fac > 50) fprintf(stderr,"WARNING: FENE force factor between Pair (%d,%d) large: %f at distance %f\n", p1->r.identity,p2->r.identity,fac,sqrt(dist2)) );

  for(i=0;i<3;i++) {
    p1->f[i] -= fac*dx[i];
    p2->f[i] += fac*dx[i];
  }

  ONEPART_TRACE(if(p1->r.identity==check_id) fprintf(stderr,"%d: OPT: FENE f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p1->f[0],p1->f[1],p1->f[2],p2->r.identity,sqrt(dist2),fac));
  ONEPART_TRACE(if(p2->r.identity==check_id) fprintf(stderr,"%d: OPT: FENE f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p2->f[0],p2->f[1],p2->f[2],p1->r.identity,sqrt(dist2),fac));

}

MDINLINE double fene_pair_energy(Particle *p1, Particle *p2, int type_num)
{
  double dx[3], dist2=0.0, energy;
  get_mi_vector(dx, p1->r.p, p2->r.p);
  dist2=sqrlen(dx);
  
  if(dist2 >= bonded_ia_params[type_num].p.fene.r2) {
    fprintf(stderr,"%d: add_fene_pair_force: ERROR: FENE Bond between Pair (%d,%d) broken: dist=%f\n",
	    this_node,p1->r.identity,p2->r.identity,sqrt(dist2)); 
    errexit();
  }

  energy = -0.5*bonded_ia_params[type_num].p.fene.k*bonded_ia_params[type_num].p.fene.r2;
  energy *= log((1.0 - dist2/bonded_ia_params[type_num].p.fene.r2));
  return energy;
}

#endif
