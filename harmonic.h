// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2004; all rights reserved unless otherwise stated.
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

/// set the parameters for the harmonic potential
MDINLINE int harmonic_set_params(int bond_type, double k, double r)
{
  if(bond_type < 0)
    return TCL_ERROR;

  make_bond_type_exist(bond_type);

  bonded_ia_params[bond_type].p.harmonic.k = k;
  bonded_ia_params[bond_type].p.harmonic.r = r;
  bonded_ia_params[bond_type].type = BONDED_IA_HARMONIC;
  bonded_ia_params[bond_type].p.harmonic.r2 = SQR(bonded_ia_params[bond_type].p.harmonic.r);
  bonded_ia_params[bond_type].num  = 1;

  /* broadcast interaction parameters */
  mpi_bcast_ia_params(bond_type, -1); 

  return TCL_OK;
}

/** Computes the HARMONIC pair force and adds this
    force to the particle forces (see \ref #inter). 
    @param p1        Pointer to first particle.
    @param p2        Pointer to second/middle particle.
    @param iaparams  bond type number of the angle interaction (see \ref #inter).
    @param dx        particle distance vector
    @param force     returns force of particle 1
    @return 0.
*/
MDINLINE int calc_harmonic_pair_force(Particle *p1, Particle *p2, Bonded_ia_parameters *iaparams, double dx[3], double force[3])
{
  int i;
  double fac;
  double dist2 = sqrlen(dx);
  double dist = sqrt(dist2);

  fac = -iaparams->p.harmonic.k*(dist - iaparams->p.harmonic.r);
  fac /= dist;

  for(i=0;i<3;i++)
    force[i] = fac*dx[i];

  ONEPART_TRACE(if(p1->p.identity==check_id) fprintf(stderr,"%d: OPT: HARMONIC f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p1->f.f[0],p1->f.f[1],p1->f.f[2],p2->p.identity,dist2,fac));
  ONEPART_TRACE(if(p2->p.identity==check_id) fprintf(stderr,"%d: OPT: HARMONIC f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p2->f.f[0],p2->f.f[1],p2->f.f[2],p1->p.identity,dist2,fac));


  return 0;
}

MDINLINE int harmonic_pair_energy(Particle *p1, Particle *p2, Bonded_ia_parameters *iaparams, double dx[3], double *_energy)
{
  double dist2 = sqrlen(dx);
  double dist = sqrt(dist2);
  *_energy = 0.5*iaparams->p.harmonic.k*SQR(dist - iaparams->p.harmonic.r);
  return 0;
}

#endif
