/* This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
   It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
   and by which you are legally bound while utilizing this file in any form or way.
   There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
   You should have received a copy of that license along with this program;
   if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
   write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
   Copyright (c) 2002-2009; all rights reserved unless otherwise stated. */
#ifndef FENE_H
#define FENE_H
#include "errorhandling.h"

/** \file fene.h
 *  Routines to calculate the FENE Energy or/and FENE force 
 *  for a particle pair.
 *  \ref forces.c
*/

/************************************************************/

/// set the parameters for the fene potential
MDINLINE int fene_set_params(int bond_type, double k, double drmax, double r0)
{
  if(bond_type < 0)
    return TCL_ERROR;

  make_bond_type_exist(bond_type);

  bonded_ia_params[bond_type].p.fene.k = k;
  bonded_ia_params[bond_type].p.fene.drmax = drmax;
  bonded_ia_params[bond_type].p.fene.r0 = r0;

  bonded_ia_params[bond_type].p.fene.drmax2 = SQR(bonded_ia_params[bond_type].p.fene.drmax);
  bonded_ia_params[bond_type].p.fene.drmax2i = 1.0/bonded_ia_params[bond_type].p.fene.drmax2;

  bonded_ia_params[bond_type].type = BONDED_IA_FENE;
  bonded_ia_params[bond_type].num  = 1;

  /* broadcast interaction parameters */
  mpi_bcast_ia_params(bond_type, -1); 
  
  return TCL_OK;
}

/// parse parameters for the fene potential
MDINLINE int inter_parse_fene(Tcl_Interp *interp, int bond_type, int argc, char **argv)
{
  double k, drmax, r0;

  if (argc != 3 && argc != 4) {
    Tcl_AppendResult(interp, "fene needs 2 or 3 parameters: "
		     "<k> <drmax> [<r0>]", (char *) NULL);
    return TCL_ERROR;
  }

  if ((! ARG_IS_D(1, k)) || (! ARG_IS_D(2, drmax)))
    {
      Tcl_AppendResult(interp, "fene needs 2 or 3 DOUBLE parameters: "
		       "<k> <drmax> [<r0>]", (char *) NULL);
      return TCL_ERROR;
    }

  if (argc == 4) {
    if (! ARG_IS_D(3, r0))
      {
	Tcl_AppendResult(interp, "fene needs 2 or 3 DOUBLE parameters: "
			 "<k> <drmax> [<r0>]", (char *) NULL);
	return TCL_ERROR;
      }
  } else {
    /* default value for r0 is 0.0. */
    r0 = 0.0;
  }
  
  CHECK_VALUE(fene_set_params(bond_type, k, drmax, r0), "bond type must be nonnegative");
}

/** Computes the FENE pair force and adds this
    force to the particle forces (see \ref #inter). 
    @param p1        Pointer to first particle.
    @param p2        Pointer to second/middle particle.
    @param iaparams  bond type number of the angle interaction (see \ref #inter).
    @param dx        particle distance vector
    @param force     returns force of particle 1
    @return true if the bond is broken
*/
MDINLINE int calc_fene_pair_force(Particle *p1, Particle *p2, Bonded_ia_parameters *iaparams, double dx[3], double force[3])
{
  int i;
  double fac, dr, len2, len;
 
  len2 = sqrlen(dx);
  len = sqrt(len2);
  dr = len - iaparams->p.fene.r0;

  if(dr >= iaparams->p.fene.drmax)
    return 1;

  fac = -iaparams->p.fene.k * dr / (len * (1.0 - dr*dr*iaparams->p.fene.drmax2i));
  
  FENE_TRACE(if(fac > 50) fprintf(stderr,"WARNING: FENE force factor between Pair (%d,%d) large: %f at distance %f\n", p1->p.identity,p2->p.identity,fac,sqrt(len2)) );
  
  for(i=0;i<3;i++)
    force[i] = fac*dx[i];

  ONEPART_TRACE(if(p1->p.identity==check_id) fprintf(stderr,"%d: OPT: FENE f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p1->f.f[0],p1->f.f[1],p1->f.f[2],p2->p.identity,sqrt(len2),fac));
  ONEPART_TRACE(if(p2->p.identity==check_id) fprintf(stderr,"%d: OPT: FENE f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p2->f.f[0],p2->f.f[1],p2->f.f[2],p1->p.identity,sqrt(len2),fac));
  
  return 0;
}

MDINLINE int fene_pair_energy(Particle *p1, Particle *p2, Bonded_ia_parameters *iaparams, double dx[3], double *_energy)
{
  double energy, dr;

  /* compute bond stretching (r-r0) */
  dr = sqrt(sqrlen(dx))-iaparams->p.fene.r0;

  /* check bond stretching */
  if(dr >= iaparams->p.fene.drmax) {
    char *errtext = runtime_error(128 + 2*TCL_INTEGER_SPACE);
    ERROR_SPRINTF(errtext,"{077 FENE bond broken between particles %d and %d} ", p1->p.identity, p2->p.identity); 
    return 1;
  }

  energy = -0.5*iaparams->p.fene.k*iaparams->p.fene.drmax2;
  energy *= log((1.0 - dr*dr*iaparams->p.fene.drmax2i));
  *_energy = energy;
  return 0;
}

#endif
