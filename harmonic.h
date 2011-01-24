/*
  Copyright (C) 2010 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany
  
  This file is part of ESPResSo.
  
  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/
#ifndef HARMONIC_H
#define HARMONIC_H
/** \file harmonic.h
 *  Routines to calculate the HARMONIC Energy or/and HARMONIC force 
 *  for a particle pair.
 *  \ref forces.c
*/

/************************************************************/

/// set the parameters for the harmonic potential
MDINLINE int harmonic_set_params(int bond_type, double k, double r,double r_cut)
{
  if(bond_type < 0)
    return TCL_ERROR;

  make_bond_type_exist(bond_type);

  bonded_ia_params[bond_type].p.harmonic.k = k;
  bonded_ia_params[bond_type].p.harmonic.r = r;
  bonded_ia_params[bond_type].p.harmonic.r_cut = r_cut;
  bonded_ia_params[bond_type].type = BONDED_IA_HARMONIC;
  bonded_ia_params[bond_type].num  = 1;

  /* broadcast interaction parameters */
  mpi_bcast_ia_params(bond_type, -1); 

  return TCL_OK;
}

/// parse parameters for the harmonic potential
MDINLINE int tclcommand_inter_parse_harmonic(Tcl_Interp *interp, int bond_type, int argc, char **argv)
{
  double k, r,r_cut;

  if (argc < 3) {
    Tcl_AppendResult(interp, "harmonic needs at least 2 parameters: "
		     "<k_harmonic> <r_harmonic> [<r_cut>]", (char *) NULL);
    return TCL_ERROR;
  }

  if ((! ARG_IS_D(1, k)) || (! ARG_IS_D(2, r))) {
    Tcl_AppendResult(interp, "harmonic needs at least 2 DOUBLE parameters: "
		     "<k_harmonic> <r_harmonic> [<r_cut>]", (char *) NULL);
    return TCL_ERROR;
  }

  if (argc < 4) {
    r_cut = -1.0;
  } else if (! ARG_IS_D(3, r_cut))  {
    Tcl_AppendResult(interp, "<r_cut> should be DOUBLE", (char *) NULL);
    return TCL_ERROR;
  }

  CHECK_VALUE(harmonic_set_params(bond_type, k, r, r_cut), "bond type must be nonnegative");
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

  if ((iaparams->p.harmonic.r_cut > 0.0) &&
      (dist > iaparams->p.harmonic.r_cut)) 
    return 1;

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

  if ((iaparams->p.harmonic.r_cut > 0.0) && 
      (dist > iaparams->p.harmonic.r_cut)) 
    return 1;

  *_energy = 0.5*iaparams->p.harmonic.k*SQR(dist - iaparams->p.harmonic.r);
  return 0;
}

#endif
