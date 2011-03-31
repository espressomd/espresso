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
#ifndef SUBT_LJ_H
#define SUBT_LJ_H
/** \file subt_lj.h
 *  Routines to subtract the LENNARD-JONES Energy and/or the LENNARD-JONES force 
 *  for a particle pair.
 *  \ref forces.c
*/

#ifdef LENNARD_JONES

/************************************************************/

/// set the parameters for the subtract LJ potential
MDINLINE int subt_lj_set_params(int bond_type, double k, double r)
{
  if(bond_type < 0)
    return TCL_ERROR;

  make_bond_type_exist(bond_type);

  bonded_ia_params[bond_type].p.subt_lj.k = k;
  bonded_ia_params[bond_type].p.subt_lj.r = r;
  bonded_ia_params[bond_type].type = BONDED_IA_SUBT_LJ;  
  bonded_ia_params[bond_type].p.subt_lj.r2 = SQR(bonded_ia_params[bond_type].p.subt_lj.r);
  bonded_ia_params[bond_type].num = 1;

  mpi_bcast_ia_params(bond_type, -1); 

  return TCL_OK;
}

/// parse parameters for the subt_lj potential
MDINLINE int tclcommand_inter_parse_subt_lj(Tcl_Interp *interp, int bond_type, int argc, char **argv)
{
  double k, r;
  if (argc != 3) {
    Tcl_AppendResult(interp, "subt_lj needs 2 dummy parameters: "
		     "<k_subt_lj> <r_subt_lj>", (char *) NULL);
    return TCL_ERROR;
  }

  if ((! ARG_IS_D(1, k)) || (! ARG_IS_D(2, r))) {
    Tcl_AppendResult(interp, "subt_lj needs 2 dummy DOUBLE parameters: "
		     "<k_subt_lj> <r_subt_lj>", (char *) NULL);
    return TCL_ERROR;
  }

  CHECK_VALUE(subt_lj_set_params(bond_type, k, r), "bond type must be nonnegative");
}

/** Computes the negative of the LENNARD-JONES pair forces 
    and adds this force to the particle forces (see \ref tclcommand_inter). 
    @param p1        Pointer to first particle.
    @param p2        Pointer to second/middle particle.
    @param iaparams  Parameters of interaction
    @param dx        change in position
    @param force     force on particles
    @return true if bond is broken
*/
MDINLINE int calc_subt_lj_pair_force(Particle *p1, Particle *p2, Bonded_ia_parameters *iaparams, double dx[3], double force[3])
{
  int i;
  double fac_lj=0.0;
  IA_parameters *ia_params;
  double r_off, frac2, frac6;
  double dist2 = sqrlen(dx);
  double dist = sqrt(dist2);

  if(dist >= iaparams->p.subt_lj.r)
    return 1;

  ia_params = get_ia_param(p1->p.type,p2->p.type);
  if(dist < ia_params->LJ_cut+ia_params->LJ_offset) { 
    r_off = dist - ia_params->LJ_offset;

    if(r_off > ia_params->LJ_capradius) {
      /* normal case: resulting force/energy smaller than capping. */
      frac2 = SQR(ia_params->LJ_sig/r_off);
      frac6 = frac2*frac2*frac2;
      fac_lj = 48.0 * ia_params->LJ_eps * frac6*(frac6 - 0.5) / (r_off * dist);			  
    }
    else if(dist > 0.0) {
      /* capped part of lj potential. */
      frac2 = SQR(ia_params->LJ_sig/ia_params->LJ_capradius);
      frac6 = frac2*frac2*frac2;
      fac_lj = 48.0 * ia_params->LJ_eps * frac6*(frac6 - 0.5) / (ia_params->LJ_capradius * dist);
    }
  } 

  for(i=0;i<3;i++)
    force[i] = -fac_lj*dx[i];

  ONEPART_TRACE(if(p1->p.identity==check_id) fprintf(stderr,"%d: OPT: SUBT_LJ f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac_lj %.3e\n",this_node,p1->f.f[0],p1->f.f[1],p1->f.f[2],p2->p.identity,sqrt(dist2),fac_lj));
  ONEPART_TRACE(if(p2->p.identity==check_id) fprintf(stderr,"%d: OPT: SUBT_LJ f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac_lj %.3e\n",this_node,p2->f.f[0],p2->f.f[1],p2->f.f[2],p1->p.identity,sqrt(dist2),fac_lj));

  return 0;
}

MDINLINE int subt_lj_pair_energy(Particle *p1, Particle *p2, Bonded_ia_parameters *iaparams, double dx[3], double *_energy)
{
  double energy=0.0;
  IA_parameters *ia_params;
  double r_off, frac2, frac6;
  double dist2 = sqrlen(dx);
  double dist = sqrt(dist2);
  
  if(dist >= iaparams->p.subt_lj.r)
    return 1;
  
  ia_params = get_ia_param(p1->p.type,p2->p.type);
  
  if(dist < ia_params->LJ_cut+ia_params->LJ_offset) {
    r_off = dist - ia_params->LJ_offset;
    if(r_off > ia_params->LJ_capradius) {
      /* normal case: resulting force/energy smaller than capping. */
      frac2 = SQR(ia_params->LJ_sig/r_off);
      frac6 = frac2*frac2*frac2;
      energy = 4.0*ia_params->LJ_eps*(SQR(frac6)-frac6+ia_params->LJ_shift);
    }   
    else if(dist > 0.0) {
      /* capped part of lj potential. */
      frac2 = SQR(ia_params->LJ_sig/ia_params->LJ_capradius);
      frac6 = frac2*frac2*frac2;
      energy = 4.0*ia_params->LJ_eps*(SQR(frac6)-frac6+ia_params->LJ_shift);
    }
  }
  *_energy = -energy;
  return 0;
}
#endif
#endif
