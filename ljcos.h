// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2003; all rights reserved unless otherwise stated.
#include "utils.h"
#include "parser.h"

#ifndef LJCOS_H
#define LJCOS_H
/** \file ljcos.h
 *  Routines to calculate the lennard jones+cosine energy and/or force 
 *  for a particle pair.
 *  \ref forces.c
 *
 *  <b>Responsible:</b>
 *  <a href="mailto:sayar@mpip-mainz.mpg.de">Mehmet</a>
*/

#ifdef LJCOS

MDINLINE int lj_cos_set_params(int part_type_a, int part_type_b,
		      double eps, double sig, double cut,
		      double offset)
{
  IA_parameters *data, *data_sym;

  double facsq;

  make_particle_type_exist(part_type_a);
  make_particle_type_exist(part_type_b);
    
  data     = get_ia_param(part_type_a, part_type_b);
  data_sym = get_ia_param(part_type_b, part_type_a);
  
  if (!data || !data_sym) {
    return TCL_ERROR;
  }

  /* LJCOS should be symmetrically */
  data_sym->LJCOS_eps    = data->LJCOS_eps    = eps;
  data_sym->LJCOS_sig    = data->LJCOS_sig    = sig;
  data_sym->LJCOS_cut    = data->LJCOS_cut    = cut;
  data_sym->LJCOS_offset = data->LJCOS_offset = offset;

  /* Calculate dependent parameters */
  facsq = driwu2*SQR(sig);
  data_sym->LJCOS_rmin = data->LJCOS_rmin = sqrt(driwu2)*sig;
  data_sym->LJCOS_alfa = data->LJCOS_alfa = PI/(SQR(data->LJCOS_cut)-facsq);
  data_sym->LJCOS_beta = data->LJCOS_beta = PI*(1.-(1./(SQR(data->LJCOS_cut)/facsq-1.)));

  /* broadcast interaction parameters */
  mpi_bcast_ia_params(part_type_a, part_type_b);
  mpi_bcast_ia_params(part_type_b, part_type_a);
  
  return TCL_OK;
}

MDINLINE int printljcosIAToResult(Tcl_Interp *interp, int i, int j)
{
  char buffer[TCL_DOUBLE_SPACE + 2*TCL_INTEGER_SPACE];
  IA_parameters *data = get_ia_param(i, j);

    Tcl_PrintDouble(interp, data->LJCOS_eps, buffer);
    Tcl_AppendResult(interp, "lj-cos ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, data->LJCOS_sig, buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, data->LJCOS_cut, buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, data->LJCOS_offset, buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, data->LJCOS_alfa, buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, data->LJCOS_beta, buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, data->LJCOS_rmin, buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);  
    
    return TCL_OK;
}

MDINLINE int ljcos_parser(Tcl_Interp * interp,
			   int part_type_a, int part_type_b,
			   int argc, char ** argv, int *change)
{
  double tmp;
  double eps, sig, cut, offset;
  if (argc < 5) {
    Tcl_AppendResult(interp, "lj-cos needs 4 parameters: "
		     "<ljcos_eps> <ljcos_sig> <ljcos_cut> <ljcos_offset>",
		     (char *) NULL);
    return TCL_ERROR;
  }

  /* copy lj-cos parameters */
  if ((! ARG_IS_D(1, eps))   ||
      (! ARG_IS_D(2, sig))   ||
      (! ARG_IS_D(3, cut))   ||
      (! ARG_IS_D(4, offset)    )) {
    Tcl_AppendResult(interp, "lj-cos needs 4 DOUBLE parameters: "
		     "<ljcos_eps> <ljcos_sig> <ljcos_cut> <ljcos_offset>",
		     (char *) NULL);
    return TCL_ERROR;
  }
  *change = 5;

  /* fix for the inconsistency in the ljcos parameters.
     There are 7 parameters for ljcos, but you read in only four of them.
     The rest is calculated in lj_cos_set_params.
     This is a problem with the blockfile format (Mehmet) 
  */

  if (argc >= 8 && ARG_IS_D(5, tmp) && ARG_IS_D(6, tmp) && ARG_IS_D(7, tmp))
    *change += 3;
  else
    Tcl_ResetResult(interp);

  if (lj_cos_set_params(part_type_a, part_type_b, eps, sig, cut, offset) == TCL_ERROR) {
    Tcl_AppendResult(interp, "particle types must be non-negative", (char *) NULL);
    return TCL_ERROR;
  }

  return TCL_OK;
}

MDINLINE void add_ljcos_pair_force(Particle *p1, Particle *p2, IA_parameters *ia_params,
				double d[3], double dist)
{
  int j;
  double r_off, frac2, frac6, fac=0.0;

  if(dist < ia_params->LJCOS_cut+ia_params->LJCOS_offset) {
    r_off = dist - ia_params->LJCOS_offset;
    /* cos part of ljcos potential. */
    if(dist > ia_params->LJCOS_rmin+ia_params->LJCOS_offset) {
      fac   = (r_off/dist) * ia_params->LJCOS_alfa * ia_params->LJCOS_eps * (sin(ia_params->LJCOS_alfa * SQR(r_off) + ia_params->LJCOS_beta));
      for(j=0;j<3;j++) {
	    /* vector d is rescaled to length LJ_capradius */
	    p1->f.f[j] += fac * d[j];
	    p2->f.f[j] -= fac * d[j];
#ifdef NPT
	    if(integ_switch == INTEG_METHOD_NPT_ISO)
	      nptiso.p_vir += fac*d[j] * d[j];
#endif
      }
    }
    /* lennard-jones part of the potential. */
    else if(dist > 0) {
      frac2 = SQR(ia_params->LJCOS_sig/r_off);
      frac6 = frac2*frac2*frac2;
      fac   = 48.0 * ia_params->LJCOS_eps * frac6*(frac6 - 0.5) / (r_off * dist);

      for(j=0;j<3;j++) {
	    p1->f.f[j] += fac * d[j];
	    p2->f.f[j] -= fac * d[j];
#ifdef NPT
	    if(integ_switch == INTEG_METHOD_NPT_ISO)
	      nptiso.p_vir += fac*d[j] * d[j];
#endif
      }
#ifdef LJ_WARN_WHEN_CLOSE
      if(fac*dist > 1000) fprintf(stderr,"%d: LJCOS-Warning: Pair (%d-%d) force=%f dist=%f\n",
				  this_node,p1->p.identity,p2->p.identity,fac*dist,dist);
#endif
    }
    /* this should not happen! */
    else {
      LJ_TRACE(fprintf(stderr, "%d: Lennard-Jones warning: Particles id1=%d id2=%d exactly on top of each other\n",this_node,p1->p.identity,p2->p.identity));

      frac2 = SQR(ia_params->LJ_sig/ia_params->LJ_capradius);
      frac6 = frac2*frac2*frac2;
      fac   = 48.0 * ia_params->LJ_eps * frac6*(frac6 - 0.5) / ia_params->LJ_capradius;

      p1->f.f[0] += fac * ia_params->LJ_capradius;
      p2->f.f[0] -= fac * ia_params->LJ_capradius;
#ifdef NPT
      if(integ_switch == INTEG_METHOD_NPT_ISO)
	nptiso.p_vir += fac*ia_params->LJ_capradius * ia_params->LJ_capradius;
#endif
    }

    ONEPART_TRACE(if(p1->p.identity==check_id) fprintf(stderr,"%d: OPT: LJ   f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p1->f.f[0],p1->f.f[1],p1->f.f[2],p2->p.identity,dist,fac));
    ONEPART_TRACE(if(p2->p.identity==check_id) fprintf(stderr,"%d: OPT: LJ   f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p2->f.f[0],p2->f.f[1],p2->f.f[2],p1->p.identity,dist,fac));

    LJ_TRACE(fprintf(stderr,"%d: LJ: Pair (%d-%d) dist=%.3f: force+-: (%.3e,%.3e,%.3e)\n",
		     this_node,p1->p.identity,p2->p.identity,dist,fac*d[0],fac*d[1],fac*d[2]));
  }
}


MDINLINE double ljcos_pair_energy(Particle *p1, Particle *p2, IA_parameters *ia_params,
				double d[3], double dist)
{
  double r_off, frac2, frac6;

  if(dist < ia_params->LJCOS_cut+ia_params->LJCOS_offset) {
    r_off = dist-ia_params->LJCOS_offset;
    /* lennard-jones part of the potential. */
    if (dist < (ia_params->LJCOS_rmin+ia_params->LJCOS_offset)) {
      //printf("this is nomal ,  %.3e \n",r_off);
      frac2 = SQR(ia_params->LJCOS_sig/r_off);
      frac6 = frac2*frac2*frac2;
      return 4.0*ia_params->LJCOS_eps*(SQR(frac6)-frac6);
    }
    /* cosine part of the potential. */
    else if (dist < (ia_params->LJCOS_cut+ia_params->LJCOS_offset)) {
      return .5*ia_params->LJCOS_eps*(cos(ia_params->LJCOS_alfa*SQR(r_off)+ia_params->LJCOS_beta)-1.);
    }
    /* this should not happen! */
    else {
      fprintf(stderr,"this is the distance, which is negative %.3e\n",r_off);
    }
  }
  return 0.0;
}

#endif
#endif
