// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2009; all rights reserved unless otherwise stated.
#ifndef STEPPOT_H
#define STEPPOT_H

/** \file steppot.h
 *  Routines to calculate the smooth step potential energy and/or force 
 *  for a particle pair.
 *  \ref forces.c
*/

#ifdef SMOOTH_STEP

MDINLINE int printSmStIAToResult(Tcl_Interp *interp, int i, int j)
{
  char buffer[TCL_DOUBLE_SPACE+TCL_INTEGER_SPACE];
  IA_parameters *data = get_ia_param(i, j);

  Tcl_AppendResult(interp, "smooth-step ", (char *) NULL);
  Tcl_PrintDouble(interp, data->SmSt_d, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  sprintf(buffer, "%d", data->SmSt_n);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->SmSt_eps, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->SmSt_k0, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->SmSt_sig, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->SmSt_cut, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);  
  
  return TCL_OK;
}

MDINLINE int smooth_step_set_params(int part_type_a, int part_type_b,
				      double d, int n, double eps,
				      double k0, double sig,
				      double cut)
{
  IA_parameters *data, *data_sym;

  make_particle_type_exist(part_type_a);
  make_particle_type_exist(part_type_b);
    
  data     = get_ia_param(part_type_a, part_type_b);
  data_sym = get_ia_param(part_type_b, part_type_a);

  if (!data || !data_sym) {
    return TCL_ERROR;
  }

  /* SmSt should be symmetrically */
  data->SmSt_eps    = data_sym->SmSt_eps    = eps;
  data->SmSt_sig    = data_sym->SmSt_sig    = sig;
  data->SmSt_cut    = data_sym->SmSt_cut    = cut;
  data->SmSt_d      = data_sym->SmSt_d      = d;
  data->SmSt_n      = data_sym->SmSt_n      = n;
  data->SmSt_k0     = data_sym->SmSt_k0     = k0;
 
  
  /* broadcast interaction parameters */
  mpi_bcast_ia_params(part_type_a, part_type_b);
  mpi_bcast_ia_params(part_type_b, part_type_a);

  return TCL_OK;
}

MDINLINE int SmSt_parser(Tcl_Interp * interp,
		       int part_type_a, int part_type_b,
		       int argc, char ** argv)
{
  /* parameters needed for LJ */
  double eps, sig, cut, d, k0;
  int n;

  /* get smooth step potential interaction type */
  if (argc < 7) {
    Tcl_AppendResult(interp, "smooth step potential needs 6 parameters: "
		     "<sigma1> <power> <epsilon> <multiplier> <sigma2> <cutoff>",
		     (char *) NULL);
    return 0;
  }

  /* copy smooth step parameters */
  if ((! ARG_IS_D(1, d))     ||
      (! ARG_IS_I(2, n))     ||
      (! ARG_IS_D(3, eps))   ||
      (! ARG_IS_D(4, k0))    ||
      (! ARG_IS_D(5, sig))   ||
      (! ARG_IS_D(6, cut)    )) {
   Tcl_AppendResult(interp, "smooth step potential needs 6 parameters: "
		     "<sigma1> <power> <epsilon> <multiplier> <sigma2> <cutoff>",
		     (char *) NULL);
    return TCL_ERROR;
  }

  if (smooth_step_set_params(part_type_a, part_type_b,
			       d, n, eps, k0, sig,
			       cut) == TCL_ERROR) {
    Tcl_AppendResult(interp, "particle types must be non-negative", (char *) NULL);
    return 0;
  }
  return 7;
}

/** Calculate smooth step force between particle p1 and p2 */
MDINLINE void add_SmSt_pair_force(Particle *p1, Particle *p2, IA_parameters *ia_params,
				  double d[3], double dist,double dist2, double force[3])
{
  int j;
  double frac, fracP, fac=0.0,er;
  if(dist < ia_params->SmSt_cut) {
      frac = ia_params->SmSt_d/dist;
      fracP = pow(frac,ia_params->SmSt_n);
      er=exp(2.*ia_params->SmSt_k0*(dist-ia_params->SmSt_sig));
      fac   =  (ia_params->SmSt_n * fracP+2.*ia_params->SmSt_eps*ia_params->SmSt_k0*dist*er/SQR(1.0+er))/dist2;

      for(j=0;j<3;j++)
	force[j] += fac * d[j];
    
  }
}

/** calculate smooth step potential energy between particle p1 and p2. */
MDINLINE double SmSt_pair_energy(Particle *p1, Particle *p2, IA_parameters *ia_params,
				 double d[3], double dist,double dist2)
{
  double frac, fracP, er;
 
  if(dist < ia_params->SmSt_cut) {
      frac = ia_params->SmSt_d/dist;
      fracP = pow(frac,ia_params->SmSt_n);
      er=exp(2.*ia_params->SmSt_k0*(dist-ia_params->SmSt_sig));
  
      return fracP+ia_params->SmSt_eps/(1.0+er);
    }
  return 0.0;
}

#endif /* ifdef SMOOTH_STEP */
#endif
