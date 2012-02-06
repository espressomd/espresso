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
#ifndef STEPPOT_H
#define STEPPOT_H

/** \file steppot.h
 *  Routines to calculate the smooth step potential energy and/or force 
 *  for a particle pair.
 *  \ref forces.c
*/

#ifdef SMOOTH_STEP

MDINLINE int tclprint_to_result_SmStIA(Tcl_Interp *interp, int i, int j)
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
  IA_parameters *data = get_ia_param_safe(part_type_a, part_type_b);

  if (!data) return TCL_ERROR;

  data->SmSt_eps    = eps;
  data->SmSt_sig    = sig;
  data->SmSt_cut    = cut;
  data->SmSt_d      = d;
  data->SmSt_n      = n;
  data->SmSt_k0     = k0;
 
  /* broadcast interaction parameters */
  mpi_bcast_ia_params(part_type_a, part_type_b);

  return TCL_OK;
}

MDINLINE int tclcommand_inter_parse_SmSt(Tcl_Interp * interp,
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

  if(ia_params->SmSt_cut <=0.) 
   return;
  if(dist > ia_params->SmSt_cut) 
   return;
   
  int j;
  double frac, fracP, fac=0.0,er;
      frac = ia_params->SmSt_d/dist;
      fracP = pow(frac,ia_params->SmSt_n);
      er=exp(2.*ia_params->SmSt_k0*(dist-ia_params->SmSt_sig));
      fac   =  (ia_params->SmSt_n * fracP+2.*ia_params->SmSt_eps*ia_params->SmSt_k0*dist*er/SQR(1.0+er))/dist2;

      for(j=0;j<3;j++)
	force[j] += fac * d[j];
    
  
}

/** calculate smooth step potential energy between particle p1 and p2. */
MDINLINE double SmSt_pair_energy(Particle *p1, Particle *p2, IA_parameters *ia_params,
				 double d[3], double dist,double dist2)
{
  if(ia_params->SmSt_cut <=0.) 
   return 0;
  if(dist > ia_params->SmSt_cut) 
  return 0;

  double frac, fracP, er;
 
      frac = ia_params->SmSt_d/dist;
      fracP = pow(frac,ia_params->SmSt_n);
      er=exp(2.*ia_params->SmSt_k0*(dist-ia_params->SmSt_sig));
  
      return fracP+ia_params->SmSt_eps/(1.0+er);
}

#endif /* ifdef SMOOTH_STEP */
#endif
